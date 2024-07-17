#include <gtest/gtest.h>

#include "scran_tests/scran_tests.hpp"

#include "tatami/tatami.hpp"
#include "tatami_stats/tatami_stats.hpp"
#include "scran_variances/model_gene_variances.hpp"

#include <cmath>

class ModelGeneVariancesTest : public ::testing::TestWithParam<int> {
protected:
    inline static std::shared_ptr<tatami::NumericMatrix> dense_row, dense_column, sparse_row, sparse_column;

    static void SetUpTestSuite() {
        int nr = 178, nc = 155;
        auto vec = scran_tests::simulate_vector(nr * nc, []{
            scran_tests::SimulationParameters sparams;
            sparams.density = 0.3;
            sparams.lower = 0;
            sparams.upper = 5;
            return sparams;
        }());

        dense_row = std::unique_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double, int>(nr, nc, std::move(vec)));
        dense_column = tatami::convert_to_dense(dense_row.get(), false);
        sparse_row = tatami::convert_to_compressed_sparse(dense_row.get(), true);
        sparse_column = tatami::convert_to_compressed_sparse(dense_row.get(), false);
    }
};

TEST_P(ModelGeneVariancesTest, UnblockedStats) {
    scran_variances::ModelGeneVariancesOptions opt;
    auto ref = scran_variances::model_gene_variances(*dense_row, opt);

    auto nthreads = GetParam();
    opt.num_threads = nthreads;

    if (nthreads == 1) {
        // Cursory checks.
        auto means = tatami_stats::sums::by_row(dense_row.get());
        for (auto& x : means) {
            x /= dense_row->ncol();
        }
        scran_tests::compare_almost_equal(ref.means, means);
        scran_tests::compare_almost_equal(ref.variances, tatami_stats::variances::by_row(dense_row.get()));

        for (auto f : ref.fitted) {
            EXPECT_TRUE(f > 0);
        }

        int nonzero = 0; 
        for (auto f : ref.residuals) {
            nonzero += (f != 0);
        }
        EXPECT_TRUE(nonzero > 0); // there is at least one non-zero residual; but we can't expect this of everyone.

    } else {
        // Checking against the same call, but parallelized.
        auto res1 = scran_variances::model_gene_variances(*dense_row, opt);
        EXPECT_EQ(ref.means, res1.means);
        EXPECT_EQ(ref.variances, res1.variances);
    }

    // Almost equal, as there are minor differences due to numerical imprecision.
    auto res2 = scran_variances::model_gene_variances(*dense_column, opt);
    scran_tests::compare_almost_equal(ref.means, res2.means);
    scran_tests::compare_almost_equal(ref.variances, res2.variances);

    auto res3 = scran_variances::model_gene_variances(*sparse_row, opt);
    scran_tests::compare_almost_equal(ref.means, res3.means);
    scran_tests::compare_almost_equal(ref.variances, res3.variances);

    auto res4 = scran_variances::model_gene_variances(*sparse_column, opt);
    scran_tests::compare_almost_equal(ref.means, res4.means);
    scran_tests::compare_almost_equal(ref.variances, res4.variances);
}

TEST_P(ModelGeneVariancesTest, BlockedStats) {
    std::vector<int> blocks(dense_row->ncol());
    for (size_t i = 0; i < blocks.size(); ++i) {
        blocks[i] = i % 3;
    }

    scran_variances::ModelGeneVariancesOptions opt;
    opt.compute_average = false;
    auto ref = scran_variances::model_gene_variances_blocked(*dense_row, blocks.data(), opt);
    EXPECT_TRUE(ref.average.means.empty());
    EXPECT_TRUE(ref.average.variances.empty());

    auto nthreads = GetParam();
    opt.num_threads = nthreads;

    if (nthreads == 1) {
        // Cursory checks.
        EXPECT_EQ(ref.per_block.size(), 3);
    } else {
        // Checking against the same call, but parallelized.
        auto res1 = scran_variances::model_gene_variances_blocked(*dense_row, blocks.data(), opt);
        for (size_t i = 0; i < 3; ++i) {
            EXPECT_EQ(ref.per_block[i].means, res1.per_block[i].means);
            EXPECT_EQ(ref.per_block[i].variances, res1.per_block[i].variances);
            EXPECT_EQ(ref.per_block[i].fitted, res1.per_block[i].fitted);
            EXPECT_EQ(ref.per_block[i].residuals, res1.per_block[i].residuals);
        }
    }

    auto res2 = scran_variances::model_gene_variances_blocked(*dense_column, blocks.data(), opt);
    for (size_t i = 0; i < 3; ++i) {
        scran_tests::compare_almost_equal(ref.per_block[i].means, res2.per_block[i].means);
        scran_tests::compare_almost_equal(ref.per_block[i].variances, res2.per_block[i].variances);
    }

    auto res3 = scran_variances::model_gene_variances_blocked(*sparse_row, blocks.data(), opt);
    for (size_t i = 0; i < 3; ++i) {
        scran_tests::compare_almost_equal(ref.per_block[i].means, res3.per_block[i].means);
        scran_tests::compare_almost_equal(ref.per_block[i].variances, res3.per_block[i].variances);
    }

    auto res4 = scran_variances::model_gene_variances_blocked(*sparse_column, blocks.data(), opt);
    for (size_t i = 0; i < 3; ++i) {
        scran_tests::compare_almost_equal(ref.per_block[i].means, res3.per_block[i].means);
        scran_tests::compare_almost_equal(ref.per_block[i].variances, res3.per_block[i].variances);
    }

    // Checking averages with equiweighting.
    opt.compute_average = true;
    opt.block_weight_policy = scran_blocks::WeightPolicy::EQUAL;
    {
        auto ares = scran_variances::model_gene_variances_blocked(*dense_row, blocks.data(), opt);
        EXPECT_EQ(ares.per_block[0].means, ref.per_block[0].means);
        EXPECT_EQ(ares.per_block[1].variances, ref.per_block[1].variances);
        EXPECT_EQ(ares.per_block[2].fitted, ref.per_block[2].fitted);
        EXPECT_EQ(ares.per_block[1].residuals, ref.per_block[1].residuals);

        std::vector<double> expected_means(dense_row->nrow()), 
            expected_variances(dense_row->nrow()),
            expected_fitted(dense_row->nrow()),
            expected_residuals(dense_row->nrow());

        for (size_t r = 0, rend = expected_means.size(); r < rend; ++r) {
            for (int b = 0; b < 3; ++b) {
                expected_means[r] += ares.per_block[b].means[r];
                expected_variances[r] += ares.per_block[b].variances[r];
                expected_fitted[r] += ares.per_block[b].fitted[r];
                expected_residuals[r] += ares.per_block[b].residuals[r];
            }

            expected_means[r] /= 3;
            expected_variances[r] /= 3;
            expected_fitted[r] /= 3;
            expected_residuals[r] /= 3;
        }

        scran_tests::compare_almost_equal(expected_means, ares.average.means);
        scran_tests::compare_almost_equal(expected_variances, ares.average.variances);
        scran_tests::compare_almost_equal(expected_fitted, ares.average.fitted);
        scran_tests::compare_almost_equal(expected_residuals, ares.average.residuals);

        // Checking limit of the variable policy.
        opt.block_weight_policy = scran_blocks::WeightPolicy::VARIABLE;
        opt.variable_block_weight_parameters.lower_bound = 0;
        opt.variable_block_weight_parameters.upper_bound = 0;

        auto vres = scran_variances::model_gene_variances_blocked(*dense_row, blocks.data(), opt);
        scran_tests::compare_almost_equal(ares.average.means, vres.average.means);
        scran_tests::compare_almost_equal(ares.average.variances, vres.average.variances);
        scran_tests::compare_almost_equal(ares.average.fitted, vres.average.fitted);
        scran_tests::compare_almost_equal(ares.average.residuals, vres.average.residuals);
    }

    // Checking averages without equiweighting.
    opt.block_weight_policy = scran_blocks::WeightPolicy::NONE;
    {
        auto ares = scran_variances::model_gene_variances_blocked(*dense_row, blocks.data(), opt);
        auto block_size = tatami_stats::tabulate_groups(blocks.data(), blocks.size());

        std::vector<double> expected_means(dense_row->nrow()), 
            expected_variances(dense_row->nrow()),
            expected_fitted(dense_row->nrow()),
            expected_residuals(dense_row->nrow());

        for (size_t r = 0, rend = expected_means.size(); r < rend; ++r) {
            for (int b = 0; b < 3; ++b) {
                expected_means[r] += ares.per_block[b].means[r] * block_size[b];
                expected_variances[r] += ares.per_block[b].variances[r] * block_size[b];
                expected_fitted[r] += ares.per_block[b].fitted[r] * block_size[b];
                expected_residuals[r] += ares.per_block[b].residuals[r] * block_size[b];
            }

            expected_means[r] /= blocks.size();
            expected_variances[r] /= blocks.size();
            expected_fitted[r] /= blocks.size();
            expected_residuals[r] /= blocks.size();
        }

        scran_tests::compare_almost_equal(expected_means, ares.average.means);
        scran_tests::compare_almost_equal(expected_variances, ares.average.variances);
        scran_tests::compare_almost_equal(expected_fitted, ares.average.fitted);
        scran_tests::compare_almost_equal(expected_residuals, ares.average.residuals);

        // Checking limit of the variable policy.
        opt.block_weight_policy = scran_blocks::WeightPolicy::VARIABLE;
        opt.variable_block_weight_parameters.lower_bound = 0;
        opt.variable_block_weight_parameters.upper_bound = 100000;

        auto vres = scran_variances::model_gene_variances_blocked(*dense_row, blocks.data(), opt);
        scran_tests::compare_almost_equal(ares.average.means, vres.average.means);
        scran_tests::compare_almost_equal(ares.average.variances, vres.average.variances);
        scran_tests::compare_almost_equal(ares.average.fitted, vres.average.fitted);
        scran_tests::compare_almost_equal(ares.average.residuals, vres.average.residuals);
    }
}

INSTANTIATE_TEST_SUITE_P(
    ModelGeneVariances,
    ModelGeneVariancesTest,
    ::testing::Values(1, 3) // number of threads
);

TEST(ModelGeneVariances, EmptyBlocks) {
    int nr = 99, nc = 201;
    auto vec = scran_tests::simulate_vector(nr * nc, []{
        scran_tests::SimulationParameters sparams;
        sparams.density = 0.1;
        sparams.seed = 69;
        return sparams;
    }());
    auto dense_row = std::unique_ptr<tatami::NumericMatrix>(new tatami::DenseRowMatrix<double, int>(nr, nc, std::move(vec)));

    std::vector<int> blocks(dense_row->ncol());
    blocks[0] = 1;
    blocks[1] = 3;

    scran_variances::ModelGeneVariancesOptions opt;
    opt.block_weight_policy = scran_blocks::WeightPolicy::NONE;

    auto res = scran_variances::model_gene_variances_blocked(*dense_row, blocks.data(), opt);
    EXPECT_FALSE(std::isnan(res.per_block[0].means[0]));
    EXPECT_FALSE(std::isnan(res.per_block[0].variances[0]));
    EXPECT_FALSE(std::isnan(res.per_block[1].means[0]));
    EXPECT_TRUE(std::isnan(res.per_block[1].variances[0]));
    EXPECT_TRUE(std::isnan(res.per_block[2].means[0]));
    EXPECT_TRUE(std::isnan(res.per_block[2].variances[0]));
    EXPECT_FALSE(std::isnan(res.per_block[3].means[0]));
    EXPECT_TRUE(std::isnan(res.per_block[3].variances[0]));

    auto expected = tatami_stats::sums::by_row(dense_row.get());
    for (auto& x : expected) {
        x /= nc;
    }
    scran_tests::compare_almost_equal(res.average.means, expected);
    scran_tests::compare_almost_equal(res.average.variances, res.per_block[0].variances);
}
