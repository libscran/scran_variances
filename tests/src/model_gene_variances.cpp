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

TEST_P(ModelGeneVariancesTest, Unblocked) {
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

TEST_P(ModelGeneVariancesTest, Blocked) {
    std::vector<int> blocks(dense_row->ncol());
    for (size_t i = 0; i < blocks.size(); ++i) {
        blocks[i] = i % 3;
    }

    scran_variances::ModelGeneVariancesOptions opt;
    opt.average_policy = scran_variances::AveragePolicy::NONE;
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
}

TEST_P(ModelGeneVariancesTest, BlockedMean) {
    std::vector<int> blocks(dense_row->ncol());
    for (size_t i = 0; i < blocks.size(); ++i) {
        blocks[i] = i % 3;
    }

    scran_variances::ModelGeneVariancesOptions opt;
    opt.average_policy = scran_variances::AveragePolicy::MEAN;

    // Checking averages with equiweighting.
    opt.block_weight_policy = scran_blocks::WeightPolicy::EQUAL;
    {
        auto ares = scran_variances::model_gene_variances_blocked(*dense_row, blocks.data(), opt);
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
    opt.block_weight_policy = scran_blocks::WeightPolicy::SIZE;
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

TEST_P(ModelGeneVariancesTest, BlockedMedian) {
    std::vector<int> blocks(dense_row->ncol());
    for (size_t i = 0; i < blocks.size(); ++i) {
        blocks[i] = i % 3;
    }

    scran_variances::ModelGeneVariancesOptions opt;
    opt.average_policy = scran_variances::AveragePolicy::QUANTILE;
    auto ares = scran_variances::model_gene_variances_blocked(*dense_row, blocks.data(), opt);

    std::vector<double> expected_means(dense_row->nrow()), buffer_means,
        expected_variances(dense_row->nrow()), buffer_variances,
        expected_fitted(dense_row->nrow()), buffer_fitted,
        expected_residuals(dense_row->nrow()), buffer_residuals;

    for (size_t r = 0, rend = expected_means.size(); r < rend; ++r) {
        buffer_means.clear();
        buffer_variances.clear();
        buffer_fitted.clear();
        buffer_residuals.clear();

        for (int b = 0; b < 3; ++b) {
            buffer_means.push_back(ares.per_block[b].means[r]);
            buffer_variances.push_back(ares.per_block[b].variances[r]);
            buffer_fitted.push_back(ares.per_block[b].fitted[r]);
            buffer_residuals.push_back(ares.per_block[b].residuals[r]);
        }

        std::sort(buffer_means.begin(), buffer_means.end());
        std::sort(buffer_variances.begin(), buffer_variances.end());
        std::sort(buffer_fitted.begin(), buffer_fitted.end());
        std::sort(buffer_residuals.begin(), buffer_residuals.end());

        expected_means[r] = buffer_means[1];
        expected_variances[r] = buffer_variances[1];
        expected_fitted[r] = buffer_fitted[1];
        expected_residuals[r] = buffer_residuals[1];
    }

    EXPECT_EQ(expected_means, ares.average.means);
    EXPECT_EQ(expected_variances, ares.average.variances);
    EXPECT_EQ(expected_fitted, ares.average.fitted);
    EXPECT_EQ(expected_residuals, ares.average.residuals);
}

INSTANTIATE_TEST_SUITE_P(
    ModelGeneVariances,
    ModelGeneVariancesTest,
    ::testing::Values(1, 3) // number of threads
);

class ModelGeneVariancesNearEmptyBlockTest : public ::testing::Test {
protected:
    int nr = 109, nc = 152;
    std::shared_ptr<tatami::NumericMatrix> dense_row;
    std::vector<int> blocks;

    void SetUp() {
        auto vec = scran_tests::simulate_vector(nr * nc, []{
            scran_tests::SimulationParameters sparams;
            sparams.density = 0.1;
            sparams.seed = 69;
            sparams.lower = 0;
            sparams.upper = 5;
            return sparams;
        }());
        dense_row.reset(new tatami::DenseMatrix<double, int, decltype(vec)>(nr, nc, std::move(vec), true));

        // Block 0 has everything in it.
        // Block 1 and 3 have one observation.
        // Block 4 has two observations.
        blocks.resize(nc);
        blocks[0] = 1;
        blocks[1] = 3;
        blocks[2] = 4;
        blocks[3] = 4;
    }
};

TEST_F(ModelGeneVariancesNearEmptyBlockTest, Mean) {
    scran_variances::ModelGeneVariancesOptions opt;
    opt.block_weight_policy = scran_blocks::WeightPolicy::SIZE;

    auto res = scran_variances::model_gene_variances_blocked(*dense_row, blocks.data(), opt);
    EXPECT_FALSE(std::isnan(res.per_block[0].means[0]));
    EXPECT_FALSE(std::isnan(res.per_block[0].variances[0]));
    EXPECT_FALSE(std::isnan(res.per_block[1].means[0]));
    EXPECT_TRUE(std::isnan(res.per_block[1].variances[0]));
    EXPECT_TRUE(std::isnan(res.per_block[2].means[0]));
    EXPECT_TRUE(std::isnan(res.per_block[2].variances[0]));
    EXPECT_FALSE(std::isnan(res.per_block[3].means[0]));
    EXPECT_TRUE(std::isnan(res.per_block[3].variances[0]));
    EXPECT_FALSE(std::isnan(res.per_block[4].means[0]));
    EXPECT_FALSE(std::isnan(res.per_block[4].variances[0]));

    // Subset to blocks 0 and 4 to ignore all blocks with fewer than two cells. This means we do [2, ncol).
    tatami::DelayedSubsetBlock<double, int> sub(dense_row, 2, nc - 2, false);
    std::vector<int> subblocks(nc - 2);
    subblocks[0] = 1;
    subblocks[1] = 1;
    auto ref = scran_variances::model_gene_variances_blocked(sub, subblocks.data(), opt);
    EXPECT_EQ(ref.average.variances, res.average.variances);
    EXPECT_EQ(ref.average.fitted, res.average.fitted);
    EXPECT_EQ(ref.average.residuals, res.average.residuals);

    // Computing the expected mean. As the weight policy is SIZE, we basically just take the row means.
    auto expected = tatami_stats::sums::by_row(dense_row.get());
    for (auto& x : expected) {
        x /= nc;
    }
    scran_tests::compare_almost_equal_containers(expected, res.average.means, {});
}

TEST_F(ModelGeneVariancesNearEmptyBlockTest, Median) {
    scran_variances::ModelGeneVariancesOptions opt;
    opt.average_policy = scran_variances::AveragePolicy::QUANTILE;
    auto res = scran_variances::model_gene_variances_blocked(*dense_row, blocks.data(), opt);

    // Subset to blocks 0 and 4 to ignore all blocks with fewer than two cells. This means we do [2, ncol).
    tatami::DelayedSubsetBlock<double, int> sub(dense_row, 2, nc - 2, false);
    std::vector<int> subblocks(nc - 2);
    subblocks[0] = 1;
    subblocks[1] = 1;

    auto ref = scran_variances::model_gene_variances_blocked(sub, subblocks.data(), opt);
    EXPECT_EQ(ref.average.variances, res.average.variances);
    EXPECT_EQ(ref.average.fitted, res.average.fitted);
    EXPECT_EQ(ref.average.residuals, res.average.residuals);

    // Computing the expected mean by indexing block levels to remove all empty blocks. 
    std::vector<int> mblocks(nc);
    mblocks[0] = 1;
    mblocks[1] = 2;
    mblocks[2] = 3;
    mblocks[3] = 3;

    auto mref = scran_variances::model_gene_variances_blocked(*dense_row, mblocks.data(), opt);
    EXPECT_EQ(mref.average.means, res.average.means);
}

TEST(ModelGeneVariances, NullAverages) {
    // Get some test coverage for the case where the Buffer::average pointers
    // null and thus should be skipped regardless of what average_policy says.

    int nr = 10, nc = 6;
    auto vec = scran_tests::simulate_vector(nr * nc, []{
        scran_tests::SimulationParameters sparams;
        sparams.density = 0.1;
        sparams.seed = 69;
        sparams.lower = 0;
        sparams.upper = 5;
        return sparams;
    }());
    tatami::DenseMatrix<double, int, decltype(vec)> mat(nr, nc, std::move(vec), true);

    int nblocks = 3;
    std::vector<int> blocks(nc);
    std::fill_n(blocks.begin() + 2, 2, 1);
    std::fill_n(blocks.begin() + 4, 2, 2);

    scran_variances::ModelGeneVariancesBlockedResults<double> output(nr, nblocks, false);
    scran_variances::ModelGeneVariancesBlockedBuffers<double> buffers;
    sanisizer::resize(buffers.per_block, nblocks);
    for (decltype(nblocks) b = 0; b < nblocks; ++b) {
        auto& current = buffers.per_block[b];
        current.means = output.per_block[b].means.data();
        current.variances = output.per_block[b].variances.data();
        current.fitted = output.per_block[b].fitted.data();
        current.residuals = output.per_block[b].residuals.data();
    }

    buffers.average.means = NULL;
    buffers.average.variances = NULL;
    buffers.average.fitted = NULL;
    buffers.average.residuals = NULL;

    scran_variances::ModelGeneVariancesOptions opt;
    opt.average_policy = scran_variances::AveragePolicy::MEAN;
    scran_variances::model_gene_variances_blocked(mat, blocks.data(), buffers, opt);
    opt.average_policy = scran_variances::AveragePolicy::QUANTILE;
    scran_variances::model_gene_variances_blocked(mat, blocks.data(), buffers, opt);

    for (decltype(nblocks) b = 0; b < nblocks; ++b) {
        EXPECT_FALSE(std::isnan(buffers.per_block[b].means[0]));
        EXPECT_FALSE(std::isnan(buffers.per_block[b].variances[1]));
        EXPECT_FALSE(std::isnan(buffers.per_block[b].fitted[2]));
        EXPECT_FALSE(std::isnan(buffers.per_block[b].residuals[4]));
    }
}
