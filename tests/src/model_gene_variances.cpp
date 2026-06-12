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
            scran_tests::SimulateVectorParameters sparams;
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
        auto comp = tatami_stats::variance(true, *dense_row, {});
        scran_tests::compare_almost_equal_containers(ref.mean, comp.mean, {});
        scran_tests::compare_almost_equal_containers(ref.variance, comp.variance, {});

        for (auto f : ref.fitted) {
            EXPECT_TRUE(f > 0);
        }

        int nonzero = 0; 
        for (auto f : ref.residual) {
            nonzero += (f != 0);
        }
        EXPECT_TRUE(nonzero > 0); // there is at least one non-zero residual; but we can't expect this of everyone.

    } else {
        // Checking against the same call, but parallelized.
        auto res1 = scran_variances::model_gene_variances(*dense_row, opt);
        EXPECT_EQ(ref.mean, res1.mean);
        EXPECT_EQ(ref.variance, res1.variance);
    }

    // Almost equal, as there are minor differences due to numerical imprecision.
    auto res2 = scran_variances::model_gene_variances(*dense_column, opt);
    scran_tests::compare_almost_equal_containers(ref.mean, res2.mean, {});
    scran_tests::compare_almost_equal_containers(ref.variance, res2.variance, {});

    auto res3 = scran_variances::model_gene_variances(*sparse_row, opt);
    scran_tests::compare_almost_equal_containers(ref.mean, res3.mean, {});
    scran_tests::compare_almost_equal_containers(ref.variance, res3.variance, {});

    auto res4 = scran_variances::model_gene_variances(*sparse_column, opt);
    scran_tests::compare_almost_equal_containers(ref.mean, res4.mean, {});
    scran_tests::compare_almost_equal_containers(ref.variance, res4.variance, {});
}

TEST_P(ModelGeneVariancesTest, Blocked) {
    const int nblocks = 3;
    const int NC = dense_row->ncol();
    std::vector<int> blocks(NC);
    for (int i = 0; i < NC; ++i) {
        blocks[i] = i % nblocks;
    }

    scran_variances::ModelGeneVariancesOptions opt;
    opt.block_average_policy = scran_variances::BlockAveragePolicy::NONE;
    auto ref = scran_variances::model_gene_variances_blocked(*dense_row, blocks.data(), nblocks, opt);
    EXPECT_TRUE(ref.average.mean.empty());
    EXPECT_TRUE(ref.average.variance.empty());

    auto nthreads = GetParam();
    opt.num_threads = nthreads;

    if (nthreads == 1) {
        // Cursory checks.
        EXPECT_EQ(ref.per_block.size(), nblocks);
    } else {
        // Checking against the same call, but parallelized.
        auto res1 = scran_variances::model_gene_variances_blocked(*dense_row, blocks.data(), nblocks, opt);
        for (int i = 0; i < nblocks; ++i) {
            EXPECT_EQ(ref.per_block[i].mean, res1.per_block[i].mean);
            EXPECT_EQ(ref.per_block[i].variance, res1.per_block[i].variance);
            EXPECT_EQ(ref.per_block[i].fitted, res1.per_block[i].fitted);
            EXPECT_EQ(ref.per_block[i].residual, res1.per_block[i].residual);
        }
    }

    auto res2 = scran_variances::model_gene_variances_blocked(*dense_column, blocks.data(), nblocks, opt);
    for (int i = 0; i < nblocks; ++i) {
        scran_tests::compare_almost_equal_containers(ref.per_block[i].mean, res2.per_block[i].mean, {});
        scran_tests::compare_almost_equal_containers(ref.per_block[i].variance, res2.per_block[i].variance, {});
    }

    auto res3 = scran_variances::model_gene_variances_blocked(*sparse_row, blocks.data(), nblocks, opt);
    for (int i = 0; i < nblocks; ++i) {
        scran_tests::compare_almost_equal_containers(ref.per_block[i].mean, res3.per_block[i].mean, {});
        scran_tests::compare_almost_equal_containers(ref.per_block[i].variance, res3.per_block[i].variance, {});
    }

    auto res4 = scran_variances::model_gene_variances_blocked(*sparse_column, blocks.data(), nblocks, opt);
    for (int i = 0; i < nblocks; ++i) {
        scran_tests::compare_almost_equal_containers(ref.per_block[i].mean, res3.per_block[i].mean, {});
        scran_tests::compare_almost_equal_containers(ref.per_block[i].variance, res3.per_block[i].variance, {});
    }
}

TEST_P(ModelGeneVariancesTest, BlockedMean) {
    const int nblocks = 3;
    const int NC = dense_row->ncol();
    std::vector<int> blocks(NC);
    for (int i = 0; i < NC; ++i) {
        blocks[i] = i % nblocks;
    }

    scran_variances::ModelGeneVariancesOptions opt;
    opt.block_average_policy = scran_variances::BlockAveragePolicy::MEAN;

    // Checking averages with equiweighting.
    opt.block_weight_policy = scran_blocks::WeightPolicy::EQUAL;
    {
        auto ares = scran_variances::model_gene_variances_blocked(*dense_row, blocks.data(), nblocks, opt);
        std::vector<double> expected_means(dense_row->nrow()), 
            expected_variances(dense_row->nrow()),
            expected_fitted(dense_row->nrow()),
            expected_residuals(dense_row->nrow());

        for (std::size_t r = 0, rend = expected_means.size(); r < rend; ++r) {
            for (int b = 0; b < nblocks; ++b) {
                expected_means[r] += ares.per_block[b].mean[r];
                expected_variances[r] += ares.per_block[b].variance[r];
                expected_fitted[r] += ares.per_block[b].fitted[r];
                expected_residuals[r] += ares.per_block[b].residual[r];
            }

            expected_means[r] /= nblocks;
            expected_variances[r] /= nblocks;
            expected_fitted[r] /= nblocks;
            expected_residuals[r] /= nblocks;
        }

        scran_tests::compare_almost_equal_containers(expected_means, ares.average.mean, {});
        scran_tests::compare_almost_equal_containers(expected_variances, ares.average.variance, {});
        scran_tests::compare_almost_equal_containers(expected_fitted, ares.average.fitted, {});
        scran_tests::compare_almost_equal_containers(expected_residuals, ares.average.residual, {});

        // Checking limit of the variable policy.
        opt.block_weight_policy = scran_blocks::WeightPolicy::VARIABLE;
        opt.variable_block_weight_parameters.lower_bound = 0;
        opt.variable_block_weight_parameters.upper_bound = 0;

        auto vres = scran_variances::model_gene_variances_blocked(*dense_row, blocks.data(), nblocks, opt);
        scran_tests::compare_almost_equal_containers(ares.average.mean, vres.average.mean, {});
        scran_tests::compare_almost_equal_containers(ares.average.variance, vres.average.variance, {});
        scran_tests::compare_almost_equal_containers(ares.average.fitted, vres.average.fitted, {});
        scran_tests::compare_almost_equal_containers(ares.average.residual, vres.average.residual, {});
    }

    // Checking averages without equiweighting.
    opt.block_weight_policy = scran_blocks::WeightPolicy::SIZE;
    {
        auto ares = scran_variances::model_gene_variances_blocked(*dense_row, blocks.data(), nblocks, opt);
        std::vector<int> block_size(nblocks);
        for (int c = 0; c < NC; ++c) {
            block_size[blocks[c]] += 1;
        }

        std::vector<double> expected_means(dense_row->nrow()), 
            expected_variances(dense_row->nrow()),
            expected_fitted(dense_row->nrow()),
            expected_residuals(dense_row->nrow());

        for (std::size_t r = 0, rend = expected_means.size(); r < rend; ++r) {
            for (int b = 0; b < nblocks; ++b) {
                expected_means[r] += ares.per_block[b].mean[r] * block_size[b];
                expected_variances[r] += ares.per_block[b].variance[r] * block_size[b];
                expected_fitted[r] += ares.per_block[b].fitted[r] * block_size[b];
                expected_residuals[r] += ares.per_block[b].residual[r] * block_size[b];
            }

            expected_means[r] /= blocks.size();
            expected_variances[r] /= blocks.size();
            expected_fitted[r] /= blocks.size();
            expected_residuals[r] /= blocks.size();
        }

        scran_tests::compare_almost_equal_containers(expected_means, ares.average.mean, {});
        scran_tests::compare_almost_equal_containers(expected_variances, ares.average.variance, {});
        scran_tests::compare_almost_equal_containers(expected_fitted, ares.average.fitted, {});
        scran_tests::compare_almost_equal_containers(expected_residuals, ares.average.residual, {});

        // Checking limit of the variable policy.
        opt.block_weight_policy = scran_blocks::WeightPolicy::VARIABLE;
        opt.variable_block_weight_parameters.lower_bound = 0;
        opt.variable_block_weight_parameters.upper_bound = 100000;

        auto vres = scran_variances::model_gene_variances_blocked(*dense_row, blocks.data(), nblocks, opt);
        scran_tests::compare_almost_equal_containers(ares.average.mean, vres.average.mean, {});
        scran_tests::compare_almost_equal_containers(ares.average.variance, vres.average.variance, {});
        scran_tests::compare_almost_equal_containers(ares.average.fitted, vres.average.fitted, {});
        scran_tests::compare_almost_equal_containers(ares.average.residual, vres.average.residual, {});
    }
}

TEST_P(ModelGeneVariancesTest, BlockedMedian) {
    const int nblocks = 3;
    const int NC = dense_row->ncol();
    std::vector<int> blocks(NC);
    for (int i = 0; i < NC; ++i) {
        blocks[i] = i % nblocks;
    }

    scran_variances::ModelGeneVariancesOptions opt;
    opt.block_average_policy = scran_variances::BlockAveragePolicy::QUANTILE;
    auto ares = scran_variances::model_gene_variances_blocked(*dense_row, blocks.data(), nblocks, opt);

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
            buffer_means.push_back(ares.per_block[b].mean[r]);
            buffer_variances.push_back(ares.per_block[b].variance[r]);
            buffer_fitted.push_back(ares.per_block[b].fitted[r]);
            buffer_residuals.push_back(ares.per_block[b].residual[r]);
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

    EXPECT_EQ(expected_means, ares.average.mean);
    EXPECT_EQ(expected_variances, ares.average.variance);
    EXPECT_EQ(expected_fitted, ares.average.fitted);
    EXPECT_EQ(expected_residuals, ares.average.residual);
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
            scran_tests::SimulateVectorParameters sparams;
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

    auto res = scran_variances::model_gene_variances_blocked(*dense_row, blocks.data(), 5, opt);
    EXPECT_FALSE(std::isnan(res.per_block[0].mean[0]));
    EXPECT_FALSE(std::isnan(res.per_block[0].variance[0]));
    EXPECT_FALSE(std::isnan(res.per_block[1].mean[0]));
    EXPECT_TRUE(std::isnan(res.per_block[1].variance[0]));
    EXPECT_TRUE(std::isnan(res.per_block[2].mean[0]));
    EXPECT_TRUE(std::isnan(res.per_block[2].variance[0]));
    EXPECT_FALSE(std::isnan(res.per_block[3].mean[0]));
    EXPECT_TRUE(std::isnan(res.per_block[3].variance[0]));
    EXPECT_FALSE(std::isnan(res.per_block[4].mean[0]));
    EXPECT_FALSE(std::isnan(res.per_block[4].variance[0]));

    // Subset to blocks 0 and 4 to ignore all blocks with fewer than two cells. This means we do [2, ncol).
    tatami::DelayedSubsetBlock<double, int> sub(dense_row, 2, nc - 2, false);
    std::vector<int> subblocks(nc - 2);
    subblocks[0] = 1;
    subblocks[1] = 1;
    auto ref = scran_variances::model_gene_variances_blocked(sub, subblocks.data(), 2, opt);
    EXPECT_EQ(ref.average.variance, res.average.variance);
    EXPECT_EQ(ref.average.fitted, res.average.fitted);
    EXPECT_EQ(ref.average.residual, res.average.residual);

    // Computing the expected mean. As the weight policy is SIZE, we basically just take the row means.
    auto expected = tatami_stats::sum(true, *dense_row, {});
    for (auto& x : expected) {
        x /= nc;
    }
    scran_tests::compare_almost_equal_containers(expected, res.average.mean, {});
}

TEST_F(ModelGeneVariancesNearEmptyBlockTest, Median) {
    scran_variances::ModelGeneVariancesOptions opt;
    opt.block_average_policy = scran_variances::BlockAveragePolicy::QUANTILE;
    auto res = scran_variances::model_gene_variances_blocked(*dense_row, blocks.data(), 5, opt);

    // Subset to blocks 0 and 4 to ignore all blocks with fewer than two cells. This means we do [2, ncol).
    tatami::DelayedSubsetBlock<double, int> sub(dense_row, 2, nc - 2, false);
    std::vector<int> subblocks(nc - 2);
    subblocks[0] = 1;
    subblocks[1] = 1;

    auto ref = scran_variances::model_gene_variances_blocked(sub, subblocks.data(), 2, opt);
    EXPECT_EQ(ref.average.variance, res.average.variance);
    EXPECT_EQ(ref.average.fitted, res.average.fitted);
    EXPECT_EQ(ref.average.residual, res.average.residual);

    // Computing the expected mean by indexing block levels to remove all empty blocks. 
    std::vector<int> mblocks(nc);
    mblocks[0] = 1;
    mblocks[1] = 2;
    mblocks[2] = 3;
    mblocks[3] = 3;

    auto mref = scran_variances::model_gene_variances_blocked(*dense_row, mblocks.data(), 4, opt);
    EXPECT_EQ(mref.average.mean, res.average.mean);
}

TEST(ModelGeneVariances, NullAverages) {
    // Get some test coverage for the case where the Buffer::average pointers
    // null and thus should be skipped regardless of what block_average_policy says.

    int nr = 10, nc = 6;
    auto vec = scran_tests::simulate_vector(nr * nc, []{
        scran_tests::SimulateVectorParameters sparams;
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

    scran_variances::ModelGeneVariancesBlockedResults<double> output(nr, nblocks, /* average = */ false, /* trend = */ true);
    scran_variances::ModelGeneVariancesBlockedBuffers<double> buffers;
    sanisizer::resize(buffers.per_block, nblocks);
    for (decltype(nblocks) b = 0; b < nblocks; ++b) {
        auto& current = buffers.per_block[b];
        current.mean = output.per_block[b].mean.data();
        current.variance = output.per_block[b].variance.data();
        current.fitted = output.per_block[b].fitted.data();
        current.residual = output.per_block[b].residual.data();
    }

    buffers.average.mean = NULL;
    buffers.average.variance = NULL;
    buffers.average.fitted = NULL;
    buffers.average.residual = NULL;

    scran_variances::ModelGeneVariancesOptions opt;
    opt.block_average_policy = scran_variances::BlockAveragePolicy::MEAN;
    scran_variances::model_gene_variances_blocked(mat, blocks.data(), 3, buffers, opt);
    opt.block_average_policy = scran_variances::BlockAveragePolicy::QUANTILE;
    scran_variances::model_gene_variances_blocked(mat, blocks.data(), 3, buffers, opt);

    for (decltype(nblocks) b = 0; b < nblocks; ++b) {
        EXPECT_FALSE(std::isnan(buffers.per_block[b].mean[0]));
        EXPECT_FALSE(std::isnan(buffers.per_block[b].variance[1]));
        EXPECT_FALSE(std::isnan(buffers.per_block[b].fitted[2]));
        EXPECT_FALSE(std::isnan(buffers.per_block[b].residual[4]));
    }
}

TEST(ModelGeneVariances, NoTrend) {
    int nr = 10, nc = 6;
    auto vec = scran_tests::simulate_vector(nr * nc, []{
        scran_tests::SimulateVectorParameters sparams;
        sparams.density = 0.1;
        sparams.seed = 69;
        sparams.lower = 0;
        sparams.upper = 5;
        return sparams;
    }());
    tatami::DenseMatrix<double, int, decltype(vec)> mat(nr, nc, std::move(vec), true);

    {
        scran_variances::ModelGeneVariancesOptions opt;
        auto ref = scran_variances::model_gene_variances(mat, opt);

        opt.trend = false;
        auto res = scran_variances::model_gene_variances(mat, opt);
        EXPECT_EQ(ref.mean, res.mean);
        EXPECT_EQ(ref.variance, res.variance);
        EXPECT_TRUE(res.fitted.empty());
        EXPECT_TRUE(res.residual.empty());
    }

    std::vector<int> block { 0, 1, 0, 1, 0, 1 };

    // Check blocked as well.
    {
        scran_variances::ModelGeneVariancesOptions opt;
        auto ref = scran_variances::model_gene_variances_blocked(mat, block.data(), 2, opt);

        opt.trend = false;
        auto res = scran_variances::model_gene_variances_blocked(mat, block.data(), 2, opt);
        EXPECT_EQ(ref.average.mean, res.average.mean);
        EXPECT_EQ(ref.average.variance, res.average.variance);
        EXPECT_TRUE(res.average.fitted.empty());
        EXPECT_TRUE(res.average.residual.empty());

        ASSERT_EQ(res.per_block.size(), 2);
        for (int b = 0; b < 2; ++b) {
            EXPECT_EQ(ref.per_block[b].mean, res.per_block[b].mean);
            EXPECT_EQ(ref.per_block[b].variance, res.per_block[b].variance);
            EXPECT_TRUE(res.per_block[b].fitted.empty());
            EXPECT_TRUE(res.per_block[b].residual.empty());
        }
    }

    // Check that an error is correctly thrown if average fit/residuals are requested without per-block values.
    {
        scran_variances::ModelGeneVariancesBlockedResults<double> output(nr, 2, /* average = */ false, /* trend = */ false);
        scran_variances::ModelGeneVariancesBlockedBuffers<double> buffers;
        sanisizer::resize(buffers.per_block, 2);
        for (int b = 0; b < 2; ++b) {
            auto& current = buffers.per_block[b];
            current.mean = output.per_block[b].mean.data();
            current.variance = output.per_block[b].variance.data();
            current.fitted = NULL;
            current.residual = NULL;
        }

        scran_variances::ModelGeneVariancesResults<double> averages(nr, /* trend = */ true);
        scran_variances::ModelGeneVariancesBuffers<double> abuffers;
        abuffers.mean = averages.mean.data();
        abuffers.variance = averages.mean.data();
        abuffers.fitted = averages.fitted.data();
        abuffers.residual = averages.residual.data();
        buffers.average = abuffers; 

        std::string msg;
        try {
            scran_variances::model_gene_variances_blocked(mat, block.data(), 2, buffers, {});
        } catch (std::exception& e) {
            msg = e.what();
        }
        EXPECT_TRUE(msg.find("per-block trend fits") != std::string::npos);
    }
}

TEST(ModelGeneVariances, BlockedMismatch) {
    tatami::DenseMatrix<double, int, std::vector<double> > mat(10, 0, std::vector<double>(), true);
    scran_variances::ModelGeneVariancesBlockedBuffers<double> buffers;
    std::string msg;
    try {
        scran_variances::model_gene_variances_blocked(mat, static_cast<const int*>(NULL), 3, buffers, {});
    } catch (std::exception& e) {
        msg = e.what();
    }
    EXPECT_TRUE(msg.find("not equal to") != std::string::npos);
}
