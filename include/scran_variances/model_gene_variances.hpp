#ifndef SCRAN_MODEL_GENE_VARIANCES_HPP
#define SCRAN_MODEL_GENE_VARIANCES_HPP

#include <algorithm>
#include <vector>
#include <limits>
#include <cstddef>
#include <cassert>
#include <optional>

#include "tatami/tatami.hpp"
#include "tatami_stats/tatami_stats.hpp"
#include "scran_blocks/scran_blocks.hpp"
#include "sanisizer/sanisizer.hpp"

#include "fit_variance_trend.hpp"
#include "utils.hpp"

/**
 * @file model_gene_variances.hpp
 * @brief Model the per-gene variances. 
 */

namespace scran_variances {

/**
 * Policy for averaging statistics across blocks.
 *
 * - `MEAN`: weighted mean, where weights are computed using `scran_blocks::compute_weights()`.
 * - `QUANTILE`: quantile, defaulting to 50%, a.k.a., the median.
 * - `NONE`: do not report any inter-block average. 
 */
enum class BlockAveragePolicy : unsigned char { MEAN, QUANTILE, NONE };

/**
 * @brief Options for `model_gene_variances()` and friends.
 */
struct ModelGeneVariancesOptions {
    /**
     * Whether to fit a mean-variance trend.
     * It can be occasionally desirable to disable the trend fitting if the trend is to be obtained from some other source,
     * e.g., using a separate set of spike-in transcripts to estimate the technical component of each gene's variance.
     */
    bool trend = true;

    /**
     * Options for fitting the mean-variance trend.
     * Ignored if `ModelGeneVariancesOptions::trend = false`.
     */
    FitVarianceTrendOptions fit_variance_trend_options;

    /**
     * Policy to use for averaging statistics across blocks.
     * Only relevant to `model_gene_variances_blocked()`.
     * Ignored for overloads accepting `ModelGeneVariancesBlockedBuffers` where all entries in `average` are `NULL`.
     */
    BlockAveragePolicy block_average_policy = BlockAveragePolicy::MEAN;

    /**
     * Policy for weighting the contribution from each block when computing the mean for each statistic.
     * Only relevant to `model_gene_variances_blocked()` when `ModelGeneVariancesOptions::average_policy = BlockAveragePolicy::MEAN`. 
     *
     * By default, we define equal weights for blocks beyond a certain size (see `ModelGeneVariancesOptions::variable_block_weight_parameters`).
     * For smaller blocks, the weight is linearly proportional to its size to avoid putting too much confidence in estimates from very small blocks.
     *
     * Other options include `scran_blocks::WeightPolicy::EQUAL`, where all blocks are equally weighted regardless of size;
     * and `scran_blocks::WeightPolicy::NONE`, where the contribution of each block is proportional to its size.
     */
    scran_blocks::WeightPolicy block_weight_policy = scran_blocks::WeightPolicy::VARIABLE;

    /**
     * Parameters for the variable block weights, including the threshold at which blocks are considered to be large enough to have equal weight.
     * Only relevant to `model_gene_variances_blocked()` when `ModelGeneVariancesOptions::average_policy = BlockAveragePolicy::MEAN`
     * and `ModelGeneVariancesOptions::block_weight_policy = scran_blocks::WeightPolicy::VARIABLE`.
     */
    scran_blocks::VariableWeightParameters variable_block_weight_parameters; 

    /**
     * @cond
     */
    // Back-compatibility only.
    bool compute_average = true;
    /**
     * @endcond
     */

    /**
     * Quantile to use as an "average" statistic across blocks.
     * Only relevant to `model_gene_variances_blocked()` when `ModelGeneVariancesOptions::average_policy = BlockAveragePolicy::QUANTILE`. 
     * Defaults to 0.5, a.k.a., the median.
     */
    double block_quantile = 0.5;

    /**
     * Number of threads to use for the variance calculations and trend fitting. 
     * The parallelization scheme is defined by `tatami::parallelize()`. 
     */
    int num_threads = 1;
};

/**
 * @brief Buffers for `model_gene_variances()` and friends.
 * @tparam Stat_ Floating-point type of the output statistics.
 *
 * In general, the pointers in this class should _not_ be set to `NULL`.
 * The only exception is for instances of this class that are used as `ModelGeneVariancesBlockedBuffers::average`,
 * where setting the pointer to `NULL` will omit calculation of the corresponding average statistic.
 */
template<typename Stat_>
struct ModelGeneVariancesBuffers {
    /**
     * Pointer to an array of length equal to the number of genes, to be filled with the mean log-expression for each gene.
     */
    Stat_* mean;

    /**
     * Pointer to an array of length equal to the number of genes, containing the variance in the log-expression for each gene.
     */
    Stat_* variance;

    /**
     * Pointer to an array of length equal to the number of genes, containing the fitted value of the mean-variance trend for each gene.
     *
     * If this or `ModelGeneVariancesBuffers::residual` is `NULL`, no trend is fitted.
     */
    Stat_* fitted;

    /**
     * Vector of length equal to the number of genes, containing the residuals of the mean-variance trend for each gene.
     *
     * If this or `ModelGeneVariancesBuffers::fitted` is `NULL`, no trend is fitted.
     */
    Stat_* residual;
};

/**
 * @cond
 */
template<typename Value_, typename Index_, typename Stat_>
void model_gene_variances_direct(
    const tatami::Matrix<Value_, Index_>& mat,
    Stat_* const output_mean,
    Stat_* const output_variance,
    const int num_threads
) {
    const auto NR = mat.nrow();
    const auto NC = mat.ncol();

    if (mat.sparse()) {
        tatami::Options topt;
        topt.sparse_extract_index = false;

        tatami::parallelize([&](int, Index_ s, Index_ l) -> void {
            auto ext = tatami::consecutive_extractor<true>(mat, true, s, l, topt);
            auto vbuffer = tatami::create_container_of_Index_size<std::vector<Value_> >(NC);
            quickstats::RssWorkspace<Stat_> work;
            for (Index_ x = 0; x < l; ++x) {
                auto out = ext->fetch(vbuffer.data(), NULL);
                const auto res = quickstats::rss(NC, out.number, out.value, work);
                output_mean[x + s] = res.mean;
                output_variance[x + s] = quickstats::rss_to_variance(NC, res.rss);
            }
        }, NR, num_threads);

    } else {
        tatami::parallelize([&](int, Index_ s, Index_ l) -> void {
            auto ext = tatami::consecutive_extractor<false>(mat, true, s, l);
            auto buffer = tatami::create_container_of_Index_size<std::vector<Value_> >(NC);
            quickstats::RssWorkspace<Stat_> work;
            for (Index_ x = 0; x < l; ++x) {
                auto out = ext->fetch(buffer.data());
                const auto res = quickstats::rss(NC, out, work);
                output_mean[x + s] = res.mean;
                output_variance[x + s] = quickstats::rss_to_variance(NC, res.rss);
            }
        }, NR, num_threads);
    }
}

template<typename Value_, typename Index_, typename Stat_>
void model_gene_variances_running(
    const tatami::Matrix<Value_, Index_>& mat,
    Stat_* const output_mean,
    Stat_* const output_variance,
    const int num_threads
) {
    const auto NR = mat.nrow();
    const auto NC = mat.ncol();
    const bool is_sparse = mat.is_sparse();

    const bool do_parallel = num_threads > 1;
    std::optional<std::vector<std::optional<std::vector<Stat_> > > > all_partial_mean, all_partial_rss;
    std::optional<std::vector<Index_> > all_partial_count;
    if (do_parallel) {
        // -1, as we'll repurpose the RSS output buffer to store the partial RSS of the first thread.
        all_partial_rss.emplace(sanisizer::cast<I<decltype(all_partial_rss->size())> >(num_threads - 1));
        all_partial_mean.emplace(sanisizer::cast<I<decltype(all_partial_mean->size())> >(num_threads));
        all_partial_count.emplace(sanisizer::cast<I<decltype(all_partial_count->size())> >(num_threads));
    }

    std::fill_n(output_variance, NR, 0);
    std::fill_n(output_mean, NR, 0);

    const int nused = tatami::parallelize([&](int thread, Index_ s, Index_ l) -> void {
        Stat_* rss_ptr;
        Stat_* mean_ptr;
        std::optional<std::vector<Stat_> > cur_rss, cur_mean;

        if (!do_parallel) {
            // Storing mean and RSS directly in the output vector to cut down two allocations if we're not working in parallel.
            rss_ptr = output_variance;
            mean_ptr = output_mean;
        } else {
            // Storing the partial RSS directly in the output vector to save ourselves an allocation if we're in the first thread.
            // We can't do the same for the mean, though, as we need to keep the partial mean and the global mean separate for the reduction.
            cur_mean.emplace(tatami::cast_Index_to_container_size<std::vector<Stat_> >(NR));
            mean_ptr = cur_mean->data();
            if (thread == 0) {
                rss_ptr = output_variance;
            } else {
                cur_rss.emplace(tatami::cast_Index_to_container_size<std::vector<Stat_> >(NR));
                rss_ptr = cur_rss->data();
            }
        }

        if (is_sparse) {
            tatami::Options topt;
            topt.sparse_ordered_index = false; // we don't care about the ordering of the non-zero elements.
            auto ext = tatami::consecutive_extractor<true>(mat, false, s, l, topt);
            auto vbuffer = tatami::create_container_of_Index_size<std::vector<Value_> >(NR);
            auto ibuffer = tatami::create_container_of_Index_size<std::vector<Index_> >(NR);
            auto nonzeros = tatami::create_container_of_Index_size<std::vector<Index_> >(NR);

            quickstats::RssRunningSparse<Index_, Value_, Stat_> runner(NR, mean_ptr, rss_ptr, nonzeros.data());
            for (Index_ x = 0; x < l; ++x) {
                auto out = ext->fetch(vbuffer.data(), ibuffer.data());
                runner.add(out.number, out.value, out.index);
            }
            runner.finish();

        } else {
            auto ext = tatami::consecutive_extractor<false>(mat, false, s, l);
            auto buffer = tatami::create_container_of_Index_size<std::vector<Value_> >(NR);

            quickstats::RssRunningDense<Value_, Stat_> runner(NR, mean_ptr, rss_ptr);
            for (Index_ x = 0; x < l; ++x) {
                auto out = ext->fetch(buffer.data());
                runner.add(out);
            }
            runner.finish();
        }

        if (do_parallel) {
            (*all_partial_count)[thread] = l;
            (*all_partial_mean)[thread] = std::move(cur_mean);
            if (thread > 0) {
                (*all_partial_rss)[thread - 1] = std::move(cur_rss);
            }
        }
    }, NC, num_threads);

    // Don't check nused > 1, as it's possible for do_parallel = true with nused = 1 if not all threads are used.
    // This would cause us to leave output_mean and output_variance empty.
    if (do_parallel) {
        const auto& ap_count = *all_partial_count;
        const auto& ap_mean = *all_partial_mean;
        const auto& ap_rss = *all_partial_rss;

        // Computing the global mean.
        for (int u = 0; u < nused; ++u) {
            const Stat_ mult = static_cast<Stat_>(ap_count[u]) / static_cast<Stat_>(NC);
            const auto& cur_mean = *(ap_mean[u]);
            for (Index_ d = 0; d < NR; ++d) {
                output_mean[d] += cur_mean[d] * mult;
            }
        }

        // Combining the RSS. We can use recenter_rss_unsafe() as we are guaranteed that cur_count > 0,
        // as parallelize() will only ever split into non-empty ranges if those ranges are used.
        for (int u = 0; u < nused; ++u) {
            const auto cur_count = ap_count[u];
            const auto& cur_mean = *(ap_mean[u]);
            if (u == 0) {
                for (Index_ d = 0; d < NR; ++d) {
                    output_variance[d] = quickstats::recenter_rss_unsafe(cur_count, output_variance[d], cur_mean[d], output_mean[d]); 
                }
            } else {
                const auto& cur_rss = *(ap_rss[u - 1]);
                for (Index_ d = 0; d < NR; ++d) {
                    output_variance[d] += quickstats::recenter_rss_unsafe(cur_count, cur_rss[d], cur_mean[d], output_mean[d]); 
                }
            }
        }
    }

    if (NC == 0) {
        // If NC == 0, nused = 0 and RssRunning::finish() doesn't even get called.
        // So, we have to fill the output means with NaNs manually.
        std::fill_n(output_mean, NR, std::numeric_limits<Stat_>::quiet_NaN());
    }
    quickstats::rss_to_variance(NR, NC, output_variance);
}
/**
 * @endcond
 */

/** 
 * Model the per-gene variances as a function of the mean in single-cell expression data.
 * We compute the mean and variance for each gene and fit a trend to the variances with respect to the means using `fit_variance_trend()`.
 * We assume that most genes at any given abundance are not highly variable, such that the fitted value of the trend is interpreted as the "uninteresting" variance - 
 * this is mostly attributed to technical variation like sequencing noise, but can also represent constitutive biological noise like transcriptional bursting.
 * Under this assumption, the residual can be treated as a measure of biologically interesting variation.
 * Genes with large residuals can then be selected for downstream analyses, e.g., with `choose_highly_variable_genes()`.
 *
 * @tparam Value_ Data type of the matrix.
 * @tparam Index_ Integer type of the row/column indices.
 * @tparam Stat_ Floating-point type of the output statistics.
 *
 * @param mat Matrix of expression values, typically after normalization and log-transformation.
 * Rows should be genes while columns should be cells.
 * @param buffers Collection of buffers in which to store the computed statistics.
 * @param options Further options.
 */
template<typename Value_, typename Index_, typename Stat_> 
void model_gene_variances(
    const tatami::Matrix<Value_, Index_>& mat, 
    const ModelGeneVariancesBuffers<Stat_> buffers,
    const ModelGeneVariancesOptions& options
) {
    if (mat.prefer_rows()) {
        model_gene_variances_direct(mat, buffers.mean, buffers.variance, options.num_threads);
    } else {
        model_gene_variances_running(mat, buffers.mean, buffers.variance, options.num_threads);
    }

    FitVarianceTrendWorkspace<Stat_> work;
    auto fopt = options.fit_variance_trend_options;
    fopt.num_threads = options.num_threads;

    if (buffers.fitted != NULL && buffers.residual != NULL) {
        const auto NR = mat.nrow();
        if (mat.ncol() >= 2) {
            fit_variance_trend(NR, buffers.mean, buffers.variance, buffers.fitted, buffers.residual, work, fopt);
        } else {
            std::fill_n(buffers.fitted, NR, std::numeric_limits<double>::quiet_NaN());
            std::fill_n(buffers.residual, NR, std::numeric_limits<double>::quiet_NaN());
        }
    }
}

/**
 * @brief Results of `model_gene_variances()`. 
 * @tparam Stat_ Floating-point type of the output statistics.
 */
template<typename Stat_>
struct ModelGeneVariancesResults {
    /**
     * @cond
     */
    ModelGeneVariancesResults() = default;

    ModelGeneVariancesResults(const std::size_t ngenes, const bool trend) :
        mean(sanisizer::cast<I<decltype(mean.size())> >(ngenes)
#ifdef SCRAN_VARIANCES_TEST_INIT
            , SCRAN_VARIANCES_TEST_INIT
#endif
        ),
        variance(sanisizer::cast<I<decltype(variance.size())> >(ngenes)
#ifdef SCRAN_VARIANCES_TEST_INIT
            , SCRAN_VARIANCES_TEST_INIT
#endif
        ),
        fitted(sanisizer::cast<I<decltype(fitted.size())> >(trend ? ngenes : 0)
#ifdef SCRAN_VARIANCES_TEST_INIT
            , SCRAN_VARIANCES_TEST_INIT
#endif
        ),
        residual(sanisizer::cast<I<decltype(residual.size())> >(trend ? ngenes : 0)
#ifdef SCRAN_VARIANCES_TEST_INIT
            , SCRAN_VARIANCES_TEST_INIT
#endif
        )
    {}
    /**
     * @endcond
     */

    /**
     * Vector of length equal to the number of genes, containing the mean log-expression for each gene.
     */
    std::vector<Stat_> mean;

    /**
     * Vector of length equal to the number of genes, containing the variance in the log-expression for each gene.
     */
    std::vector<Stat_> variance;

    /**
     * Vector of length equal to the number of genes, containing the fitted value of the mean-variance trend for each gene.
     *
     * This will be empty if `ModelGeneVariancesOptions::trend = false`.
     */
    std::vector<Stat_> fitted;

    /**
     * Vector of length equal to the number of genes, containing the residuals of the mean-variance trend for each gene.
     *
     * This will be empty if `ModelGeneVariancesOptions::trend = false`.
     */
    std::vector<Stat_> residual;
};

/** 
 * Overload of `model_gene_variances()` that allocates space for the output statistics.
 *
 * @tparam Stat_ Floating-point type of the output statistics.
 * @tparam Value_ Data type of the matrix.
 * @tparam Index_ Integer type of the row/column indices.
 *
 * @param mat Matrix of expression values, typically after normalization and log-transformation.
 * Rows should be genes while columns should be cells.
 * @param options Further options.
 *
 * @return Results of the variance modelling.
 */
template<typename Stat_ = double, typename Value_, typename Index_>
ModelGeneVariancesResults<Stat_> model_gene_variances(const tatami::Matrix<Value_, Index_>& mat, const ModelGeneVariancesOptions& options) {
    ModelGeneVariancesResults<Stat_> output(mat.nrow(), options.trend); // cast is safe, as any tatami Index_ can always fit into a size_t.

    ModelGeneVariancesBuffers<Stat_> buffers;
    buffers.mean = output.mean.data();
    buffers.variance = output.variance.data();

    if (options.trend) {
        buffers.fitted = output.fitted.data();
        buffers.residual = output.residual.data();
    } else {
        buffers.fitted = NULL;
        buffers.residual = NULL;
    }

    model_gene_variances(mat, std::move(buffers), options);
    return output;
}

/**
 * @brief Buffers for `model_gene_variances_blocked()`.
 * @tparam Stat_ Floating-point type of the output statistics.
 */
template<typename Stat_>
struct ModelGeneVariancesBlockedBuffers {
    /**
     * Vector of length equal to the number of blocks,
     * where each entry contains the buffers to store the variance modelling results for a single block.
     */
    std::vector<ModelGeneVariancesBuffers<Stat_> > per_block;

    /**
     * Buffers to store the average across blocks for all statistics in `per_block`.
     *
     * Any of the pointers in this object may be `NULL`, in which case the corresponding average is not computed.
     *
     * If `average.fitted` or `average.residual` is not `NULL`, all entries of `per_block` should have non-`NULL` pointers for their own `fitted` and `residual`.
     * (That is, a mean-variance trend should have been fitted in each block.) 
     * Otherwise, an error will be thrown.
     */
    ModelGeneVariancesBuffers<Stat_> average;
};

/**
 * @brief Results of `model_gene_variances_blocked()`.
 * @tparam Stat_ Floating-point type of the output statistics.
 */
template<typename Stat_>
struct ModelGeneVariancesBlockedResults {
    /**
     * @cond
     */
    ModelGeneVariancesBlockedResults() = default;

    ModelGeneVariancesBlockedResults(const std::size_t ngenes, const std::size_t nblocks, const bool do_average, const bool do_trend) :
        average(do_average ? ngenes : 0, do_trend)
    {
        per_block.reserve(nblocks);
        for (I<decltype(nblocks)> b = 0; b < nblocks; ++b) {
            per_block.emplace_back(ngenes, do_trend);
        }
    }
    /**
     * @endcond
     */

    /**
     * Vector of length equal to the number of blocks, where each entry contains the variance modelling results for a single block.
     */
    std::vector<ModelGeneVariancesResults<Stat_> > per_block;

    /**
     * Average across blocks for all statistics in `per_block`.
     * This is only populated if `ModelGeneVariancesOptions::average_policy` is not `BlockAveragePolicy::NONE`.
     */
    ModelGeneVariancesResults<Stat_> average;
};

/**
 * @cond
 */
template<typename Index_, typename Stat_>
void finalize_direct_block_means(
    const std::size_t num_blocks,
    const std::vector<Index_>& block_sizes,
    std::vector<Stat_>& cur_means, 
    const Index_ i,
    const std::vector<ModelGeneVariancesBuffers<Stat_> >& output
) {
    for (std::size_t b = 0; b < num_blocks; ++b) {
        auto& mean = cur_means[b];
        if (block_sizes[b] == 0) {
            mean = std::numeric_limits<Stat_>::quiet_NaN();
        } else {
            mean /= block_sizes[b];
        }
        output[b].mean[i] = mean;
    }
}

template<typename Value_, typename Index_, typename Block_, typename Stat_>
void model_gene_variances_blocked_direct(
    const tatami::Matrix<Value_, Index_>& mat, 
    const Block_* const block, 
    const std::size_t num_blocks,
    const std::vector<Index_>& block_sizes,
    const std::vector<ModelGeneVariancesBuffers<Stat_> >& output, 
    int num_threads
) {
    const auto NR = mat.nrow(), NC = mat.ncol();
    assert(sanisizer::is_equal(num_blocks, output.size()));
    assert(sanisizer::is_equal(num_blocks, block_sizes.size()));

    if (mat.sparse()) {
        tatami::parallelize([&](int, Index_ s, Index_ l) -> void {
            auto ext = tatami::consecutive_extractor<true>(mat, true, s, l);
            auto vbuffer = tatami::create_container_of_Index_size<std::vector<Value_> >(NC);
            auto ibuffer = tatami::create_container_of_Index_size<std::vector<Index_> >(NC);
            auto cur_means = sanisizer::create<std::vector<Stat_> >(num_blocks);
            auto cur_rss = sanisizer::create<std::vector<Stat_> >(num_blocks);
            auto cur_non_zeros = sanisizer::create<std::vector<Index_> >(num_blocks);

            for (Index_ x = 0; x < l; ++x) {
                auto range = ext->fetch(vbuffer.data(), ibuffer.data());

                // Computing the mean first.
                for (Index_ i = 0; i < range.number; ++i) {
                    const auto g = block[range.index[i]];
                    cur_means[g] += range.value[i];
                    ++cur_non_zeros[g];
                }
                finalize_direct_block_means(num_blocks, block_sizes, cur_means, static_cast<Index_>(s + x), output);

                // Now computing the RSS.
                for (Index_ i = 0; i < range.number; ++i) {
                    const auto g = block[range.index[i]];
                    const auto delta = range.value[i] - cur_means[g];
                    cur_rss[g] += delta * delta;
                }
                for (std::size_t b = 0; b < num_blocks; ++b) {
                    const Stat_ my_rss = cur_rss[b] + cur_means[b] * cur_means[b] * (block_sizes[b] - cur_non_zeros[b]);
                    output[b].variance[s + x] = quickstats::rss_to_variance(block_sizes[b], my_rss);
                }

                std::fill(cur_means.begin(), cur_means.end(), 0);
                std::fill(cur_rss.begin(), cur_rss.end(), 0);
                std::fill(cur_non_zeros.begin(), cur_non_zeros.end(), 0);
            }
        }, NR, num_threads);

    } else {
        tatami::parallelize([&](int, Index_ s, Index_ l) -> void {
            auto ext = tatami::consecutive_extractor<false>(mat, true, s, l);
            auto buffer = tatami::create_container_of_Index_size<std::vector<Value_> >(NC);
            auto cur_means = sanisizer::create<std::vector<Stat_> >(num_blocks);
            auto cur_rss = sanisizer::create<std::vector<Stat_> >(num_blocks);

            for (Index_ x = 0; x < l; ++x) {
                auto ptr = ext->fetch(buffer.data());

                // Computing the mean first.
                for (Index_ j = 0; j < NC; ++j) {
                    cur_means[block[j]] += ptr[j];
                }
                finalize_direct_block_means(num_blocks, block_sizes, cur_means, static_cast<Index_>(s + x), output);

                // Now computing the RSS.
                for (Index_ j = 0; j < NC; ++j) {
                    const auto g = block[j];
                    const auto delta = ptr[j] - cur_means[g];
                    cur_rss[g] += delta * delta;
                }
                for (std::size_t b = 0; b < num_blocks; ++b) {
                    output[b].variance[s + x] = quickstats::rss_to_variance(block_sizes[b], cur_rss[b]);
                }

                std::fill(cur_means.begin(), cur_means.end(), 0);
                std::fill(cur_rss.begin(), cur_rss.end(), 0);
            }
        }, NR, num_threads);
    }
}

template<typename Value_, typename Index_, typename Group_, typename Stat_>
void model_gene_variances_blocked_running(
    const tatami::Matrix<Value_, Index_>& mat,
    const Group_* const block, 
    const std::size_t num_blocks,
    const std::vector<Index_>& block_sizes,
    const std::vector<ModelGeneVariancesBuffers<Stat_> >& output, 
    int num_threads
) {
    const auto NR = mat.nrow(), NC = mat.ncol();
    const bool is_sparse = mat.is_sparse();

    const bool do_parallel = num_threads > 1;
    std::optional<std::vector<std::optional<std::vector<std::vector<Stat_> > > > > all_partial_mean, all_partial_rss;
    std::optional<std::vector<std::optional<std::vector<Index_> > > > all_partial_count;
    if (do_parallel) {
        // -1, as we'll repurpose the RSS output buffer to store the partial RSS of the first thread.
        all_partial_rss.emplace(sanisizer::cast<I<decltype(all_partial_rss->size())> >(num_threads - 1));
        all_partial_mean.emplace(sanisizer::cast<I<decltype(all_partial_mean->size())> >(num_threads));
        all_partial_count.emplace(sanisizer::cast<I<decltype(all_partial_count->size())> >(num_threads));
    }

    assert(sanisizer::is_equal(num_blocks, output.size()));
    for (std::size_t b = 0; b < num_blocks; ++b) {
        std::fill_n(output[b].mean, NR, 0);
        std::fill_n(output[b].variance, NR, 0);
    }

    const int nused = tatami::parallelize([&](int thread, Index_ s, Index_ l) -> void {
        auto tmp_mean_ptrs = sanisizer::create<std::vector<Stat_*> >(num_blocks);
        auto tmp_rss_ptrs = sanisizer::create<std::vector<Stat_*> >(num_blocks);
        std::optional<std::vector<std::vector<Stat_> > > cur_mean, cur_rss;

        if (!do_parallel) {
            // Storing mean and RSS directly in the output vector to cut down two allocations if we're not working in parallel.
            for (std::size_t b = 0; b < num_blocks; ++b) {
                tmp_mean_ptrs[b] = output[b].mean;
                tmp_rss_ptrs[b] = output[b].variance;
            }

        } else {
            // Storing the partial RSS directly in the output vectors to save ourselves an allocation if we're in the first thread.
            // We can't do the same for the mean, though, as we need to keep the partial mean and the global mean separate for the reduction.
            cur_mean.emplace(sanisizer::cast<I<decltype(cur_mean->size())> >(num_blocks));
            for (std::size_t b = 0; b < num_blocks; ++b) {
                tatami::resize_container_to_Index_size((*cur_mean)[b], NR);
                tmp_mean_ptrs[b] = (*cur_mean)[b].data();
            }

            if (thread == 0) {
                for (std::size_t b = 0; b < num_blocks; ++b) {
                    tmp_rss_ptrs[b] = output[b].variance;
                }
            } else {
                cur_rss.emplace(sanisizer::cast<I<decltype(cur_rss->size())> >(num_blocks));
                for (std::size_t b = 0; b < num_blocks; ++b) {
                    tatami::resize_container_to_Index_size((*cur_rss)[b], NR);
                    tmp_rss_ptrs[b] = (*cur_rss)[b].data();
                }
            }
        }

        if (is_sparse) {
            auto ext = tatami::consecutive_extractor<true>(mat, false, s, l);
            auto vbuffer = tatami::create_container_of_Index_size<std::vector<Value_> >(NR);
            auto ibuffer = tatami::create_container_of_Index_size<std::vector<Index_> >(NR);

            auto nonzeros = sanisizer::create<std::vector<std::vector<Index_> > >(num_blocks);
            std::vector<quickstats::RssRunningSparse<Index_, Value_, Stat_> > runners;
            sanisizer::reserve(runners, num_blocks);
            for (std::size_t b = 0; b < num_blocks; ++b) {
                tatami::resize_container_to_Index_size(nonzeros[b], NR);
                runners.emplace_back(NR, tmp_mean_ptrs[b], tmp_rss_ptrs[b], nonzeros[b].data());
            }

            for (Index_ x = 0; x < l; ++x) {
                auto out = ext->fetch(vbuffer.data(), ibuffer.data());
                runners[block[x + s]].add(out.number, out.value, out.index);
            }

            for (std::size_t b = 0; b < num_blocks; ++b) {
                runners[b].finish();
            }

        } else {
            auto ext = tatami::consecutive_extractor<false>(mat, false, s, l);
            auto buffer = tatami::create_container_of_Index_size<std::vector<Value_> >(NR);

            std::vector<quickstats::RssRunningDense<Value_, Stat_> > runners;
            sanisizer::reserve(runners, num_blocks);
            for (std::size_t b = 0; b < num_blocks; ++b) {
                runners.emplace_back(NR, tmp_mean_ptrs[b], tmp_rss_ptrs[b]);
            }

            for (Index_ x = 0; x < l; ++x) {
                auto out = ext->fetch(buffer.data());
                runners[block[x + s]].add(out);
            }

            for (std::size_t b = 0; b < num_blocks; ++b) {
                runners[b].finish();
            }
        }

        if (do_parallel) {
            auto cur_count = sanisizer::create<std::vector<Index_> >(num_blocks);
            for (Index_ x = 0; x < l; ++x) {
                cur_count[block[x + s]] += 1;
            }
            (*all_partial_count)[thread] = std::move(cur_count);

            (*all_partial_mean)[thread] = std::move(cur_mean);
            if (thread > 0) {
                (*all_partial_rss)[thread - 1] = std::move(cur_rss);
            }
        }
    }, NC, num_threads);

    if (!do_parallel) {
        for (std::size_t b = 0; b < num_blocks; ++b) {
            quickstats::rss_to_variance(NR, block_sizes[b], output[b].variance);
        }
        if (NC == 0) {
            // If NC == 0, nused = 0 and RssRunning::finish() doesn't even get called.
            // So, we have to fill the output means with NaNs manually.
            for (std::size_t b = 0; b < num_blocks; ++b) {
                std::fill_n(output[b].mean, NR, std::numeric_limits<Stat_>::quiet_NaN());
            }
        }

    } else {
        const auto& ap_mean = *all_partial_mean;
        const auto& ap_rss = *all_partial_rss;

        // Computing the global mean.
        assert(sanisizer::is_equal(num_blocks, block_sizes.size()));
        for (std::size_t b = 0; b < num_blocks; ++b) {
            const auto cur_output = output[b].mean;
            if (block_sizes[b] == 0) {
                std::fill_n(cur_output, NR, std::numeric_limits<Stat_>::quiet_NaN());
                continue;
            }
            for (int u = 0; u < nused; ++u) {
                const auto cur_count = (*((*all_partial_count)[u]))[b];
                if (cur_count == 0) {
                    continue;
                }
                const auto& cur_mean = (*(ap_mean[u]))[b];
                const Stat_ mult = static_cast<Stat_>(cur_count) / static_cast<Stat_>(block_sizes[b]);
                for (Index_ d = 0; d < NR; ++d) {
                    cur_output[d] += cur_mean[d] * mult;
                }
            }
        }

        // Combining the RSS. 
        for (std::size_t b = 0; b < num_blocks; ++b) {
            const auto& cur_global = output[b].mean;
            const auto cur_output = output[b].variance;
            for (int u = 0; u < nused; ++u) {
                const auto cur_count = (*((*all_partial_count)[u]))[b];
                if (cur_count == 0) {
                    continue;
                }
                const auto& cur_mean = (*(ap_mean[u]))[b];
                if (u == 0) { 
                    for (Index_ d = 0; d < NR; ++d) {
                        // At this point, we know cur_count > 0, so we can use the unsafe function.
                        cur_output[d] = quickstats::recenter_rss_unsafe(cur_count, cur_output[d], cur_mean[d], cur_global[d]); 
                    }
                } else {
                    const auto& cur_rss = (*(ap_rss[u - 1]))[b];
                    for (Index_ d = 0; d < NR; ++d) {
                        cur_output[d] += quickstats::recenter_rss_unsafe(cur_count, cur_rss[d], cur_mean[d], cur_global[d]); 
                    }
                }
            }
            quickstats::rss_to_variance(NR, block_sizes[b], cur_output);
        }
    }
}

template<typename Stat_, typename Index_>
void extract_blocked_weights(
    const std::size_t num_blocks,
    const std::vector<Stat_>& block_weights,
    const std::vector<Index_>& block_sizes,
    const Index_ min_size,
    std::vector<Stat_>& tmp_weights
) {
    assert(sanisizer::is_equal(num_blocks, block_weights.size()));
    assert(sanisizer::is_equal(num_blocks, block_sizes.size()));
    tmp_weights.clear();
    for (std::size_t b = 0; b < num_blocks; ++b) {
        if (block_sizes[b] < min_size) { // skip blocks with insufficient cells.
            continue;
        }
        tmp_weights.push_back(block_weights[b]);
    }
}

template<typename Stat_, typename Index_, class Function_>
void extract_blocked_pointers(
    const std::size_t num_blocks,
    const std::vector<ModelGeneVariancesBuffers<Stat_> >& per_block, 
    const std::vector<Index_>& block_sizes,
    const Index_ min_size,
    const Function_ fun,
    std::vector<Stat_*>& tmp_pointers
) {
    assert(sanisizer::is_equal(num_blocks, per_block.size()));
    assert(sanisizer::is_equal(num_blocks, block_sizes.size()));
    tmp_pointers.clear();
    for (std::size_t b = 0; b < num_blocks; ++b) {
        if (block_sizes[b] < min_size) { // skip blocks with insufficient cells.
            continue;
        }
        tmp_pointers.push_back(fun(per_block[b]));
    }
}
/**
 * @endcond
 */

/** 
 * Model the per-feature variances from a log-expression matrix with blocking.
 * The mean and variance of each gene is computed separately for all cells in each block,
 * and a separate trend is fitted to each block to obtain residuals (see `model_gene_variances()`).
 * This ensures that sample and batch effects do not confound the variance estimates.
 *
 * We also compute the average of each statistic across blocks, using the policy described in `ModelGeneVariancesOptions::average_policy`.
 * This is either a quantile (i.e., median, by default) or weighted mean of values for each gene.
 * Weights are determined by `ModelGeneVariancesOptions::block_weight_policy` and are based on the size of each block.
 * The average residual is particularly useful for feature selection with `choose_highly_variable_genes()`.
 *
 * @tparam Value_ Data type of the matrix.
 * @tparam Index_ Integer type of the row/column indices.
 * @tparam Block_ Integer type of the block IDs.
 * @tparam Stat_ Floating-point type of the output statistics.
 *
 * @param mat Matrix of expression values, typically after normalization and log-transformation.
 * Rows should be genes while columns should be cells.
 * @param[in] block Pointer to an array of length equal to the number of cells.
 * Each entry should be a 0-based block identifier in \f$[0, B)\f$ where \f$B\f$ is the total number of blocks.
 * @param num_blocks Total number of blocks, a.k.a., \f$B\f$.
 * @param[out] buffers Collection of pointers of arrays in which to store the output statistics.
 * The length of `ModelGeneVariancesBlockedResults::per_block` should be equal to the number of blocks.
 * @param options Further options.
 */
template<typename Value_, typename Index_, typename Block_, typename Stat_>
void model_gene_variances_blocked(
    const tatami::Matrix<Value_, Index_>& mat, 
    const Block_* const block, 
    const std::size_t num_blocks,
    const ModelGeneVariancesBlockedBuffers<Stat_>& buffers,
    const ModelGeneVariancesOptions& options
) {
    if (!sanisizer::is_equal(num_blocks, buffers.per_block.size())) {
        throw std::runtime_error("length of 'buffers.per_block' is not equal to 'num_blocks'");
    }
    assert(mat.ncol() == 0 || sanisizer::is_less_than(*std::max_element(block, block + mat.ncol()), num_blocks));

    // Just compute the block sizes here for simplicity.
    auto block_sizes = sanisizer::create<std::vector<Index_> >(num_blocks);
    const auto NC = mat.ncol();
    for (Index_ c = 0; c < NC; ++c) {
        block_sizes[block[c]] += 1;
    }

    if (mat.prefer_rows()) {
        model_gene_variances_blocked_direct(mat, block, num_blocks, block_sizes, buffers.per_block, options.num_threads);
    } else {
        model_gene_variances_blocked_running(mat, block, num_blocks, block_sizes, buffers.per_block, options.num_threads);
    }

    FitVarianceTrendWorkspace<Stat_> work;
    auto fopt = options.fit_variance_trend_options;
    fopt.num_threads = options.num_threads;
    bool all_trends_fitted = true;

    const auto NR = mat.nrow();
    for (std::size_t b = 0; b < num_blocks; ++b) {
        const auto& current = buffers.per_block[b];
        if (current.fitted == NULL || current.residual == NULL) {
            all_trends_fitted = false;
            continue;
        }
        if (block_sizes[b] >= 2) {
            fit_variance_trend(NR, current.mean, current.variance, current.fitted, current.residual, work, fopt);
        } else {
            std::fill_n(current.fitted, NR, std::numeric_limits<double>::quiet_NaN());
            std::fill_n(current.residual, NR, std::numeric_limits<double>::quiet_NaN());
        }
    }

    const auto ave_means = buffers.average.mean;
    const auto ave_variances = buffers.average.variance;
    const auto ave_fitted = buffers.average.fitted;
    const auto ave_residuals = buffers.average.residual;

    if ((ave_fitted || ave_residuals) && !all_trends_fitted) {
        throw std::runtime_error("cannot compute average fitted values/residuals without per-block trend fits");
    }

    std::vector<Stat_*> tmp_pointers;
    tmp_pointers.reserve(num_blocks);

    if (options.block_average_policy == BlockAveragePolicy::MEAN) {
        const auto block_weight = scran_blocks::compute_weights<Stat_>(block_sizes, options.block_weight_policy, options.variable_block_weight_parameters);
        std::vector<Stat_> tmp_weights;
        tmp_weights.reserve(num_blocks);

        if (ave_means) {
            extract_blocked_weights(num_blocks, block_weight, block_sizes, static_cast<Index_>(1), tmp_weights);
            extract_blocked_pointers(num_blocks, buffers.per_block, block_sizes, static_cast<Index_>(1), [](const auto& x) -> Stat_* { return x.mean; }, tmp_pointers);
            scran_blocks::parallel_weighted_means(NR, tmp_pointers, tmp_weights.data(), ave_means, /* skip_nan = */ false);
        }

        // Skip blocks without enough cells to compute the variance.
        extract_blocked_weights(num_blocks, block_weight, block_sizes, static_cast<Index_>(2), tmp_weights);

        if (ave_variances) {
            extract_blocked_pointers(num_blocks, buffers.per_block, block_sizes, static_cast<Index_>(2), [](const auto& x) -> Stat_* { return x.variance; }, tmp_pointers);
            scran_blocks::parallel_weighted_means(NR, tmp_pointers, tmp_weights.data(), ave_variances, /* skip_nan = */ false);
        }

        if (ave_fitted) {
            extract_blocked_pointers(num_blocks, buffers.per_block, block_sizes, static_cast<Index_>(2), [](const auto& x) -> Stat_* { return x.fitted; }, tmp_pointers);
            scran_blocks::parallel_weighted_means(NR, tmp_pointers, tmp_weights.data(), ave_fitted, /* skip_nan = */ false);
        }

        if (ave_residuals) {
            extract_blocked_pointers(num_blocks, buffers.per_block, block_sizes, static_cast<Index_>(2), [](const auto& x) -> Stat_* { return x.residual; }, tmp_pointers);
            scran_blocks::parallel_weighted_means(NR, tmp_pointers, tmp_weights.data(), ave_residuals, /* skip_nan = */ false);
        }

    } else if (options.block_average_policy == BlockAveragePolicy::QUANTILE) {
        if (ave_means) {
            extract_blocked_pointers(num_blocks, buffers.per_block, block_sizes, static_cast<Index_>(1), [](const auto& x) -> Stat_* { return x.mean; }, tmp_pointers);
            scran_blocks::parallel_quantiles(NR, tmp_pointers, options.block_quantile, ave_means, /* skip_nan = */ false);
        }

        // Again, skip blocks without enough cells to compute the variance.

        if (ave_variances) {
            extract_blocked_pointers(num_blocks, buffers.per_block, block_sizes, static_cast<Index_>(2), [](const auto& x) -> Stat_* { return x.variance; }, tmp_pointers);
            scran_blocks::parallel_quantiles(NR, tmp_pointers, options.block_quantile, ave_variances, /* skip_nan = */ false);
        }

        if (ave_fitted) {
            extract_blocked_pointers(num_blocks, buffers.per_block, block_sizes, static_cast<Index_>(2), [](const auto& x) -> Stat_* { return x.fitted; }, tmp_pointers);
            scran_blocks::parallel_quantiles(NR, tmp_pointers, options.block_quantile, ave_fitted, /* skip_nan = */ false);
        }

        if (ave_residuals) {
            extract_blocked_pointers(num_blocks, buffers.per_block, block_sizes, static_cast<Index_>(2), [](const auto& x) -> Stat_* { return x.residual; }, tmp_pointers);
            scran_blocks::parallel_quantiles(NR, tmp_pointers, options.block_quantile, ave_residuals, /* skip_nan = */ false);
        }
    }
}

/** 
 * Overload of `model_gene_variances_blocked()` that allocates space for the output statistics.
 *
 * @tparam Stat_ Floating-point type of the output statistics.
 * @tparam Value_ Data type of the matrix.
 * @tparam Index_ Integer type of the row/column indices.
 * @tparam Block_ Integer type of the block IDs.
 *
 * @param mat Matrix of expression values, typically after normalization and log-transformation.
 * Rows should be genes while columns should be cells.
 * @param[in] block Pointer to an array of length equal to the number of cells, containing 0-based block identifiers.
 * @param num_blocks Total number of blocks.
 * @param options Further options.
 *
 * @return Results of the variance modelling in each block.
 * An average for each statistic is also computed if `ModelGeneVariancesOptions::average_policy` is not `BlockAveragePolicy::NONE`.
 */
template<typename Stat_ = double, typename Value_, typename Index_, typename Block_>
ModelGeneVariancesBlockedResults<Stat_> model_gene_variances_blocked(
    const tatami::Matrix<Value_, Index_>& mat,
    const Block_* const block,
    const std::size_t num_blocks,
    const ModelGeneVariancesOptions& options
) {
    const bool do_average = options.compute_average /* for back-compatibility */ && options.block_average_policy != BlockAveragePolicy::NONE;
    ModelGeneVariancesBlockedResults<Stat_> output(
        mat.nrow(), // cast is safe, any tatami Index_ can always fit into a size_t.
        num_blocks,
        do_average,
        options.trend
    );

    ModelGeneVariancesBlockedBuffers<Stat_> buffers;
    sanisizer::resize(buffers.per_block, num_blocks);
    for (std::size_t b = 0; b < num_blocks; ++b) {
        auto& current = buffers.per_block[b];
        current.mean = output.per_block[b].mean.data();
        current.variance = output.per_block[b].variance.data();

        if (options.trend) {
            current.fitted = output.per_block[b].fitted.data();
            current.residual = output.per_block[b].residual.data();
        } else {
            current.fitted = NULL;
            current.residual = NULL;
        }
    }

    if (!do_average) {
        buffers.average.mean = NULL;
        buffers.average.variance = NULL;
        buffers.average.fitted = NULL;
        buffers.average.residual = NULL;
    } else {
        buffers.average.mean = output.average.mean.data();
        buffers.average.variance = output.average.variance.data();

        if (options.trend) {
            buffers.average.fitted = output.average.fitted.data();
            buffers.average.residual = output.average.residual.data();
        } else {
            buffers.average.fitted = NULL;
            buffers.average.residual = NULL;
        }
    }

    model_gene_variances_blocked(mat, block, num_blocks, buffers, options);
    return output;
}

}

#endif
