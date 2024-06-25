#ifndef SCRAN_MODEL_GENE_VARIANCES_H
#define SCRAN_MODEL_GENE_VARIANCES_H

#include "tatami/tatami.hpp"
#include "tatami_stats/tatami_stats.hpp"

#include "block_weights.hpp"
#include "average_vectors.hpp"
#include "fit_variance_trend.hpp"

#include <algorithm>
#include <vector>
#include <limits>

/**
 * @file model_gene_variances.hpp
 * @brief Model the per-gene variances. 
 *
 * This scans through a log-transformed normalized expression matrix and computes per-feature means and variances.
 * It then fits a trend to the variances with respect to the means using `fit_variance_trend::compute()`.
 * We assume that most genes at any given abundance are not highly variable, such that the fitted value of the trend is interpreted as the "uninteresting" variance - 
 * this is mostly attributed to technical variation like sequencing noise, but can also represent constitutive biological noise like transcriptional bursting.
 * Under this assumption, the residual can be treated as a quantification of biologically interesting variation, and can be used to identify relevant features for downstream analyses.
 */

namespace scran {

/**
 * @file scran::model_gene_variances
 * @brief Model the per-gene variances. 
 */
namespace model_gene_variances {

/**
 * @brief Default parameters for variance modelling.
 */
struct Options {
    /**
     * Options for fitting the mean-variance trend.
     */
    fit_variance_trend::Options fit_variance_trend_options;

    /**
     * Weighting policy to use for averaging statistics across blocks.
     * Only relevant for `compute_blocked()` overloads where averaged outputs are requested.
     */
    block_weights::Policy block_weight_policy = block_weights::Policy::VARIABLE;

    /**
     * Parameters for the variable block weights.
     * Only relevant for `compute_blocked()` overloads where averaged outputs are requested and `Options::block_weight_policy = block_weights::Policy::VARIABLE`.
     */
    block_weights::VariableParameters variable_block_weight_parameters; 

    /**
     * Whether to compute the average of each statistic across blocks.
     * Note that this only affects the `compute_blocked()` method that returns a `BlockResults` object.
     */
    bool compute_average = true;

    /**
     * Number of threads to use. 
     */
    int num_threads = 1;
};

/**
 * @cond
 */
namespace internal {

template<typename Value_, typename Index_, typename Stat_, typename Block_> 
void compute_dense_row(const tatami::Matrix<Value_, Index_>* mat, std::vector<Stat_*>& means, std::vector<Stat_*>& variances, const Block_* block, const std::vector<Index_>& block_size, int num_threads) {
    bool blocked = (block != NULL);
    auto nblocks = block_size.size();
    auto NR = mat->nrow(), NC = mat->ncol();

    tatami::parallelize([&](size_t, Index_ start, Index_ length) -> void {
        std::vector<Stat_> tmp_means(blocked ? nblocks : 0);
        std::vector<Stat_> tmp_vars(blocked ? nblocks : 0);

        std::vector<Value_> buffer(NC);
        auto ext = tatami::consecutive_extractor<false>(mat, true, start, length);
        for (Index_ r = start, end = start + length; r < end; ++r) {
            auto ptr = ext->fetch(buffer.data());

            if (blocked) {
                tatami_stats::grouped_variances::direct(ptr, NC, block, nblocks, block_size.data(), tmp_means.data(), tmp_vars.data(), false, static_cast<Index_*>(NULL));
                for (size_t b = 0; b < nblocks; ++b) {
                    means[b][r] = tmp_means[b];
                    variances[b][r] = tmp_vars[b];
                }
            } else {
                auto stat = tatami_stats::variances::direct(ptr, NC, false);
                means[0][r] = stat.first;
                variances[0][r] = stat.second;
            }
        }
    }, NR, num_threads);
}

template<typename Value_, typename Index_, typename Stat_, typename Block_> 
void compute_sparse_row(const tatami::Matrix<Value_, Index_>* mat, std::vector<Stat_*>& means, std::vector<Stat_*>& variances, const Block_* block, const std::vector<Index_>& block_size, int num_threads) {
    bool blocked = (block != NULL);
    auto nblocks = block_size.size();
    auto NR = mat->nrow(), NC = mat->ncol();

    tatami::parallelize([&](size_t, Index_ start, Index_ length) -> void {
        std::vector<Stat_> tmp_means(nblocks);
        std::vector<Stat_> tmp_vars(nblocks);
        std::vector<Index_> tmp_nzero(nblocks);

        std::vector<Value_> vbuffer(NC);
        std::vector<Index_> ibuffer(NC);
        tatami::Options opt;
        opt.sparse_ordered_index = false;
        auto ext = tatami::consecutive_extractor<true>(mat, true, start, length, opt);

        for (Index_ r = start, end = start + length; r < end; ++r) {
            auto range = ext->fetch(vbuffer.data(), ibuffer.data());

            if (blocked) {
                tatami_stats::grouped_variances::direct(range.value, range.index, range.number, block, nblocks, block_size.data(), tmp_means.data(), tmp_vars.data(), tmp_nzero.data(), false, static_cast<Index_*>(NULL));
                for (size_t b = 0; b < nblocks; ++b) {
                    means[b][r] = tmp_means[b];
                    variances[b][r] = tmp_vars[b];
                }
            } else {
                auto stat = tatami_stats::variances::direct(range.value, range.number, NC, false);
                means[0][r] = stat.first;
                variances[0][r] = stat.second;
            }
        }
    }, NR, num_threads);
}

template<typename Value_, typename Index_, typename Stat_, typename Block_> 
void compute_dense_column(const tatami::Matrix<Value_, Index_>* mat, std::vector<Stat_*>& means, std::vector<Stat_*>& variances, const Block_* block, const std::vector<Index_>& block_size, int num_threads) {
    bool blocked = (block != NULL);
    auto nblocks = block_size.size();
    auto NR = mat->nrow(), NC = mat->ncol();

    tatami::parallelize([&](size_t thread, Index_ start, Index_ length) -> void {
        std::vector<Value_> buffer(length);
        auto ext = tatami::consecutive_extractor<false>(mat, false, 0, NC, start, length);

        std::vector<tatami_stats::LocalOutputBuffer<Stat_> > local_var_output;
        local_var_output.reserve(nblocks);
        std::vector<tatami_stats::LocalOutputBuffer<Stat_> > local_mean_output;
        local_mean_output.reserve(nblocks);
        std::vector<tatami_stats::variances::RunningDense<Stat_, Value_, Index_> > runners;
        runners.reserve(nblocks);

        for (size_t b = 0; b < nblocks; ++b) {
            local_var_output.emplace_back(thread, start, length, variances[b]);
            local_mean_output.emplace_back(thread, start, length, means[b]);
            runners.emplace_back(length, local_mean_output.back().data(), local_var_output.back().data(), false);
        }

        if (blocked) {
            for (Index_ c = 0; c < NC; ++c) {
                auto ptr = ext->fetch(buffer.data());
                runners[block[c]].add(ptr);
            }
        } else {
            for (Index_ c = 0; c < NC; ++c) {
                auto ptr = ext->fetch(buffer.data());
                runners[0].add(ptr);
            }
        }

        for (size_t b = 0; b < nblocks; ++b) {
            runners[b].finish();
            local_var_output[b].transfer();
            local_mean_output[b].transfer();
        }
    }, NR, num_threads);
}

template<typename Value_, typename Index_, typename Stat_, typename Block_> 
void compute_sparse_column(const tatami::Matrix<Value_, Index_>* mat, std::vector<Stat_*>& means, std::vector<Stat_*>& variances, const Block_* block, const std::vector<Index_>& block_size, int num_threads) {
    bool blocked = (block != NULL);
    auto nblocks = block_size.size();
    auto NR = mat->nrow(), NC = mat->ncol();
    std::vector<std::vector<Index_> > nonzeros(nblocks, std::vector<Index_>(NR));

    tatami::parallelize([&](size_t thread, Index_ start, Index_ length) -> void {
        std::vector<Value_> vbuffer(length);
        std::vector<Index_> ibuffer(length);
        tatami::Options opt;
        opt.sparse_ordered_index = false;
        auto ext = tatami::consecutive_extractor<true>(mat, false, 0, NC, start, length, opt);

        std::vector<tatami_stats::LocalOutputBuffer<Stat_> > local_var_output;
        local_var_output.reserve(nblocks);
        std::vector<tatami_stats::LocalOutputBuffer<Stat_> > local_mean_output;
        local_mean_output.reserve(nblocks);
        std::vector<tatami_stats::variances::RunningSparse<Stat_, Value_, Index_> > runners;
        runners.reserve(nblocks);

        for (size_t b = 0; b < nblocks; ++b) {
            local_var_output.emplace_back(thread, start, length, variances[b]);
            local_mean_output.emplace_back(thread, start, length, means[b]);
            runners.emplace_back(length, local_mean_output.back().data(), local_var_output.back().data(), false, start);
        }

        if (blocked) {
            for (Index_ c = 0; c < NC; ++c) {
                auto range = ext->fetch(vbuffer.data(), ibuffer.data());
                runners[block[c]].add(range.value, range.index, range.number);
            }
        } else {
            for (Index_ c = 0; c < NC; ++c) {
                auto range = ext->fetch(vbuffer.data(), ibuffer.data());
                runners[0].add(range.value, range.index, range.number);
            }
        }

        for (size_t b = 0; b < nblocks; ++b) {
            runners[b].finish();
            local_var_output[b].transfer();
            local_mean_output[b].transfer();
        }
    }, NR, num_threads);
}

template<typename Value_, typename Index_, typename Stat_, typename Block_> 
void compute(const tatami::Matrix<Value_, Index_>* mat, std::vector<Stat_*>& means, std::vector<Stat_*>& variances, const Block_* block, const std::vector<Index_>& block_size, int num_threads) {
    if (mat->prefer_rows()) {
        if (mat->sparse()) {
            compute_sparse_row(mat, means, variances, block, block_size, num_threads);
        } else {
            compute_dense_row(mat, means, variances, block, block_size, num_threads);
        }
    } else {
        if (mat->sparse()) {
            compute_sparse_column(mat, means, variances, block, block_size, num_threads);
        } else {
            compute_dense_column(mat, means, variances, block, block_size, num_threads);
        }
    }
}

}
/**
 * @endcond
 */

/** 
 * Compute and model the per-feature variances from a log-expression matrix with blocking.
 * The mean and variance of each gene is computed separately for all cells in each block, and a separate trend is fitted to each block to obtain residuals.
 * This ensures that sample and batch effects do not confound the variance estimates.
 *
 * We also compute the average of each statistic across blocks, using the weighting strategy specified in `Options::block_weight_policy`.
 * The average residual is particularly useful for feature selection with `choose_highly_variable_genes::compute()`.
 *
 * @tparam Value_ Data type of the matrix.
 * @tparam Index_ Integer type for the row/column indices.
 * @tparam Block_ Integer type to hold the block IDs.
 * @tparam Stat_ Floating-point type for the output statistics.
 *
 * @param mat Pointer to a feature-by-cells **tatami** matrix containing log-expression values.
 * @param[in] block Pointer to an array of length equal to the number of cells.
 * Each entry should be a 0-based block identifier in \f$[0, B)\f$ where \f$B\f$ is the total number of blocks.
 * `block` can also be a `nullptr`, in which case all cells are assumed to belong to the same block.
 * @param[out] means Vector of length equal to the number of blocks, containing pointers to output arrays of length equal to the number of rows in `mat`.
 * Each vector stores the mean of each feature in the corresponding block of cells.
 * @param[out] variances Vector of length equal to the number of blocks, containing pointers to output arrays of length equal to the number of rows in `mat`.
 * Each vector stores the variance of each feature in the corresponding block of cells.
 * @param[out] fitted Vector of length equal to the number of blocks, containing pointers to output arrays of length equal to the number of rows in `mat`.
 * Each vector stores the fitted value of the trend for each feature in the corresponding block of cells.
 * @param[out] residuals Vector of length equal to the number of blocks, containing pointers to output arrays of length equal to the number of rows in `mat`.
 * Each vector stores the residual from the trend for each feature in the corresponding block of cells.
 * @param[out] ave_means Pointer to an array of length equal to the number of rows in `mat`, storing the average mean across blocks for each gene.
 * If `nullptr`, the average calculation is skipped.
 * @param[out] ave_variances Pointer to an array of length equal to the number of rows in `mat`, storing the average variance across blocks for each gene.
 * If `nullptr`, the average calculation is skipped.
 * @param[out] ave_fitted Pointer to an array of length equal to the number of rows in `mat`, storing the average fitted value across blocks for each gene.
 * If `nullptr`, the average calculation is skipped.
 * @param[out] ave_residuals Pointer to an array of length equal to the number of rows in `mat`, storing the average residual across blocks for each gene.
 * If `nullptr`, the average calculation is skipped.
 * @param options Further options.
 */
template<typename Value_, typename Index_, typename Block_, typename Stat_>
void compute_blocked(
    const tatami::Matrix<Value_, Index_>* mat, 
    const Block_* block, 
    std::vector<Stat_*> means, 
    std::vector<Stat_*> variances,
    std::vector<Stat_*> fitted, 
    std::vector<Stat_*> residuals,
    Stat_* ave_means,
    Stat_* ave_variances,
    Stat_* ave_fitted,
    Stat_* ave_residuals,
    const Options& options)
{
    Index_ NR = mat->nrow(), NC = mat->ncol();
    std::vector<Index_> block_size;

    if (block) {
        block_size = tatami_stats::tabulate_groups(block, NC);
        internal::compute(mat, means, variances, block, block_size, options.num_threads);
    } else {
        block_size.push_back(NC); // everything is one big block.
        internal::compute(mat, means, variances, block, block_size, options.num_threads);
    }

    auto fopt = options.fit_variance_trend_options;
    for (size_t b = 0, nblocks = block_size.size(); b < nblocks; ++b) {
        if (block_size[b] >= 2) {
            fit_variance_trend::compute(NR, means[b], variances[b], fitted[b], residuals[b], fopt);
        } else {
            std::fill(fitted[b], fitted[b] + NR, std::numeric_limits<double>::quiet_NaN());
            std::fill(residuals[b], residuals[b] + NR, std::numeric_limits<double>::quiet_NaN());
        }
    }

    if (ave_means || ave_variances || ave_fitted || ave_residuals) {
        std::vector<double> block_weight = block_weights::compute(block_size, options.block_weight_policy, options.variable_block_weight_parameters);

        // Check whether there might be NaNs for very small blocks, 
        // in which case we need to ignore them.
        Index_ min_block = NC;
        for (size_t b = 0, nblocks = block_size.size(); b < nblocks; ++b) {
            if (block_size[b] < min_block) {
                min_block = block_size[b];
            }
        }

        if (ave_means) {
            average_vectors::compute_weighted(NR, means, block_weight.data(), ave_means, min_block < 1);
        }
        if (ave_variances) {
            average_vectors::compute_weighted(NR, variances, block_weight.data(), ave_variances, min_block < 2);
        }
        if (ave_fitted) {
            average_vectors::compute_weighted(NR, fitted, block_weight.data(), ave_fitted, min_block < 2);
        }
        if (ave_residuals) {
            average_vectors::compute_weighted(NR, residuals, block_weight.data(), ave_residuals, min_block < 2);
        }
    }
}

/** 
 * This overload of `compute_blocked()` omits the calculation of the averaged statistics across blocks.
 *
 * @tparam Value_ Data type of the matrix.
 * @tparam Index_ Integer type for the row/column indices.
 * @tparam Block_ Integer type to hold the block IDs.
 * @tparam Stat_ Floating-point type for the output statistics.
 *
 * @param mat Pointer to a feature-by-cells **tatami** matrix containing log-expression values.
 * @param[in] block Pointer to an array of length equal to the number of cells, containing 0-based block identifiers.
 * This can also be a `nullptr`, in which case all cells are assumed to belong to the same block.
 * @param[out] means Vector of length equal to the number of blocks, containing pointers to output arrays of length equal to the number of rows in `mat`.
 * Each vector stores the mean of each feature in the corresponding block of cells.
 * @param[out] variances Vector of length equal to the number of blocks, containing pointers to output arrays of length equal to the number of rows in `mat`.
 * Each vector stores the variance of each feature in the corresponding block of cells.
 * @param[out] fitted Vector of length equal to the number of blocks, containing pointers to output arrays of length equal to the number of rows in `mat`.
 * Each vector stores the fitted value of the trend for each feature in the corresponding block of cells.
 * @param[out] residuals Vector of length equal to the number of blocks, containing pointers to output arrays of length equal to the number of rows in `mat`.
 * Each vector stores the residual from the trend for each feature in the corresponding block of cells.
 */
template<typename Value_, typename Index_, typename Block_, typename Stat_>
void compute_blocked(
    const tatami::Matrix<Value_, Index_>* mat, 
    const Block_* block, 
    std::vector<Stat_*> means, 
    std::vector<Stat_*> variances,
    std::vector<Stat_*> fitted, 
    std::vector<Stat_*> residuals,
    const Options& options)
{
    compute_blocked(
        mat, 
        block, 
        std::move(means), 
        std::move(variances), 
        std::move(fitted), 
        std::move(residuals),
        static_cast<Stat_*>(NULL),
        static_cast<Stat_*>(NULL),
        static_cast<Stat_*>(NULL),
        static_cast<Stat_*>(NULL),
        options
    );
}

/** 
 * Compute and model the per-feature variances from a log-expression matrix, where all cells are assumed to belong to a single block.
 * This returns the mean and variance for each feature, as well as the fitted value and residuals from the mean-variance trend fitted across features.
 *
 * @tparam Value_ Data type of the matrix.
 * @tparam Index_ Integer type for the row/column indices.
 * @tparam Stat_ Floating-point type for the output statistics.
 *
 * @param mat Pointer to a feature-by-cells **tatami** matrix containing log-expression values.
 * @param[out] means Pointer to an output array of length equal to the number of rows in `mat`, used to store the mean of each feature.
 * @param[out] variances Pointer to an output array of length equal to the number of rows in `mat`, used to store the variance of each feature.
 * @param[out] fitted Pointer to an output array of length equal to the number of rows in `mat`, used to store the fitted value of the trend.
 * @param[out] residuals Pointer to an output array of length equal to the number of rows in `mat`, used to store the residual from the trend for each feature.
 */
template<typename Value_, typename Index_, typename Stat_> 
void compute(const tatami::Matrix<Value_, Index_>* mat, Stat_* means, Stat_* variances, Stat_* fitted, Stat_* residuals, const Options& options) {
    compute_blocked(
        mat, 
        static_cast<Index_*>(NULL), 
        std::vector<Stat_*>{ means }, 
        std::vector<Stat_*>{ variances }, 
        std::vector<Stat_*>{ fitted }, 
        std::vector<Stat_*>{ residuals },
        options
    );
}

/**
 * @brief Results of variance modelling without blocks.
 *
 * @tparam Stat_ Floating-point type for the output statistics.
 *
 * Meaningful instances of this object should generally be constructed by calling `model_gene_variances::compute()`.
 * Empty instances can be default-constructed as placeholders.
 */
template<typename Stat_>
struct Results {
    /**
     * @cond
     */
    Results() = default;

    Results(size_t ngenes) : means(ngenes), variances(ngenes), fitted(ngenes), residuals(ngenes) {}
    /**
     * @endcond
     */

    /**
     * Vector of length equal to the number of genes, containing the mean log-expression for each gene.
     */
    std::vector<Stat_> means;

    /**
     * Vector of length equal to the number of genes, containing the variance in the log-expression for each gene.
     */
    std::vector<Stat_> variances;

    /**
     * Vector of length equal to the number of genes, containing the fitted value of the mean-variance trend for each gene.
     */
    std::vector<Stat_> fitted;

    /**
     * Vector of length equal to the number of genes, containing the residuals of the mean-variance trend for each gene.
     */
    std::vector<Stat_> residuals;
};

/** 
 * Overload of `compute()` that allocates space for the output statistics.
 *
 * @tparam Value_ Data type of the matrix.
 * @tparam Index_ Integer type for the row/column indices.
 * @tparam Stat_ Floating-point type for the output statistics.
 *
 * @param mat Pointer to a feature-by-cells **tatami** matrix containing log-expression values.
 * @param options Further options.
 *
 * @return Results of the variance modelling.
 */
template<typename Stat_ = double, typename Value_ = double, typename Index_ = int>
Results<Stat_> compute(const tatami::Matrix<Value_, Index_>* mat, const Options& options) {
    Results<Stat_> output(mat->nrow());
    compute(
        mat, 
        output.means.data(), 
        output.variances.data(), 
        output.fitted.data(), 
        output.residuals.data(),
        options
    );
    return output;
}

/**
 * @brief Results of variance modelling with blocks.
 *
 * @tparam Stat_ Floating-point type for the output statistics.
 *
 * Meaningful instances of this object should generally be constructed by calling the `ModelGeneVariances::run_blocked()` method.
 * Empty instances can be default-constructed as placeholders.
 */
template<typename Stat_>
struct BlockResults {
    /**
     * @cond
     */
    BlockResults() = default;

    BlockResults(size_t ngenes, int nblocks, bool compute_average) : 
        per_block(nblocks, Results<Stat_>(ngenes)),
        average(compute_average ? ngenes : 0) {}
    /**
     * @endcond
     */

    /**
     * Vector of length equal to the number of blocks, where each entry contains the variance modelling results for a single block.
     */
    std::vector<Results<Stat_> > per_block;

    /**
     * Average across blocks for all statistics in `per_block`.
     * This is only populated if `Options::compute_average = true`.
     */
    Results<Stat_> average;
};

/** 
 * Overload of `compute_blocked()` that allocates space for the output statistics.
 *
 * @tparam Stat_ Floating-point type for the output statistics.
 * @tparam Value_ Data type of the matrix.
 * @tparam Index_ Integer type for the row/column indices.
 * @tparam Block_ Integer type, to hold the block IDs.
 *
 * @param mat Pointer to a feature-by-cells **tatami** matrix containing log-expression values.
 * @param[in] block Pointer to an array of length equal to the number of cells, containing 0-based block identifiers.
 * This may also be a `nullptr` in which case all cells are assumed to belong to the same block.
 * @param options Further options.
 *
 * @return Results of the variance modelling in each block.
 * An average for each statistic is also computed if `Options::compute_average = true.
 */
template<typename Stat_ = double, typename Value_ = double, typename Index_ = int, typename Block_ = int>
BlockResults<Stat_> compute_blocked(const tatami::Matrix<Value_, Index_>* mat, const Block_* block, const Options& options) {
    size_t nblocks = (block ? tatami_stats::total_groups(block, mat->ncol()) : 1);
    BlockResults<Stat_> output(mat->nrow(), nblocks, options.compute_average);

    std::vector<double*> mean_ptr, var_ptr, fit_ptr, resid_ptr;
    mean_ptr.reserve(nblocks);
    var_ptr.reserve(nblocks);
    fit_ptr.reserve(nblocks);
    resid_ptr.reserve(nblocks);

    for (size_t b = 0; b < nblocks; ++b) {
        mean_ptr.push_back(output.per_block[b].means.data());
        var_ptr.push_back(output.per_block[b].variances.data());
        fit_ptr.push_back(output.per_block[b].fitted.data());
        resid_ptr.push_back(output.per_block[b].residuals.data());
    }

    if (options.compute_average) {
        compute_blocked(
            mat, 
            block, 
            std::move(mean_ptr), 
            std::move(var_ptr), 
            std::move(fit_ptr), 
            std::move(resid_ptr),
            output.average.means.data(), 
            output.average.variances.data(), 
            output.average.fitted.data(), 
            output.average.residuals.data(),
            options
        );
    } else {
        compute_blocked(
            mat, 
            block, 
            std::move(mean_ptr), 
            std::move(var_ptr), 
            std::move(fit_ptr), 
            std::move(resid_ptr),
            options
        );
    }

    return output;
}

}

}

#endif
