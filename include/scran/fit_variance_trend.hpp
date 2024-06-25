#ifndef SCRAN_FIT_VARIANCE_TREND_H
#define SCRAN_FIT_VARIANCE_TREND_H

#include <algorithm>
#include <vector>
#include "WeightedLowess/WeightedLowess.hpp"

/**
 * @file fit_variance_trend.hpp
 * @brief Fit a mean-variance trend to log-count data.
 */

namespace scran {

/**
 * @namespace scran::fit_variance_trend
 * @brief Fit a mean-variance trend to log-count data.
 *
 * We fit a trend to the per-feature variances against the means, both of which are computed from log-normalized expression data.
 * We use a LOWESS smoother in several steps:
 *
 * 1. Filter out low-abundance genes, to ensure the span of the smoother is not skewed by many low-abundance genes.
 * 2. Take the quarter-root of the variances, to squeeze the trend towards 1.
 * This makes the trend more "linear" to improve the performance of the LOWESS smoother;
 * it also reduces the chance of obtaining negative fitted values.
 * 3. Apply the LOWESS smoother to the quarter-root variances.
 * This is done using the implementation in the **WeightedLowess** library.
 * 4. Reverse the quarter-root transformation to obtain the fitted values for all non-low-abundance genes.
 * 5. Extrapolate linearly from the left-most fitted value to the origin to obtain fitted values for the previously filtered genes.
 * This is empirically justified by the observation that mean-variance trends of log-expression data are linear at very low abundances.
 */
namespace fit_variance_trend {

/**
 * @brief Parameter defaults for trend fitting.
 */
struct Options {
    /**
     * Minimum mean log-expression for trend fitting.
     * Genes with lower means are not used in trend fitting, and their fitted values are defined by extrapolating the left edge of the fitted trend is extrapolated to the origin.
     * Only used if `Options::mean_filter = true`.
     */
    double minimum_mean = 0.1;

    /**
     * Should any filtering be performed on the mean log-expression of each gene (see `Options::minimum_mean`)?
     * This may need to be disabled if the trend is not being fitted on statistics computed from log-expression values.
     */
    bool mean_filter = true;

    /**
     * Should any transformation of the variances be performed prior to LOWESS smoothing?
     * This may need to be disabled if `FitVarianceTrend` is not being used on statistics computed from log-expression values.
     */
    bool transform = true;

    /**
     * Span for the LOWESS smoother, as a proportion of the total number of points.
     * This is only used if `Options::use_fixed_width = false`.
     */
    double span = 0.3;

    /**
     * Should a fixed-width constraint be applied to the LOWESS smoother?
     * This forces each window to be a minimum width (see `Options::fixed_width`) and avoids problems with large differences in density.
     * For example, the default smoother performs poorly at high abundances where there are few genes.
     */
    bool use_fixed_width = false;

    /**
     * Width of the window to use when `Options::use_fixed_width = true`.
     * This should be relative to the range of `mean` values in `compute()`;
     * the default value is chosen based on the typical range in single-cell RNA-seq data.
     */
    double fixed_width = 1;

    /**
     * Minimum number of observations in each window when `Options::use_fixed_width = true`.
     * This ensures that each window contains at least a given number of observations;
     * if it does not, it is extended using the standard LOWESS logic until the minimum number is achieved.
     */
    int minimum_window_count = 200;
};




/**
 * @brief Fit a mean-variance trend to log-count data.
 *
 * @tparam Float_ Floating-point type for the statistics.
 *
 * @param n Number of features.
 * @param[in] mean Pointer to an array of length `n`, containing the means for all features.
 * @param[in] variance Pointer to an array of length `n`, containing the variances for all features.
 * @param[out] fitted Pointer to an array of length `n`, to store the fitted values.
 * @param[out] residuals Pointer to an array of length `n`, to store the residuals.
 * @param options Further options.
 */
template<typename Float_>
void compute(size_t n, const Float_* mean, const Float_* variance, Float_* fitted, Float_* residuals, const Options& options) {
    std::vector<Float_> xbuffer(n), ybuffer(n);

    auto quad = [](Float_ x) -> Float_ {
        return x * x * x * x;
    };

    size_t counter = 0;
    Float_ min_mean = options.minimum_mean;
    for (size_t i = 0; i < n; ++i) {
        if (!options.mean_filter || mean[i] >= min_mean) {
            xbuffer[counter] = mean[i];
            if (options.transform) {
                ybuffer[counter] = std::pow(variance[i], 0.25); // Using the same quarter-root transform that limma::voom uses.
            } else {
                ybuffer[counter] = variance[i];
            }
            ++counter;
        }
    }

    if (counter < 2) {
        throw std::runtime_error("not enough observations above the minimum mean");
    }


    // Determining the left edge. This needs to be done before
    // SortBy::permute on the xbuffer.
    size_t left_index = std::min_element(xbuffer.begin(), xbuffer.begin() + counter) - xbuffer.begin();
    Float_ left_x = xbuffer[left_index];

    WeightedLowess::SortBy sorter(counter, xbuffer.data());
    std::vector<uint8_t> work;
    sorter.permute(xbuffer.data(), work);
    sorter.permute(ybuffer.data(), work);

    WeightedLowess::Options<Float_> smooth_opt;
    if (options.use_fixed_width) {
        smooth_opt.span = options.minimum_window_count;
        smooth_opt.span_as_proportion = false;
        smooth_opt.minimum_width = options.fixed_width;
    } else {
        smooth_opt.span = options.span;
    }

    std::vector<Float_> fbuffer(counter), rbuffer(counter);
    WeightedLowess::compute(counter, xbuffer.data(), ybuffer.data(), fbuffer.data(), rbuffer.data(), smooth_opt);

    sorter.unpermute(rbuffer.data(), work);
    sorter.unpermute(fbuffer.data(), work);

    // Identifying the left-most fitted value.
    Float_ left_fitted = (options.transform ? quad(fbuffer[left_index]) : fbuffer[left_index]);

    counter = 0;
    for (size_t i = 0; i < n; ++i) {
        if (!options.mean_filter || mean[i] >= min_mean) {
            fitted[i] = (options.transform ? quad(fbuffer[counter]) : fbuffer[counter]);
            ++counter;
        } else {
            fitted[i] = mean[i] / left_x * left_fitted; // draw a y = x line to the origin from the left of the fitted trend.
        }
        residuals[i] = variance[i] - fitted[i];
    }
    return;
}

/**
 * @brief Results of `fit_variance_trend::compute()`.
 *
 * Meaningful instances of this object should generally be constructed by calling the `fit_variance_trend::compute()` function.
 * Empty instances can be default-constructed as placeholders.
 *
 * @tparam Float_ Floating-point type for the statistics.
 */
template<typename Float_>
struct Results {
    /**
     * @cond
     */
    Results() {}

    Results(size_t n) : fitted(n), residuals(n) {}
    /**
     * @endcond
     */

    /**
     * Vector of length equal to the number of features, containing fitted values from the trend.
     */
    std::vector<Float_> fitted;

    /**
     * Vector of length equal to the number of features, containing residuals from the trend.
     */
    std::vector<Float_> residuals;
};

/**
 * Overload of `fit_variance_trend::compute()` that allocates the output vectors.
 *
 * @tparam Float_ Floating-point type for the statistics.
 *
 * @param n Number of features.
 * @param[in] mean Pointer to an array of length `n`, containing the means for all features.
 * @param[in] variance Pointer to an array of length `n`, containing the variances for all features.
 * @param options Further options.
 * 
 * @return Result of the trend fit, containing the fitted values and residuals for each gene. 
 */
template<typename Float_>
Results<Float_> compute(size_t n, const Float_* mean, const Float_* variance, const Options& options) {
    Results<Float_> output(n);
    compute(n, mean, variance, output.fitted.data(), output.residuals.data(), options);
    return output;
}

}

}

#endif
