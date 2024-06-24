#ifndef SCRAN_CHOOSE_HIGHLY_VARIABLE_GENES_HPP
#define SCRAN_CHOOSE_HIGHLY_VARIABLE_GENES_HPP

#include <vector>
#include <algorithm>
#include <numeric>
#include <cstdint>

/**
 * @file choose_highly_variable_genes.hpp
 * @brief Choose highly variable genes for downstream analyses.
 */

namespace scran {

/**
 * @namespace choose_highly_variable_genes
 * @brief Choose highly variable genes for downstream analyses.
 */
namespace choose_highly_variable_genes {

struct Options {
    /**
     * Number of top genes to choose.
     * This should be positive.
     */
    size_t top = 4000;

    /**
     * Whether larger statistics correspond to higher variances.
     */
    bool larger = true;

    /**
     * Whether to keep all genes with statistics that are tied with the `top`-th gene.
     * If `false`, ties are arbitrarily broken but the number of retained genes will not be greater than `top`.
     */
    bool keep_ties = true;
};

/**
 * @tparam Stat_ Type of the variance statistic.
 * @tparam Bool_ Type to be used as a boolean.
 *
 * @param n Number of genes.
 * @param[in] statistic Pointer to an array of length `n` containing the per-gene variance statistics.
 * @param[out] output Pointer to an array of length `n`. 
 * On output, this is filled with `true` if the gene is to be retained and `false` otherwise.
 * @param options Further options.
 */
template<typename Stat_, typename Bool_>
void compute(size_t n, const Stat_* statistic, Bool_* output, const Options& options) {
    if (options.top >= n) {
        std::fill_n(output, n, true);
        return;
    } else if (options.top == 0) {
        std::fill_n(output, n, false);
        return;
    }

    std::vector<size_t> collected(n);
    std::iota(collected.begin(), collected.end(), static_cast<size_t>(0));
    auto cBegin = collected.begin(), cMid = cBegin + options.top - 1, cEnd = collected.end();
    if (options.larger) {
        std::nth_element(cBegin, cMid, cEnd, [&](size_t l, size_t r) -> bool { return statistic[l] > statistic[r]; });
    } else {
        std::nth_element(cBegin, cMid, cEnd, [&](size_t l, size_t r) -> bool { return statistic[l] < statistic[r]; });
    }

    if (!options.keep_ties) {
        std::fill_n(output, n, false);
        for (size_t i = 0; i < options.top; ++i) {
            output[collected[i]] = true; 
        }
    } else {
        auto threshold = statistic[collected[options.top - 1]];
        if (options.larger) {
            for (size_t i = 0; i < n; ++i) {
                output[i] = statistic[i] >= threshold; 
            }
        } else {
            for (size_t i = 0; i < n; ++i) {
                output[i] = statistic[i] <= threshold; 
            }
        }
    }
}

/**
 * @tparam Stat_ Type of the variance statistic.
 * @tparam Bool_ Type to be used as a boolean.
 *
 * @param n Number of genes.
 * @param[in] statistic Pointer to an array of length `n` containing the per-gene variance statistics.
 * @param options Further options.
 *
 * @return A vector of booleans of length `n`, indicating whether each gene is to be retained.
 */
template<typename Bool_ = uint8_t, typename Stat_>
std::vector<Bool_> compute(size_t n, const Stat_* statistic, const Options& options) {
    std::vector<Bool_> output(n);
    compute(n, statistic, output.data(), options);
    return output;
}

}

}
#endif
