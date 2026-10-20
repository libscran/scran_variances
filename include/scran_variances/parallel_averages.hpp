#ifndef SCRAN_VARIANCES_PARALLEL_AVERAGES_HPP
#define SCRAN_VARIANCES_PARALLEL_AVERAGES_HPP

#include <vector>
#include <cstddef>
#include <cassert>
#include <algorithm>
#include <cmath>
#include <limits>

#include "quickstats/quickstats.hpp"
#include "utils.hpp"

namespace scran_variances {

template<typename Index_, typename Block_, class Function_, typename Weight_, typename Output_>
void parallel_weighted_means(
    const Index_ n,
    const Block_ nblocks,
    const Function_ fun,
    const std::vector<Weight_>& block_weights,
    const std::vector<Index_>& block_sizes,
    const Index_ min_size,
    Output_* const out
) {
    assert(sanisizer::is_equal(nblocks, block_weights.size()));
    assert(sanisizer::is_equal(nblocks, block_sizes.size()));
    std::fill_n(out, n, 0);

    Weight_ total = 0;
    for (I<decltype(nblocks)> b = 0; b < nblocks; ++b) {
        // By removing everything below 'min_size', all values in each entry of 'in' should be non-NaN.
        if (block_sizes[b] < min_size) {
            continue;
        }

        // Removing zero weight so that we don't have to worry about '0 * Inf ==> NaN'.
        const Weight_ weight = block_weights[b];
        if (weight == 0) {
            continue;
        }

        auto current = fun(b);
        if (weight == 1) { 
            for (I<decltype(n)> i = 0; i < n; ++i) {
                assert(!std::isnan(current[i]));
                out[i] += current[i];
            }
        } else {
            for (I<decltype(n)> i = 0; i < n; ++i) {
                assert(!std::isnan(current[i]));
                out[i] += current[i] * weight;
            }
        }
        total += weight;
    }

    if (total) {
        for (I<decltype(n)> i = 0; i < n; ++i) {
            out[i] /= total;
        }
    } else {
        std::fill_n(out, n, std::numeric_limits<Output_>::quiet_NaN());
    }
}

template<typename Function_, typename Index_, typename Output_>
void parallel_quantiles(
    const std::size_t n,
    const std::size_t nblocks,
    const Function_ fun,
    const std::vector<Index_>& block_sizes,
    const Index_ min_size,
    const double quantile,
    Output_* const out
) {
    assert(sanisizer::is_equal(nblocks, block_sizes.size()));

    typedef I<decltype(fun(0)[0])> Stat_;
    std::vector<const Stat_*> tmp_ptrs;
    tmp_ptrs.reserve(nblocks);
    for (I<decltype(nblocks)> b = 0; b < nblocks; ++b) {
        if (block_sizes[b] >= min_size) { // skip blocks with insufficient cells, as these will have all-NaN stats.
            tmp_ptrs.push_back(fun(b));
        }
    }

    const auto remaining_nblocks = tmp_ptrs.size();
    std::vector<Stat_> tmp_buffer;
    tmp_buffer.reserve(remaining_nblocks);
    quickstats::SingleQuantileFixedNumber<Output_> calc((sanisizer::Cast(remaining_nblocks)), quantile);
    for (I<decltype(n)> g = 0; g < n; ++g) {
        tmp_buffer.clear();
        for (const auto ptr : tmp_ptrs) {
            assert(!std::isnan(ptr[g])); // there should be no NaNs anymore once all blocks with insufficient cells have been removed.
            tmp_buffer.push_back(ptr[g]);
        }
        out[g] = calc(tmp_buffer.data());
    }
}

}

#endif
