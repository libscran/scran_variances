#include <gtest/gtest.h>
#include <vector>

#include "scran_variances/parallel_averages.hpp"
#include "scran_tests/scran_tests.hpp"

TEST(ParallelWeightedMeans, Basic) {
    std::vector<std::vector<double> > stuff {
        std::vector<double>{1, 2, 3, 4, 5},
        std::vector<double>{2, 1, 2, 6, 8},
        std::vector<double>{3, 7, 3, 3, 2}
    };

    std::vector<double> weights { 0.5, 2, 1.5 };
    std::vector<int> sizes { 2, 1, 3 };

    {
        std::vector<double> output(5);
        scran_variances::parallel_weighted_means(
            5, 
            3,
            [&](int b) -> const double* { return stuff[b].data(); },
            weights, 
            sizes,
            1,
            output.data()
        );
        std::vector<double> expected{ 2.250, 3.375, 2.500, 4.625, 5.375 };
        scran_tests::compare_almost_equal_containers(expected, output, {});
    }

    // Skipping blocks below the minimum size.
    {
        std::vector<double> output(5);
        scran_variances::parallel_weighted_means(
            5,
            3,
            [&](int b) -> const double* { return stuff[b].data(); },
            weights, 
            sizes,
            2,
            output.data()
        );
        std::vector<double> expected{ 2.5, 5.75, 3.0, 3.25, 2.75 };
        scran_tests::compare_almost_equal_containers(expected, output, {});
    }

    // Skipping blocks with a weight of zero.
    {
        std::vector<double> weights { 0, 3, 1 };
        std::vector<double> output(5);
        scran_variances::parallel_weighted_means(
            5,
            3,
            [&](int b) -> const double* { return stuff[b].data(); },
            weights, 
            sizes,
            1,
            output.data()
        );
        std::vector<double> expected{ 2.25, 2.5, 2.25, 5.25, 6.5 };
        scran_tests::compare_almost_equal_containers(expected, output, {});
    }

    // All blocks skipped.
    {
        std::vector<double> weights { 0, 1, 0 };
        std::vector<double> output(5);
        scran_variances::parallel_weighted_means(
            5, 
            3,
            [&](int b) -> const double* { return stuff[b].data(); },
            weights, 
            sizes,
            2,
            output.data()
        );
        std::vector<double> expected(5, std::numeric_limits<double>::quiet_NaN());
        scran_tests::compare_almost_equal_containers(expected, output, {});
    }
}

TEST(ParallelQuantiles, Basic) {
    std::vector<std::vector<double> > stuff {
        std::vector<double>{1, 2, 3, 4, 5},
        std::vector<double>{2, 1, 2, 6, 8},
        std::vector<double>{3, 7, 3, 3, 2}
    };

    {
        std::vector<double> out(5);
        scran_variances::parallel_quantiles(
            5,
            3,
            [&](int b) -> const double* { return stuff[b].data(); },
            std::vector<int>{ 3, 1, 2 },
            1,
            0.5,
            out.data()
        );
        std::vector<double> ref {2, 2, 3, 4, 5};
        scran_tests::compare_almost_equal_containers(ref, out, {});
    }

    // Skip blocks without enough cells.
    {
        std::vector<double> out(5);
        scran_variances::parallel_quantiles(
            5,
            3,
            [&](int b) -> const double* { return stuff[b].data(); },
            std::vector<int>{ 3, 1, 2 },
            2,
            0.5,
            out.data()
        );
        std::vector<double> ref {2, 4.5, 3, 3.5, 3.5};
        scran_tests::compare_almost_equal_containers(ref, out, {});
    }

    // Works with a minimum quantile.
    {
        std::vector<double> out(5);
        scran_variances::parallel_quantiles(
            5,
            3,
            [&](int b) -> const double* { return stuff[b].data(); },
            std::vector<int>{ 3, 1, 2 },
            1,
            0,
            out.data()
        );
        std::vector<double> ref {1, 1, 2, 3, 2};
        scran_tests::compare_almost_equal_containers(ref, out, {});
    }
}
