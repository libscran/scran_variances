#include <gtest/gtest.h>
#include <random>
#include "simulate_vector.h"
#include "compare_almost_equal.h"
#include "fit_variance_trend.hpp"

TEST(FitVarianceTrendTest, Basic) {
    auto x = simulate_vector(21, /* lower = */ 0.0, /* upper = */ 1.0);
    auto y = x;
    for (auto& y0 : y) { y0 *= 2; } // just to add some variety.

    scran::fit_variance_trend::Options opt;
    opt.transform = false;
    auto output = scran::fit_variance_trend::compute(x.size(), x.data(), y.data(), opt);
    compare_almost_equal(output.fitted, y); // should be an exact fit for a straight line.

    // Again, with the transformation, no filtering.
    y = x;
    for (auto& y0 : y) { y0 *= y0 * y0 * y0; } 

    opt.transform = true;
    opt.mean_filter = false;
    output = scran::fit_variance_trend::compute(x.size(), x.data(), y.data(), opt);
    compare_almost_equal(output.fitted, y);
}

TEST(FitVarianceTrendTest, Extrapolation) {
    std::vector<double> x { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
    std::vector<double> y { 0, 0, 0, 0, 0, 6, 7, 8, 9, 10 };
    
    scran::fit_variance_trend::Options opt;
    opt.transform = false;
    opt.minimum_mean = 5.5;

    auto output = scran::fit_variance_trend::compute(x.size(), x.data(), y.data(), opt);
    compare_almost_equal(output.fitted, x); // should be y = x, extrapolation to the filtered elements.

    // Same result if we shuffle it.
    auto x2 = x;
    std::reverse(x2.begin(), x2.end());
    auto y2 = y;
    std::reverse(y2.begin(), y2.end());

    auto output2 = scran::fit_variance_trend::compute(x2.size(), x2.data(), y2.data(), opt);
    compare_almost_equal(output2.fitted, x2);
}

TEST(FitVarianceTrendTest, Residuals) {
    auto x = simulate_vector(101, /* lower = */ 0.0, /* upper = */ 1.0);
    auto y = simulate_vector(101, /* lower = */ 0.1, /* upper = */ 2.0);

    scran::fit_variance_trend::Options opt;
    auto output = scran::fit_variance_trend::compute(x.size(), x.data(), y.data(), opt);

    // Just repeating the residual calculation here, after transformation.
    std::vector<double> ref(x.size());
    for (size_t i = 0; i < x.size(); ++i) {
        ref[i] = y[i] - output.fitted[i];
    }
    EXPECT_EQ(output.residuals, ref);

    // And again, without transformation.
    opt.transform = false;
    output = scran::fit_variance_trend::compute(x.size(), x.data(), y.data(), opt);

    for (size_t i = 0; i < x.size(); ++i) {
        ref[i] = y[i] - output.fitted[i];
    }
    EXPECT_EQ(output.residuals, ref);
}

TEST(FitVarianceTrendTest, FixedMode) {
    auto x = simulate_vector(101, /* lower = */ 0.0, /* upper = */ 1.0);
    auto y = simulate_vector(101, /* lower = */ 0.1, /* upper = */ 2.0);

    scran::fit_variance_trend::Options opt;
    auto output = scran::fit_variance_trend::compute(x.size(), x.data(), y.data(), opt);

    opt.use_fixed_width = true;
    opt.minimum_window_count = 10;
    opt.fixed_width = 0.2;
    auto foutput = scran::fit_variance_trend::compute(x.size(), x.data(), y.data(), opt);

    EXPECT_NE(output.residuals, foutput.residuals);

    // They eventually converge when both all window widths are at their maximum;
    // either because of a large span, or because we need to get a minimum number of counts. 
    scran::fit_variance_trend::Options opt2;
    opt2.span = 1;
    auto output2 = scran::fit_variance_trend::compute(x.size(), x.data(), y.data(), opt2);

    opt2.use_fixed_width = true;
    opt2.minimum_window_count = 200;
    auto foutput2 = scran::fit_variance_trend::compute(x.size(), x.data(), y.data(), opt2);
    EXPECT_EQ(output2.residuals, foutput2.residuals);

    opt2.minimum_window_count = 0;
    opt2.fixed_width = 10;
    foutput2 = scran::fit_variance_trend::compute(x.size(), x.data(), y.data(), opt2);
    EXPECT_EQ(output2.residuals, foutput2.residuals);
}
