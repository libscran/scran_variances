#include <gtest/gtest.h>
#include <random>

#include "scran_tests/scran_tests.hpp"

#include "scran_variances/fit_variance_trend.hpp"

TEST(FitVarianceTrendTest, Basic) {
    auto x = scran_tests::simulate_vector(21, []{
        scran_tests::SimulationParameters sparams;
        sparams.lower = 0;
        sparams.upper = 1;
        return sparams;
    }());
    auto y = x;
    for (auto& y0 : y) { y0 *= 2; } // just to add some variety.

    scran_variances::FitVarianceTrendOptions opt;
    opt.transform = false;
    auto output = scran_variances::fit_variance_trend(x.size(), x.data(), y.data(), opt);
    scran_tests::compare_almost_equal(output.fitted, y); // should be an exact fit for a straight line.

    // Again, with the transformation, no filtering.
    y = x;
    for (auto& y0 : y) { y0 *= y0 * y0 * y0; } 

    opt.transform = true;
    opt.mean_filter = false;
    output = scran_variances::fit_variance_trend(x.size(), x.data(), y.data(), opt);
    scran_tests::compare_almost_equal(output.fitted, y);
}

TEST(FitVarianceTrendTest, Extrapolation) {
    std::vector<double> x { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
    std::vector<double> y { 0, 0, 0, 0, 0, 6, 7, 8, 9, 10 };
    
    scran_variances::FitVarianceTrendOptions opt;
    opt.transform = false;
    opt.minimum_mean = 5.5;

    auto output = scran_variances::fit_variance_trend(x.size(), x.data(), y.data(), opt);
    scran_tests::compare_almost_equal(output.fitted, x); // should be y = x, extrapolation to the filtered elements.

    // Same result if we shuffle it.
    auto x2 = x;
    std::reverse(x2.begin(), x2.end());
    auto y2 = y;
    std::reverse(y2.begin(), y2.end());

    auto output2 = scran_variances::fit_variance_trend(x2.size(), x2.data(), y2.data(), opt);
    scran_tests::compare_almost_equal(output2.fitted, x2);
}

TEST(FitVarianceTrendTest, Residuals) {
    scran_tests::SimulationParameters sparams;
    sparams.lower = 0;
    sparams.upper = 1;
    sparams.seed = 42;
    auto x = scran_tests::simulate_vector(101, sparams);
    sparams.seed = 69;
    auto y = scran_tests::simulate_vector(101, sparams);

    scran_variances::FitVarianceTrendOptions opt;
    auto output = scran_variances::fit_variance_trend(x.size(), x.data(), y.data(), opt);

    // Just repeating the residual calculation here, after transformation.
    std::vector<double> ref(x.size());
    for (size_t i = 0; i < x.size(); ++i) {
        ref[i] = y[i] - output.fitted[i];
    }
    EXPECT_EQ(output.residuals, ref);

    // And again, without transformation.
    opt.transform = false;
    output = scran_variances::fit_variance_trend(x.size(), x.data(), y.data(), opt);

    for (size_t i = 0; i < x.size(); ++i) {
        ref[i] = y[i] - output.fitted[i];
    }
    EXPECT_EQ(output.residuals, ref);
}

TEST(FitVarianceTrendTest, Filtering) {
    auto x = scran_tests::simulate_vector(1001, []{
        scran_tests::SimulationParameters sparams;
        sparams.lower = 0;
        sparams.upper = 1;
        sparams.seed = 420;
        return sparams;
    }());
    auto y = scran_tests::simulate_vector(1001, []{ 
        scran_tests::SimulationParameters sparams;
        sparams.lower = 0.1;
        sparams.upper = 2;
        sparams.seed = 80085;
        return sparams;
    }());

    scran_variances::FitVarianceTrendOptions opt;
    auto ref = scran_variances::fit_variance_trend(x.size(), x.data(), y.data(), opt);

    opt.mean_filter = false;
    auto output_unfilt = scran_variances::fit_variance_trend(x.size(), x.data(), y.data(), opt);
    EXPECT_NE(ref.residuals, output_unfilt.residuals); // check that there is a difference.

    std::vector<double> submean, subvar, subfit, subresid;
    for (size_t i = 0; i < x.size(); ++i) {
        if (x[i] >= opt.minimum_mean) {
            submean.push_back(x[i]);
            subvar.push_back(y[i]);
            subfit.push_back(ref.fitted[i]);
            subresid.push_back(ref.residuals[i]);
        }
    }

    auto output_manual = scran_variances::fit_variance_trend(submean.size(), submean.data(), subvar.data(), opt);
    EXPECT_EQ(output_manual.residuals, subresid);
    EXPECT_EQ(output_manual.fitted, subfit);
}

TEST(FitVarianceTrendTest, MinWidth) {
    auto x = scran_tests::simulate_vector(101, []{
        scran_tests::SimulationParameters sparams;
        sparams.lower = 0;
        sparams.upper = 1;
        sparams.seed = 12345;
        return sparams;
    }());
    auto y = scran_tests::simulate_vector(101, []{ 
        scran_tests::SimulationParameters sparams;
        sparams.lower = 0.1;
        sparams.upper = 2;
        sparams.seed = 6789;
        return sparams;
    }());

    scran_variances::FitVarianceTrendOptions opt;
    auto output = scran_variances::fit_variance_trend(x.size(), x.data(), y.data(), opt);

    opt.use_minimum_width = true;
    opt.minimum_window_count = 10;
    opt.minimum_width = 0.2;
    auto foutput = scran_variances::fit_variance_trend(x.size(), x.data(), y.data(), opt);

    EXPECT_NE(output.residuals, foutput.residuals);

    // They eventually converge when both all window widths are at their maximum;
    // either because of a large span, or because we need to get a minimum number of counts. 
    scran_variances::FitVarianceTrendOptions opt2;
    opt2.span = 1;
    auto output2 = scran_variances::fit_variance_trend(x.size(), x.data(), y.data(), opt2);

    opt2.use_minimum_width = true;
    opt2.minimum_window_count = 200;
    auto foutput2 = scran_variances::fit_variance_trend(x.size(), x.data(), y.data(), opt2);
    EXPECT_EQ(output2.residuals, foutput2.residuals);

    opt2.minimum_window_count = 0;
    opt2.minimum_width = 10;
    foutput2 = scran_variances::fit_variance_trend(x.size(), x.data(), y.data(), opt2);
    EXPECT_EQ(output2.residuals, foutput2.residuals);
}
