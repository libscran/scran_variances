#include <gtest/gtest.h>

#include "choose_highly_variable_genes.hpp"
#include "simulate_vector.h"

class ChooseHvgsTest : public ::testing::TestWithParam<std::tuple<int, int> > {};

TEST_P(ChooseHvgsTest, Basic) {
    auto p = GetParam();
    size_t ngenes = std::get<0>(p);
    size_t ntop = std::get<1>(p);

    auto x = simulate_vector(ngenes, /* lower = */ 0, /* upper = */ 10.0, /* seed = */ ngenes * ntop + 42);

    scran::choose_highly_variable_genes::Options opt;
    opt.top = ntop;
    auto output = scran::choose_highly_variable_genes::compute(ngenes, x.data(), opt);

    // Checking that everything is in order.
    EXPECT_LE(std::min(ngenes, ntop), std::accumulate(output.begin(), output.end(), 0)); 

    double min_has = 100, max_lost = 0;
    for (size_t o = 0; o < ngenes; ++o) {
        if (output[o]) {
            min_has = std::min(min_has, x[o]);
        } else {
            max_lost = std::max(max_lost, x[o]);
        }
    }
    EXPECT_TRUE(min_has > max_lost);

    // Checking that it works for smaller values.
    opt.larger = false;
    auto output_low = scran::choose_highly_variable_genes::compute(ngenes, x.data(), opt);
    EXPECT_LE(std::min(ngenes, ntop), std::accumulate(output_low.begin(), output_low.end(), 0)); 

    double max_has = 0, min_lost = 100;
    for (size_t o = 0; o < ngenes; ++o) {
        if (output_low[o]) {
            max_has = std::max(max_has, x[o]);
        } else {
            min_lost = std::min(min_lost, x[o]);
        }
    }
    EXPECT_TRUE(max_has < min_lost);
}

INSTANTIATE_TEST_SUITE_P(
    ChooseHvgs,
    ChooseHvgsTest,
    ::testing::Combine(
        ::testing::Values(11, 111, 1111), // number of values
        ::testing::Values(5, 50, 500) // number of tops
    )
);

TEST(ChooseHvgs, Ties) {
    std::vector<double> x{ 1.1, 2.2, 1.1, 3.3, 0.0, 4.4, 4.4, 3.3 };

    {
        // Ignoring ties.
        scran::choose_highly_variable_genes::Options opt;
        opt.keep_ties = false;
        opt.top = 3;
        auto output = scran::choose_highly_variable_genes::compute(x.size(), x.data(), opt);
        EXPECT_EQ(opt.top, std::accumulate(output.begin(), output.end(), 0));

        // Keeping all ties.
        opt.keep_ties = true;
        auto output2 = scran::choose_highly_variable_genes::compute(x.size(), x.data(), opt);
        EXPECT_LT(opt.top, std::accumulate(output2.begin(), output2.end(), 0));
    }

    {
        // Ignoring ties.
        scran::choose_highly_variable_genes::Options opt;
        opt.larger = false;
        opt.keep_ties = false;
        opt.top = 2;
        auto output = scran::choose_highly_variable_genes::compute(x.size(), x.data(), opt);
        EXPECT_EQ(opt.top, std::accumulate(output.begin(), output.end(), 0));

        // Keeping all ties.
        opt.keep_ties = true;
        auto output2 = scran::choose_highly_variable_genes::compute(x.size(), x.data(), opt);
        EXPECT_LT(opt.top, std::accumulate(output2.begin(), output2.end(), 0));
    }
}

TEST(ChooseHvgs, Extremes) {
    std::vector<double> x{ 1.1, 2.2, 1.1, 3.3, 0.0, 4.4, 4.4, 3.3 };
    scran::choose_highly_variable_genes::Options opt;
    opt.top = 0;
    auto output = scran::choose_highly_variable_genes::compute(x.size(), x.data(), opt);
    EXPECT_EQ(0, std::accumulate(output.begin(), output.end(), 0));
}
