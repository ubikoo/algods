//
// test-union-find.cpp
//
// Copyright (c) 2020 Carlos Braga
// This program is free software; you can redistribute it and/or modify it
// under the terms of the MIT License. See accompanying LICENSE.md or
// https://opensource.org/licenses/MIT.
//

#include "catch2/catch.hpp"
#include <cstdio>
#include <cmath>
#include <vector>
#include "test-unionfind.h"

static double percolation_threshold = 0.59274621;

/// ---- Test UnionFind -------------------------------------------------------
void perc_statistics(std::vector<double> &data, double &avrg, double &sdev)
{
    avrg = 0.0;
    for (auto it : data) {
        avrg += it;
    }
    avrg /= static_cast<double>(data.size());

    sdev = 0.0;
    for (auto it : data) {
        sdev += (it - avrg) * (it - avrg);
    }
    sdev /= static_cast<double>(data.size());
    sdev = std::sqrt(sdev);
}

TEST_CASE("UnionFind") {
    SECTION("Test UnionFind") {
        size_t ndim = 64;
        size_t trials = 16384;
        std::vector<double> perc_threshold(trials);

        #pragma omp parallel for default(none) \
                shared(ndim, trials, perc_threshold) schedule(static)
        for (size_t i = 0; i < trials; ++i) {
            Percolation perc(ndim);
            perc_threshold[i] = perc.compute();
#if 1
            if (i%1024 == 0) {
                std::printf("trial %4lu, threshold %lf\n",
                            i, perc_threshold[i]);
            }
#endif
        } // omp parallel for

        double perc_mean, perc_sdev;
        perc_statistics(perc_threshold, perc_mean, perc_sdev);
        double perc_conf_lo =
            perc_mean - 1.96*perc_sdev / std::sqrt(static_cast<double>(trials));
        double perc_conf_hi =
            perc_mean + 1.96*perc_sdev / std::sqrt(static_cast<double>(trials));

        std::printf("mean                     = %.15lf\n"
                    "stddev                   = %.15lf\n"
                    "95%% confidence interval = [%.15lf, %.15lf]\n",
                    perc_mean, perc_sdev, perc_conf_lo, perc_conf_hi);
        REQUIRE((percolation_threshold - perc_mean) <= perc_conf_lo);
        REQUIRE((percolation_threshold - perc_mean) <= perc_conf_hi);
    }

    SECTION("Test UnionFindTempl") {
        size_t ndim = 64;
        size_t trials = 16384;
        std::vector<double> perc_threshold(trials);

        #pragma omp parallel for default(none) \
                shared(ndim, trials, perc_threshold) schedule(static)
        for (size_t i = 0; i < trials; ++i) {
            PercolationTempl perc(ndim);
            perc_threshold[i] = perc.compute();
#if 1
            if (i%1024 == 0) {
                std::printf("trial %4lu, threshold %lf\n",
                            i, perc_threshold[i]);
            }
#endif
        } // omp parallel for

        double perc_mean, perc_sdev;
        perc_statistics(perc_threshold, perc_mean, perc_sdev);
        double perc_conf_lo =
            perc_mean - 1.96*perc_sdev / std::sqrt(static_cast<double>(trials));
        double perc_conf_hi =
            perc_mean + 1.96*perc_sdev / std::sqrt(static_cast<double>(trials));

        std::printf("mean                     = %.15lf\n"
                    "stddev                   = %.15lf\n"
                    "95%% confidence interval = [%.15lf, %.15lf]\n",
                    perc_mean, perc_sdev, perc_conf_lo, perc_conf_hi);
        REQUIRE((percolation_threshold - perc_mean) <= perc_conf_lo);
        REQUIRE((percolation_threshold - perc_mean) <= perc_conf_hi);
    }
}
