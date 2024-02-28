//
// test-histogram3d.cpp
//
// Copyright (c) 2020 Carlos Braga
// This program is free software; you can redistribute it and/or modify it
// under the terms of the MIT License. See accompanying LICENSE.md or
// https://opensource.org/licenses/MIT.
//

#include "catch2/catch.hpp"
#include <random>
#include <string>
#include <fstream>
#include <iostream>
#include "algods/numeric/histogram.h"

static const size_t kNumSamples = 4194304;
static const size_t kNumBins = 8;
static const double kWeight = 1.0;

void TestHistogram3d(const std::string &filename, Histogram3d &histogram)
{
    // Clone the histogram
    Histogram3d histogram_new = histogram.clone();

    // Compare histograms.
    {
        REQUIRE(histogram_new.xbins() == histogram.xbins());
        REQUIRE(histogram_new.ybins() == histogram.ybins());
        REQUIRE(histogram_new.zbins() == histogram.zbins());
        REQUIRE(histogram_new.size() == histogram.size());
        for (size_t i = 0; i < histogram.xbins(); ++i) {
            for (size_t j = 0; j < histogram.ybins(); ++j) {
                for (size_t k = 0; k < histogram.zbins(); ++k) {
                    REQUIRE(histogram_new.count(i, j, k) ==
                            histogram.count(i, j, k));
                    REQUIRE(histogram_new.value(i, j, k) ==
                            histogram.value(i, j, k));
                }
            }
        }
    }

    // Read/Write histogram.
    std::ofstream ofs(filename);
    histogram.write(ofs);
    ofs.close();

    std::ifstream ifs(filename);
    histogram_new.read(ifs);
    ifs.close();

    // Compare histograms.
    {
        REQUIRE(histogram_new.xbins() == histogram.xbins());
        REQUIRE(histogram_new.ybins() == histogram.ybins());
        REQUIRE(histogram_new.zbins() == histogram.zbins());
        REQUIRE(histogram_new.size() == histogram.size());
        for (size_t i = 0; i < histogram.xbins(); ++i) {
            for (size_t j = 0; j < histogram.ybins(); ++j) {
                for (size_t k = 0; k < histogram.zbins(); ++k) {
                    REQUIRE(histogram_new.count(i, j, k) ==
                            histogram.count(i, j, k));
                    REQUIRE(histogram_new.value(i, j, k) ==
                            histogram.value(i, j, k));
                }
            }
        }
    }

}

TEST_CASE("Test Histogram3d") {
    SECTION("Uniform distribution") {
        // Compute the histogram.
        std::random_device seed;
        std::mt19937 rng(seed());
        std::uniform_real_distribution<double> dist(0.0, 1.0);
        Histogram3d histogram(kNumBins, kNumBins, kNumBins,
            0.0, 1.0, 0.0, 1.0, 0.0, 1.0);
        for(size_t n = 0; n < kNumSamples; ++n) {
            histogram.bin(dist(rng), dist(rng), dist(rng), kWeight);
        }

        // Test the histogram api.
        TestHistogram3d("/tmp/out.histogram3d_uniform", histogram);
    }

    SECTION("Normal distribution") {
        // Compute the histogram.
        std::random_device seed;
        std::mt19937 rng(seed());
        std::normal_distribution<double> dist(0.0, 1.0);
        Histogram3d histogram(kNumBins, kNumBins, kNumBins,
            0.0, 1.0, 0.0, 1.0, 0.0, 1.0);
        for(size_t n = 0; n < kNumSamples; ++n) {
            histogram.bin(dist(rng), dist(rng), dist(rng), kWeight);
        }

        // Test the histogram api.
        TestHistogram3d("/tmp/out.histogram3d_normal", histogram);
    }
}
