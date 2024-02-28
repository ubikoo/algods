//
// test-histogram2d.cpp
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
static const size_t kNumBins = 24;
static const double kWeight = 1.0;

void TestHistogram2d(const std::string &filename, Histogram2d &histogram)
{
    // Clone the histogram
    Histogram2d histogram_new = histogram.clone();

    // Compare histograms.
    {
        REQUIRE(histogram_new.xbins() == histogram.xbins());
        REQUIRE(histogram_new.ybins() == histogram.ybins());
        REQUIRE(histogram_new.size() == histogram.size());
        for (size_t i = 0; i < histogram.xbins(); ++i) {
            for (size_t j = 0; j < histogram.ybins(); ++j) {
                REQUIRE(histogram_new.count(i, j) == histogram.count(i, j));
                REQUIRE(histogram_new.value(i, j) == histogram.value(i, j));
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
        REQUIRE(histogram_new.size() == histogram.size());
        for (size_t i = 0; i < histogram.xbins(); ++i) {
            for (size_t j = 0; j < histogram.ybins(); ++j) {
                REQUIRE(histogram_new.count(i, j) == histogram.count(i, j));
                REQUIRE(histogram_new.value(i, j) == histogram.value(i, j));
            }
        }
    }
}

TEST_CASE("Test Histogram2d") {
    SECTION("Uniform distribution") {
        // Compute the histogram.
        std::random_device seed;
        std::mt19937 rng(seed());
        std::uniform_real_distribution<double> dist(0.0, 1.0);
        Histogram2d histogram(kNumBins, kNumBins, 0.0, 1.0, 0.0, 1.0);
        for(size_t n = 0; n < kNumSamples; ++n) {
            histogram.bin(dist(rng), dist(rng), kWeight);
        }

        // Test the histogram api.
        TestHistogram2d("/tmp/out.histogram2d_uniform", histogram);
    }

    SECTION("Normal distribution") {
        // Compute the histogram.
        std::random_device seed;
        std::mt19937 rng(seed());
        std::normal_distribution<double> dist(0.0, 1.0);
        Histogram2d histogram(kNumBins, kNumBins, 0.0, 1.0, 0.0, 1.0);
        for(size_t n = 0; n < kNumSamples; ++n) {
            histogram.bin(dist(rng), dist(rng), kWeight);
        }

        // Test the histogram api.
        TestHistogram2d("/tmp/out.histogram2d_normal", histogram);
    }
}
