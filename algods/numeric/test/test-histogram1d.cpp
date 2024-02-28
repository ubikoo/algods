//
// test-histogram1d.cpp
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
static const size_t kNumBins = 576;
static const double kWeight = 1.0;

void TestHistogram1d(const std::string &filename, Histogram1d &histogram)
{
    // Clone the histogram
    Histogram1d histogram_new = histogram.clone();

    // Compare histograms.
    {
        REQUIRE(histogram_new.xbins() == histogram.xbins());
        REQUIRE(histogram_new.size() == histogram.size());
        for (size_t i = 0; i < histogram.xbins(); ++i) {
            REQUIRE(histogram_new.count(i) == histogram.count(i));
            REQUIRE(histogram_new.value(i) == histogram.value(i));
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
        REQUIRE(histogram_new.size() == histogram.size());
        for (size_t i = 0; i < histogram.xbins(); ++i) {
            REQUIRE(histogram_new.count(i) == histogram.count(i));
            REQUIRE(histogram_new.value(i) == histogram.value(i));
        }
    }
}

TEST_CASE("Test Histogram1d") {
    SECTION("Uniform distribution") {
        // Compute the histogram.
        std::random_device seed;
        std::mt19937 rng(seed());
        std::uniform_real_distribution<double> dist(0.0, 1.0);
        Histogram1d histogram(kNumBins, 0.0, 1.0);
        for(size_t n = 0; n < kNumSamples; ++n) {
            histogram.bin(dist(rng), kWeight);
        }

        // Test the histogram api.
        TestHistogram1d("/tmp/out.histogram1d_uniform", histogram);
    }

    SECTION("Normal distribution") {
        // Compute the histogram.
        std::random_device seed;
        std::mt19937 rng(seed());
        std::normal_distribution<double> dist(0.0, 1.0);
        Histogram1d histogram(kNumBins, 0.0, 1.0);
        for(size_t n = 0; n < kNumSamples; ++n) {
            histogram.bin(dist(rng), kWeight);
        }

        // Test the histogram api.
        TestHistogram1d("/tmp/out.histogram1d_normal", histogram);
    }
}
