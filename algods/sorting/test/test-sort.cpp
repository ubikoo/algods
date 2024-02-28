//
// test-sort.cpp
//
// Copyright (c) 2020 Carlos Braga
// This program is free software; you can redistribute it and/or modify it
// under the terms of the MIT License. See accompanying LICENSE.md or
// https://opensource.org/licenses/MIT.
//

#include "catch2/catch.hpp"
#include "test-sort.h"
#include "test-insertsort.h"
#include "test-mergesort.h"
#include "test-quicksort.h"
#include "test-selectsort.h"
#include "test-shellsort.h"

/// ---- Test Sort Algorithm ---------------------------------------------------
///
TEST_CASE("Sort") {
    size_t n_small = 8192;
    size_t n_large = 262144;

    // ----------------------------------------------------------------------
    SECTION("Test Select Sort") {
        test_selectsort<long double>(n_small);
        test_selectsort<double>(n_small);
        test_selectsort<float>(n_small);

        test_selectsort<long long>(n_small);
        test_selectsort<long>(n_small);
        test_selectsort<int>(n_small);
        test_selectsort<short>(n_small);
        test_selectsort<char>(n_small);

        test_selectsort<unsigned long long>(n_small);
        test_selectsort<unsigned long>(n_small);
        test_selectsort<unsigned int>(n_small);
        test_selectsort<unsigned short>(n_small);
        test_selectsort<unsigned char>(n_small);

        test_selectsort<int64_t>(n_small);
        test_selectsort<int32_t>(n_small);
        test_selectsort<int16_t>(n_small);
        test_selectsort<int8_t>(n_small);

        test_selectsort<uint64_t>(n_small);
        test_selectsort<uint32_t>(n_small);
        test_selectsort<uint16_t>(n_small);
        test_selectsort<uint8_t>(n_small);
    }

    // ----------------------------------------------------------------------
    SECTION("Test Insert Sort") {
        test_insertsort<long double>(n_small);
        test_insertsort<double>(n_small);
        test_insertsort<float>(n_small);

        test_insertsort<long long>(n_small);
        test_insertsort<long>(n_small);
        test_insertsort<int>(n_small);
        test_insertsort<short>(n_small);
        test_insertsort<char>(n_small);

        test_insertsort<unsigned long long>(n_small);
        test_insertsort<unsigned long>(n_small);
        test_insertsort<unsigned int>(n_small);
        test_insertsort<unsigned short>(n_small);
        test_insertsort<unsigned char>(n_small);

        test_insertsort<int64_t>(n_small);
        test_insertsort<int32_t>(n_small);
        test_insertsort<int16_t>(n_small);
        test_insertsort<int8_t>(n_small);

        test_insertsort<uint64_t>(n_small);
        test_insertsort<uint32_t>(n_small);
        test_insertsort<uint16_t>(n_small);
        test_insertsort<uint8_t>(n_small);
    }

    // ----------------------------------------------------------------------
    SECTION("Test Shell Sort") {
        test_shellsort<long double>(n_large);
        test_shellsort<double>(n_large);
        test_shellsort<float>(n_large);

        test_shellsort<long long>(n_large);
        test_shellsort<long>(n_large);
        test_shellsort<int>(n_large);
        test_shellsort<short>(n_large);
        test_shellsort<char>(n_large);

        test_shellsort<unsigned long long>(n_large);
        test_shellsort<unsigned long>(n_large);
        test_shellsort<unsigned int>(n_large);
        test_shellsort<unsigned short>(n_large);
        test_shellsort<unsigned char>(n_large);

        test_shellsort<int64_t>(n_large);
        test_shellsort<int32_t>(n_large);
        test_shellsort<int16_t>(n_large);
        test_shellsort<int8_t>(n_large);

        test_shellsort<uint64_t>(n_large);
        test_shellsort<uint32_t>(n_large);
        test_shellsort<uint16_t>(n_large);
        test_shellsort<uint8_t>(n_large);
    }

    // ----------------------------------------------------------------------
    SECTION("Test Merge Sort") {
        test_mergesort<long double>(n_large);
        test_mergesort<double>(n_large);
        test_mergesort<float>(n_large);

        test_mergesort<long long>(n_large);
        test_mergesort<long>(n_large);
        test_mergesort<int>(n_large);
        test_mergesort<short>(n_large);
        test_mergesort<char>(n_large);

        test_mergesort<unsigned long long>(n_large);
        test_mergesort<unsigned long>(n_large);
        test_mergesort<unsigned int>(n_large);
        test_mergesort<unsigned short>(n_large);
        test_mergesort<unsigned char>(n_large);

        test_mergesort<int64_t>(n_large);
        test_mergesort<int32_t>(n_large);
        test_mergesort<int16_t>(n_large);
        test_mergesort<int8_t>(n_large);

        test_mergesort<uint64_t>(n_large);
        test_mergesort<uint32_t>(n_large);
        test_mergesort<uint16_t>(n_large);
        test_mergesort<uint8_t>(n_large);
    }

    // ----------------------------------------------------------------------
    SECTION("Test Quick Sort") {
        test_quicksort<long double>(n_large);
        test_quicksort<double>(n_large);
        test_quicksort<float>(n_large);

        test_quicksort<long long>(n_large);
        test_quicksort<long>(n_large);
        test_quicksort<int>(n_large);
        test_quicksort<short>(n_large);
        test_quicksort<char>(n_large);

        test_quicksort<unsigned long long>(n_large);
        test_quicksort<unsigned long>(n_large);
        test_quicksort<unsigned int>(n_large);
        test_quicksort<unsigned short>(n_large);
        test_quicksort<unsigned char>(n_large);

        test_quicksort<int64_t>(n_large);
        test_quicksort<int32_t>(n_large);
        test_quicksort<int16_t>(n_large);
        test_quicksort<int8_t>(n_large);

        test_quicksort<uint64_t>(n_large);
        test_quicksort<uint32_t>(n_large);
        test_quicksort<uint16_t>(n_large);
        test_quicksort<uint8_t>(n_large);
    }
}
