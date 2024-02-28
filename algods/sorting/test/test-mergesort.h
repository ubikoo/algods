//
// test-mergesort.h
//
// Copyright (c) 2020 Carlos Braga
// This program is free software; you can redistribute it and/or modify it
// under the terms of the MIT License. See accompanying LICENSE.md or
// https://opensource.org/licenses/MIT.
//

#ifndef TEST_MERGESORT_H_
#define TEST_MERGESORT_H_

#include "algods/sorting/sort.h"

/// ---- Test sort read and write ----------------------------------------------
///
template<typename T>
void test_mergesort(size_t n)
{
    std::printf("%s %s\n", __PRETTY_FUNCTION__, typeid(T).name());

    // Create a vector of type T and the index vector
    std::vector<T> v = generate_vector<T>(n);

    // Test Sort
    Sort<T, std::less_equal, MergeSort>(v);
    for (size_t i = 0; i < v.size()-1; ++i) {
        REQUIRE(v[i] <= v[i+1]);                 // sort in decreasing order
    }

    Sort<T, std::greater_equal, MergeSort>(v);
    for (size_t i = 0; i < v.size()-1; ++i) {
        REQUIRE(v[i] >= v[i+1]);                 // sort in increasing order
    }

    // Test ArgSort
    std::vector<size_t> ltix =
        ArgSort<T, std::less_equal, MergeSort>(v);
    for (size_t i = 0; i < v.size()-1; ++i) {
        REQUIRE(v[ltix[i]] <= v[ltix[i+1]]);     // sort in decreasing order
    }

    std::vector<size_t> gtix =
        ArgSort<T, std::greater_equal, MergeSort>(v);
    for (size_t i = 0; i < v.size()-1; ++i) {
        REQUIRE(v[gtix[i]] >= v[gtix[i+1]]);     // sort in increasing order
    }
}

#endif // TEST_MERGESORT_H_
