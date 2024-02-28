//
// selectsort.h
//
// Copyright (c) 2020 Carlos Braga
// This program is free software; you can redistribute it and/or modify it
// under the terms of the MIT License. See accompanying LICENSE.md or
// https://opensource.org/licenses/MIT.
//

#ifndef SORTING_SELECTSORT_H_
#define SORTING_SELECTSORT_H_

#include <cassert>
#include <vector>
#include <algorithm>
#include <functional>

///
/// @brief SelectSort implements the array selection sort algorithm.
///
/// API SelectSort<Item>:
///      SelectSort()       create an insertion sort object
///      ~SelectSort()      destroy an insertion sort object
///
///      void sort() (std::vector<Item> &arr, Compare<Item> comp)
///          sort an array, keeping the invariant comp(a[j],a[i]) for j < i.
///          The comparison predicate defines the order of items in the array.
///          The default is std::less<Item> and sorts the array in increasing
///          order, (a[j] < a[i]) for j < i.
///
///
/// @note The binary predicate is declared as a std::function to allow the use of
/// comparison callbacks of different types - lambdas, functors and function
/// pointers.
///
/// @see Algorithms, 4th Edition by Robert Sedgewick, Kevin Wayne.
/// https://stackoverflow.com/questions/14677997
/// https://stackoverflow.com/questions/36030589
///
template<typename Item>
class SelectSort {
public:
    void operator()(
        std::vector<Item> &arr,
        std::function<bool(const Item &, const Item &)> const &comp) {
        assert(arr.size() > 0 && "invalid size");

        // For each item i, select the lowest order item k in subarray a[i..hi]
        // and swap with i to keep the invariant.
        for (size_t i = 0; i < arr.size(); ++i) {
            size_t k = i;
            for (size_t j = i+1; j < arr.size(); ++j) {
                if (comp(arr[j], arr[k])) {
                    k = j;
                }
            }
            std::swap<Item>(arr[i], arr[k]);
        }
    }
};

#endif // SORTING_SELECTSORT_H_
