//
// shellsort.h
//
// Copyright (c) 2020 Carlos Braga
// This program is free software; you can redistribute it and/or modify it
// under the terms of the MIT License. See accompanying LICENSE.md or
// https://opensource.org/licenses/MIT.
//

#ifndef SORTING_SHELLSORT_H_
#define SORTING_SHELLSORT_H_

#include <cassert>
#include <vector>
#include <algorithm>
#include <functional>

///
/// @brief ShellSort implements the array shell sort algorithm.
///
/// ShellSort is a simple extension of insertion sort where instead of insert
/// sort with steps of 1, we insert sort with steps of h.
///
/// For example, when sorting a N size array in increasing order, if the smallest
/// element is the last, one would need N-1 exchanges to get the element to its
/// initial position. By taking larger steps one is able to bring the element to
/// its destination faster.
///
/// By taking a series of h-sorted sequences with smaller h-steps, less number of
/// exchanges would need to be performed to achieve the final order. The last
/// h-step would be 1 at which point the sort is equivalent to insertion sort.
/// However, at that point, the array is partially sorted already so less number
/// of exchanges are needed.
///
/// API ShellSort<Item>:
///      ShellSort()       create an insertion sort object
///      ~ShellSort()      destroy an insertion sort object
///
///      void sort() (std::vector<Item> &arr, Compare<Item> comp)
///          sort an array, keeping the invariant comp(a[j],a[i]) for j < i.
///          The comparison predicate defines the order of items in the array.
///          The default is std::less<Item> and sorts the array in increasing
///          order, (a[j] < a[i]) for j < i.
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
class ShellSort {
public:
    void operator()(
        std::vector<Item> &arr,
        std::function<bool(const Item &, const Item &)> const &comp) {
        assert(arr.size() > 0 && "invalid size");

        //
        // Use an increment sequence of decreasing values(3^k - 1)/2, starting
        // at the smallest increment greater than or equal to floor(n/3) and
        // decreasing to 1.
        //
        // Compute the largest h-step that can be used for an array of size n.
        // For each h-step, h-sort the array and repeat until the final sort with
        // h-step = 1.
        //
        // Example:
        // h = 4
        // i = 4,5,6,7,8,...,n
        // j = (4,0) (5,1) (6,2) (7,3)
        //     (8,4,0) (9,5,1) (12,6,2) (11,7,3)
        //     (12,8,4,0) (13,9,5,1) (14,12,6,2) (15,11,7,3)
        //     (16,12,8,4,0) (17,13,9,5,1) (18,14,12,6,2) (19,15,11,7,3)
        //     ...
        //     (n,n-4,n-8,...,min(0, h - n%h))
        //
        const size_t n = arr.size();
        size_t h = 1;
        while (h < n/3) {
            h = 3*h + 1;    // h = 1, 4, 13, 40,..., >= (n/3)
        }

        while (h >= 1) {
            // h-sort the array
            for (size_t i = h; i < n; ++i) {
                // insert a[i] among a[i-h], a[i-2h], a[i-3h]
                for (size_t j = i; j >= h && comp(arr[j], arr[j-h]); j -= h) {
                    std::swap<Item>(arr[j], arr[j-h]);
                }
            }
            h /= 3;
        }
    }
};

#endif // SORTING_SHELLSORT_H_
