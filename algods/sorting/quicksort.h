//
// quicksort.h
//
// Copyright (c) 2020 Carlos Braga
// This program is free software; you can redistribute it and/or modify it
// under the terms of the MIT License. See accompanying LICENSE.md or
// https://opensource.org/licenses/MIT.
//

#ifndef SORTING_QUICKSORT_H_
#define SORTING_QUICKSORT_H_

#include <cassert>
#include <vector>
#include <stack>
#include <algorithm>
#include <functional>
#include "algods/sorting/partition.h"

///
/// @brief QuickSort implements the array quicksort algorithm.
///
/// API QuickSort<Item>:
///      QuickSort()       create an insertion sort object
///      ~QuickSort()      destroy an insertion sort object
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
class QuickSort {
public:
    void operator()(
        std::vector<Item> &arr,
        std::function<bool(const Item &, const Item &)> const &comp) {
        assert(arr.size() > 0 && "invalid size");

        // Random shuffle the array to avoid the worst case.
        shuffle(arr);

        // Use a queue of partitions to sort.
        // Initial partition is the entire array.
        using part_pair_t = std::pair<size_t, size_t>;
        using part_stack_t = std::stack<part_pair_t>;

        part_stack_t part_stack;
        part_stack.push(part_pair_t(0, arr.size()-1));
        while (!part_stack.empty()) {
            // Pop the partition ranges.
            part_pair_t p = part_stack.top();
            part_stack.pop();

            // Partition the array and swap a[lo] with partitioning element a[k].
            size_t lo = p.first;
            size_t hi = p.second;
            if (lo >= hi) {
                continue;
            }

            std::function<bool(const Item &)> pred =
                [&] (const Item &item) -> bool { return comp(item, arr[lo]); };
            size_t k = Partition(arr, lo, hi, pred);
            std::swap(arr[lo], arr[k]);

            //
            // Push the two subpartitions into the queue. We are doing arithmetic
            // on unsigned integers, so we need to take care that k-1 does not
            // overflow when k = 0.
            //
            // lo value is guaranteed to be lo >= 0 since 0 is its initial value.
            // hi value is guaranteed to be hi <= n-1 since n-1 is its initial
            // value.
            //
            // If the partitioning index k is equal to lo, then either
            //      - the subarray a[lo..hi] has only one element
            //      - the order of (a[j], a[lo]) is false for j >= lo.
            // Either way, there is no left subpartition [lo..k-1].
            //
            // Conversely, if the partitioning index k is equal to hi, then there
            // is no false-subpartition to add because either
            //      - the subarray a[lo..hi] has only one element or
            //      - the order of (a[j], a[lo]) is true for j <= hi.
            // Either way, there is not right subpartition [k+1..hi]
            //
            if (k > lo) {
                part_stack.push(part_pair_t(lo, k-1));
            }
            if (k < hi) {
                part_stack.push(part_pair_t(k+1, hi));
            }
        }
    }

private:
    // Random shuffle an array using Knuth's shuffle algorithm. For each element
    // 0 <= i < n-1, sample a random index uniformly distributed on the closed
    // interval [i, arr.size()-1] and swap.
    void shuffle(std::vector<Item> &arr) {
        std::random_device rd;      // rng device
        std::mt19937 rng(rd());     // rng engine seeded with rd
        for (size_t i = 0; i < arr.size(); ++i) {
            std::uniform_int_distribution<size_t> dist(i, arr.size()-1);
            std::swap(arr[i], arr[dist(rng)]);
        }
    }
};

#endif // SORTING_QUICKSORT_H_
