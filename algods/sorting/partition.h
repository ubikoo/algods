//
// partition.h
//
// Copyright (c) 2020 Carlos Braga
// This program is free software; you can redistribute it and/or modify it
// under the terms of the MIT License. See accompanying LICENSE.md or
// https://opensource.org/licenses/MIT.
//

#ifndef SORTING_PARTITION_H_
#define SORTING_PARTITION_H_

#include <cassert>
#include <vector>
#include <algorithm>
#include <functional>

///
/// @brief Partition the subarray a[lo..hi] using the given predicate such
///     p(a[k]) = true,     k <= j
///     p(a[k]) = false,    k > j
/// and return the partitioning index j. For example, if
///     p(a[k]) := (a[k] <= a_p),
/// then the subarray a[lo..hi] would be partitioned such that
///     a[lo..j-1] <= a[j] <= a[j+1..hi],   a[j] = a_p
///
/// Index j is a lower bound of the false-partition or the rightmost element of
/// the true-partition on the left, a[lo..j] = true.
///
/// Index i is a upper bound of the true-partition or the leftmost element on the
/// false-partition on the right, a[i..hi] = false.
///
/// If the returned index j is equal to lo(j==lo) then either
///     - we reached the lo boundary of the subarray withouth finding a
///     single true element - the subarray is all false, or
///     - the subarray contains only one item and could not be
///     partitioned - the result is undefined.
///
/// The meaning of(j==lo) depends on purpose of partitioning. Check the output
/// result depending on the what is being used for.
///
/// @see Algorithms, 4th Edition by Robert Sedgewick, Kevin Wayne.
///
template<typename Item>
size_t Partition(std::vector<Item> &arr,
                  size_t lo,
                  size_t hi,
                  std::function<bool(const Item &)> const &pred)
{
    assert(arr.size() > 0 && "invalid array size");
    assert(lo >= 0 && lo < arr.size() && "invalid lower limit");
    assert(hi >= 0 && hi < arr.size() && "invalid upper limit");
    assert(lo <= hi && "invalid range");

    size_t i = lo - 1;
    size_t j = hi + 1;

    while (true) {
        // Find the item on lo to swap.
        while (pred(arr[++i])) {
            if (i == hi) {
                break;
            }
        }

        // Find the item on hi to swap. If pred(a[j]) compares a[j] with a[lo]
        // (cf. QuickSort), then (j==lo) check is redundant since a[lo] acts as
        // a marker.
        while (!pred(arr[--j])) {
            if (j == lo) {
                break;
            }
        }

        // Check if the pointers cross.
        if (i >= j) {
            return j;
        }

        // The pointers haven't crossed yet, swap the elments and continue.
        std::swap(arr[i], arr[j]);
    }
}

#endif // SORTING_PARTITION_H_
