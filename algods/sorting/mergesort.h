//
// mergesort.h
//
// Copyright (c) 2020 Carlos Braga
// This program is free software; you can redistribute it and/or modify it
// under the terms of the MIT License. See accompanying LICENSE.md or
// https://opensource.org/licenses/MIT.
//

#ifndef SORTING_MERGESORT_H_
#define SORTING_MERGESORT_H_

#include <cassert>
#include <vector>
#include <algorithm>
#include <functional>
#include "algods/sorting/merge.h"
#include "algods/sorting/is-ordered.h"

///
/// @brief MergeSort implements the array merge sort algorithm. The implementation
/// uses a bottom-up merge algorithm. Start by merging pairs of individual
/// items(lo, mid, hi)
///     [0),1] [1),2] ...
/// then merge the sorted pairs
///     [0,1),2,3] [4,5),6,7] ...
/// until the final merge of the two halves
///     [0,1,2,3,4,5,6,...,mid),mid+1,...,n]
/// where mid = 2^(k-1), k=ln(n).
///
/// API MergeSort<Item>:
///      MergeSort()       create an insertion sort object
///      ~MergeSort()      destroy an insertion sort object
///
///      void sort() (Vector<Item> &arr, Compare<Item> comp)
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
class MergeSort {
public:
    void operator()(
        std::vector<Item> &arr,
        std::function<bool(const Item &, const Item &)> const &comp) {
        assert(arr.size() > 0 && "invalid size");

        // Do lg(n) passes of pairwise merges.
        const size_t n = arr.size();
        std::vector<Item> aux(arr);
        for (size_t len = 1; len < n; len = len+len) {
            for (size_t lo = 0; lo < n-len; lo += len+len) {
                size_t mid = lo + len-1;
                size_t hi = std::min(lo + len+len-1, n-1);
                Merge<Item>(arr, aux, lo, mid, hi, comp);
            }
        }
    }
};

#endif // SORTING_MERGESORT_H_
