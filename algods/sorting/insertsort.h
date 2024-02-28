//
// insertsort.h
//
// Copyright (c) 2020 Carlos Braga
// This program is free software; you can redistribute it and/or modify it
// under the terms of the MIT License. See accompanying LICENSE.md or
// https://opensource.org/licenses/MIT.
//

#ifndef SORTING_INSERTSORT_H_
#define SORTING_INSERTSORT_H_

#include <cassert>
#include <vector>
#include <algorithm>
#include <functional>

///
/// @brief InsertSort implements the array insert sort algorithm.
///
/// API InsertSort<Item>:
///      InsertSort()       create an insertion sort object
///      ~InsertSort()      destroy an insertion sort object
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
class InsertSort {
public:
    void operator()(
        std::vector<Item> &arr,
        std::function<bool(const Item &, const Item &)> const &comp) {
        assert(arr.size() > 0 && "invalid size");

        //
        // For each item j <= i, if:
        //      comp(a[j-1], a[j]) == false
        // then the pair j-1 and j violate the invariant, i.e.,
        //      comp(a[j], a[j-1]) == true,
        // Swap items j and j-1 and continue until item j-1 for which
        //      comp(a[j-1], a[j]) == true
        //
        for (size_t i = 1; i < arr.size(); ++i) {
            for (size_t j = i; j > 0 && comp(arr[j], arr[j-1]); --j) {
                 std::swap<Item>(arr[j], arr[j-1]);
            }
        }
    }
};

#endif // SORTING_INSERTSORT_H_
