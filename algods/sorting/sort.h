//
// sort.h
//
// Copyright (c) 2020 Carlos Braga
// This program is free software; you can redistribute it and/or modify it
// under the terms of the MIT License. See accompanying LICENSE.md or
// https://opensource.org/licenses/MIT.
//

#ifndef SORTING_SORT_H_
#define SORTING_SORT_H_

#include <cassert>
#include <vector>
#include <algorithm>
#include <functional>
#include "algods/sorting/insertsort.h"
#include "algods/sorting/mergesort.h"
#include "algods/sorting/quicksort.h"
#include "algods/sorting/selectsort.h"
#include "algods/sorting/shellsort.h"

///
/// @brief Sort is an interface with a function templated with a data type Item,
/// a comparison function Compare<Item> specifying the sorting order, and a class
/// Operate<Item> specifying the sorting algorithm to use.
///
/// Default comparator is std::less and the default sorting algorithm is given by
/// the SelectSort function object.
/// Item is the type of the object being sorted. It is primarly defined as a
/// primitive type. Compare is a binary comparison functor defining the sorting
/// order. It is of the same type as its default type, std::less.
///
/// Operate is a template argument specifying the sorting class. The sorting
/// operator Operate itself is a template taking the array type Item and the
/// binary comparison functor Compare as arguments.
///
/// @example Sort an array of doubles in decreasing order (arr[j] > arr[i]),
///     for j < i, using SelectSort algorithm:
///         Sort<double, std::greater, SelectSort>(arr);
///
template<
    typename Item,
    template<typename> class Compare = std::less,
    template<typename> class Operate = SelectSort>
inline void Sort(
    std::vector<Item> &arr,
    Compare<Item> comp = Compare<Item>(),
    Operate<Item> op = Operate<Item>())
{
    op(arr, comp);
}

///
/// @brief Specialization of the main sorting interface. Take the array of data
/// type Item, the comparison functor Compare<Item> specifying the order, and the
/// sorting class Operate<Item,Compare> specifying the sorting algorithm to use.
/// Instead of sorting the array, return an index vector that sorts the array.
///
/// @example Compute the indices that sort an array of doubles in decreasing order
///     (arr[j] > arr[i]), j < i, using SelectSort algorithm:
///         ArgSort<double, std::greater, SelectSort>(arr);
///
template<
    typename Item,
    template<typename> class Compare = std::less,
    template<typename> class Operate = SelectSort>
inline std::vector<size_t> ArgSort(
    std::vector<Item> &arr,
    Compare<Item> comp = Compare<Item>(),
    Operate<size_t> op = Operate<size_t>())
{
    assert(arr.size() > 0 && "invalid size");

    // Create an index vector the size of the array.
    std::vector<size_t> idx(arr.size());
    for (size_t i = 0; i < arr.size(); ++i) {
        idx[i] = i;
    }

    // Sort the vector using a lambda expression as comparison predicate. The
    // predicate captures the array by reference and compares the items
    // corresponding to the index values.
    std::function<bool(const Item &, const Item &)> idx_comp =
        [&] (const size_t &i, const size_t &j) -> bool {
            return comp(arr[i], arr[j]);
        };
    op(idx, idx_comp);
    return idx;
}

#endif // SORTING_SORT_H_
