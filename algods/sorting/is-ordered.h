//
// is-ordered.h
//
// Copyright (c) 2020 Carlos Braga
// This program is free software; you can redistribute it and/or modify it
// under the terms of the MIT License. See accompanying LICENSE.md or
// https://opensource.org/licenses/MIT.
//

#ifndef SORTING_IS_ORDERED_H_
#define SORTING_IS_ORDERED_H_

#include <cassert>
#include <vector>
#include <algorithm>
#include <functional>

///
/// @brief Return true if all array elements are ordered according
/// the binary comparator:
///     comp(a[i], a[i+1]) = true, i=lo..hi.
/// Return false otherwise.
///
/// @see Algorithms, 4th Edition by Robert Sedgewick, Kevin Wayne.
///
template<typename Item>
bool IsOrdered(
    std::vector<Item> &arr,
    size_t lo,
    size_t hi,
    std::function<bool(const Item &, const Item &)> const &comp)
{
    assert(arr.size() > 0 && "invalid array size");
    assert(lo >= 0 && lo < arr.size() && "invalid lower limit");
    assert(hi >= 0 && hi < arr.size() && "invalid upper limit");
    assert(lo <= hi && "invalid range");

    for (size_t i = lo; i < hi; ++i) {
        if (!comp(arr[i], arr[i+1])) {
            return false;
        }
    }

    return true;
}

#endif // SORTING_IS_ORDERED_H_
