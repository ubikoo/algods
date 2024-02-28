//
// array.h
//
// Copyright (c) 2020 Carlos Braga
// This program is free software; you can redistribute it and/or modify it
// under the terms of the MIT License. See accompanying LICENSE.md or
// https://opensource.org/licenses/MIT.
//

#ifndef NUMERIC_ARRAY_H_
#define NUMERIC_ARRAY_H_

#include "algods/numeric/array1d.h"
#include "algods/numeric/array2d.h"
#include "algods/numeric/array3d.h"

template<typename T>
using Vector = Array1d<T>;
template<typename T>
using Matrix = Array2d<T>;
template<typename T>
using Tensor = Array3d<T>;

#endif // NUMERIC_ARRAY_H_
