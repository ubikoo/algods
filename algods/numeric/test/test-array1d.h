//
// test-array1d.h
//
// Copyright (c) 2020 Carlos Braga
// This program is free software; you can redistribute it and/or modify it
// under the terms of the MIT License. See accompanying LICENSE.md or
// https://opensource.org/licenses/MIT.
//

#ifndef TEST_NUMERIC_ARRAY1D_H_
#define TEST_NUMERIC_ARRAY1D_H_

#include <string>
#include <fstream>
#include <iostream>
#include "algods/numeric/array.h"

/// ---- Test vector read and write -------------------------------------------
template<typename T>
bool TestVectorReadWrite(size_t n)
{
    std::cout << __PRETTY_FUNCTION__ << " " << typeid(T).name() << "\n";

    // Create a vector of type T
    Vector<T> v1(n);
    for (size_t i = 0; i < n; ++i) {
        v1(i) = static_cast<T>(i);
    }

    // Write/Read vector of type T
    std::string prefix("/tmp/out.vector.");
    std::string filename = prefix + typeid(T).name();

    std::ofstream ofs(filename);
    if (!ofs) {
        std::cerr << "Failed to write file: " << filename << std::endl;
        return false;
    }
    v1.write(ofs);
    ofs.close();

    Vector<T> v2(n);
    std::ifstream ifs(filename);
    v2.read(ifs);
    ifs.close();

    // Check error
    assert(v1.size() == v2.size());
    for (size_t i = 0; i < n; ++i) {
        double err = std::fabs(
            static_cast<double>(v2(i)) - static_cast<double>(v1(i)));
        if (err != 0.0) {
            std::cerr << "err " << err
                << " v1 " << v1(i)
                << " v2 " << v2(i)
                << "\n";
            return false;
        }
    }
    return true;
}

/// ---- Test vector clone ------------------------------------------------------
template<typename T>
bool TestVectorClone(size_t n)
{
    std::cout << __PRETTY_FUNCTION__ << " " << typeid(T).name() << "\n";

    // Create a vector of type T
    Vector<T> v1(n);
    for (size_t i = 0; i < n; ++i) {
        v1(i) = static_cast<T>(i);
    }

    // Create a clone of the vector.
    Vector<T> v2 = std::move(v1.clone());

    // Check error
    assert(v1.size() == v2.size());
    for (size_t i = 0; i < n; ++i) {
        double err = std::fabs(static_cast<double>(v2(i)) -
                               static_cast<double>(v1(i)));
        if (err != 0.0) {
            std::cerr << "err " << err
                << " v1 " << v1(i)
                << " v2 " << v2(i)
                << "\n";
            return false;
        }
    }
    return true;
}

#endif // TEST_NUMERIC_ARRAY1D_H_
