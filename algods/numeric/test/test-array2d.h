//
// test-array2d.h
//
// Copyright (c) 2020 Carlos Braga
// This program is free software; you can redistribute it and/or modify it
// under the terms of the MIT License. See accompanying LICENSE.md or
// https://opensource.org/licenses/MIT.
//

#ifndef TEST_NUMERIC_ARRAY2D_H_
#define TEST_NUMERIC_ARRAY2D_H_

#include <string>
#include <fstream>
#include <iostream>
#include "algods/numeric/array.h"

/// ---- Test matrix read and write -------------------------------------------
template<typename T>
bool TestMatrixReadWrite(size_t n)
{
    std::cout << __PRETTY_FUNCTION__ << " " << typeid(T).name() << "\n";

    // Create a matrix of type T
    Matrix<T> v1(n,n);
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            v1(i,j) = static_cast<T>(i+j);
        }
    }

    // Write/Read matrix of type T
    std::string prefix("/tmp/out.matrix.");
    std::string filename = prefix + typeid(T).name();

    std::ofstream ofs(filename);
    if (!ofs) {
        std::cerr << "Failed to write file: " << filename << std::endl;
        return false;
    }
    v1.write(ofs);
    ofs.close();

    Matrix<T> v2(n,n);
    std::ifstream ifs(filename);
    v2.read(ifs);
    ifs.close();

    // Check error
    assert(v1.size() == v2.size());
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            double err = std::fabs(static_cast<double>(v2(i,j)) -
                                   static_cast<double>(v1(i,j)));
            if (err != 0.0) {
                std::cerr << "err " << err
                    << " v1 " << v1(i,j)
                    << " v2 " << v2(i,j)
                    << "\n";
                return false;
            }
        }
    }
    return true;
}

/// ---- Test matrix clone ------------------------------------------------------
template<typename T>
bool TestMatrixClone(size_t n)
{
    std::cout << __PRETTY_FUNCTION__ << " " << typeid(T).name() << "\n";

    // Create a matrix of type T
    Matrix<T> v1(n,n);
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            v1(i,j) = static_cast<T>(i+j);
        }
    }

    // Create a clone of the matrix.
    Matrix<T> v2 = std::move(v1.clone());

    // Check error
    assert(v1.size() == v2.size());
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            double err = std::fabs(static_cast<double>(v2(i,j)) -
                                   static_cast<double>(v1(i,j)));
            if (err != 0.0) {
                std::cerr << "err " << err
                    << " v1 " << v1(i,j)
                    << " v2 " << v2(i,j)
                    << "\n";
                return false;
            }
        }
    }
    return true;
}

#endif // TEST_NUMERIC_ARRAY2D_H_
