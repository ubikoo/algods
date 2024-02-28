//
// test-array3d.h
//
// Copyright (c) 2020 Carlos Braga
// This program is free software; you can redistribute it and/or modify it
// under the terms of the MIT License. See accompanying LICENSE.md or
// https://opensource.org/licenses/MIT.
//

#ifndef TEST_NUMERIC_ARRAY3D_H_
#define TEST_NUMERIC_ARRAY3D_H_

#include <string>
#include <fstream>
#include <iostream>
#include "algods/numeric/array.h"

/// ---- Test tensor read and write -------------------------------------------
template<typename T>
bool TestTensorReadWrite(size_t n)
{
    std::cout << __PRETTY_FUNCTION__ << " " << typeid(T).name() << "\n";

    // Create a tensor of type T
    Tensor<T> v1(n,n,n);
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            for (size_t k = 0; k < n; ++k) {
                v1(i,j,k) = static_cast<T>(i+j+k);
            }
        }
    }

    // Write/Read tensor of type T
    std::string prefix("/tmp/out.tensor.");
    std::string filename = prefix + typeid(T).name();

    std::ofstream ofs(filename);
    if (!ofs) {
        std::cerr << "Failed to write file: " << filename << std::endl;
        return false;
    }
    v1.write(ofs);
    ofs.close();

    Tensor<T> v2(n,n,n);
    std::ifstream ifs(filename);
    v2.read(ifs);
    ifs.close();

    // Check error
    assert(v1.size() == v2.size());
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            for (size_t k = 0; k < n; ++k) {
                double err = std::fabs(static_cast<double>(v2(i,j,k)) -
                                       static_cast<double>(v1(i,j,k)));
                if (err != 0.0) {
                    std::cerr << "err " << err
                        << " v1 " << v1(i,j,k)
                        << " v2 " << v2(i,j,k)
                        << "\n";
                    return false;
                }
            }
        }
    }
    return true;
}

/// ---- Test tensor clone ------------------------------------------------------
template<typename T>
bool TestTensorClone(size_t n)
{
    std::cout << __PRETTY_FUNCTION__ << " " << typeid(T).name() << "\n";

    // Create a tensor of type T
    Tensor<T> v1(n,n,n);
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            for (size_t k = 0; k < n; ++k) {
                v1(i,j,k) = static_cast<T>(i+j+k);
            }
        }
    }

    // Create a clone of the tensor.
    Tensor<T> v2 = std::move(v1.clone());

    // Check error
    assert(v1.size() == v2.size());
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            for (size_t k = 0; k < n; ++k) {
                double err = std::fabs(static_cast<double>(v2(i,j,k)) -
                                       static_cast<double>(v1(i,j,k)));
                if (err != 0.0) {
                    std::cerr << "err " << err
                        << " v1 " << v1(i,j,k)
                        << " v2 " << v2(i,j,k)
                        << "\n";
                    return false;
                }
           }
        }
    }
    return true;
}

#endif // TEST_NUMERIC_ARRAY3D_H_
