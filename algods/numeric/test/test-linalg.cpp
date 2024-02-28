//
// test-linalg.cpp
//
// Copyright (c) 2020 Carlos Braga
// This program is free software; you can redistribute it and/or modify it
// under the terms of the MIT License. See accompanying LICENSE.md or
// https://opensource.org/licenses/MIT.
//

#include "catch2/catch.hpp"
#include <array>
#include <iostream>
#include "test-linalg.h"

TEST_CASE("Linear Algebra") {
    const size_t n = 5;
    const std::array<size_t, n> ndim = {256, 512, 1024, 2048, 4096};
    const std::array<size_t, n> jdim = {64, 128, 256, 512, 1024};
    const std::array<size_t, n> iter = {64, 16, 8, 2, 1};

    // Test Gauss elimination solver
    SECTION("Gauss solver") {
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < iter[i]; ++j) {
                std::cout << "TestLinalgGauss<double> dim = "
                    << ndim[i] << ", run " << iter[i] << "\n";
                TestLinalgGauss<double>(ndim[i]);
            }
        }

        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < iter[i]; ++j) {
                std::cout << "TestLinalgGauss<float> dim = "
                    << ndim[i] << ", run " << iter[i] << "\n";
                TestLinalgGauss<float>(ndim[i]);
            }
        }
    }

    // Test Jacobi eigenvalue solver
    SECTION("Jacobi solver") {
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < iter[i]; ++j) {
                std::cout << "TestLinalgJacobi<double> dim = "
                     << jdim[i] << ", run " << iter[i] << "\n";
                TestLinalgJacobi<double>(jdim[i]);
            }
        }

        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < iter[i]; ++j) {
                std::cout << "TestLinalgJacobi<float> dim = "
                    << jdim[i] << ", run " << iter[i] << "\n";
                TestLinalgJacobi<float>(jdim[i]);
            }
        }
    }
}
