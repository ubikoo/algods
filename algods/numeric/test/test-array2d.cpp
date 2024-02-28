//
// test-array2d.cpp
//
// Copyright (c) 2020 Carlos Braga
// This program is free software; you can redistribute it and/or modify it
// under the terms of the MIT License. See accompanying LICENSE.md or
// https://opensource.org/licenses/MIT.
//

#include "catch2/catch.hpp"
#include "test-array2d.h"

TEST_CASE("Matrix") {
    SECTION("Test Matrix API") {
        const size_t n = 8;            // minimum matrix size

        // Test matrix read/write
        REQUIRE(TestMatrixReadWrite<long double>(n*n*n));
        REQUIRE(TestMatrixReadWrite<double>(n*n*n));
        REQUIRE(TestMatrixReadWrite<float>(n*n));

        REQUIRE(TestMatrixReadWrite<long long>(n*n*n));
        REQUIRE(TestMatrixReadWrite<long>(n*n*n));
        REQUIRE(TestMatrixReadWrite<int>(n*n));
        REQUIRE(TestMatrixReadWrite<short>(n));

        REQUIRE(TestMatrixReadWrite<unsigned long long>(n*n*n));
        REQUIRE(TestMatrixReadWrite<unsigned long>(n*n*n));
        REQUIRE(TestMatrixReadWrite<unsigned int>(n*n));
        REQUIRE(TestMatrixReadWrite<unsigned short>(n));

        REQUIRE(TestMatrixReadWrite<int64_t>(n*n*n));
        REQUIRE(TestMatrixReadWrite<int32_t>(n*n));
        REQUIRE(TestMatrixReadWrite<int16_t>(n));

        REQUIRE(TestMatrixReadWrite<uint64_t>(n*n*n));
        REQUIRE(TestMatrixReadWrite<uint32_t>(n*n));
        REQUIRE(TestMatrixReadWrite<uint16_t>(n));

        // Test matrix copy/assign
        REQUIRE(TestMatrixClone<long double>(n*n*n));
        REQUIRE(TestMatrixClone<double>(n*n*n));
        REQUIRE(TestMatrixClone<float>(n*n));

        REQUIRE(TestMatrixClone<long long>(n*n*n));
        REQUIRE(TestMatrixClone<long>(n*n*n));
        REQUIRE(TestMatrixClone<int>(n*n));
        REQUIRE(TestMatrixClone<short>(n));
        REQUIRE(TestMatrixClone<char>(n));

        REQUIRE(TestMatrixClone<unsigned long long>(n*n*n));
        REQUIRE(TestMatrixClone<unsigned long>(n*n*n));
        REQUIRE(TestMatrixClone<unsigned int>(n*n));
        REQUIRE(TestMatrixClone<unsigned short>(n));
        REQUIRE(TestMatrixClone<unsigned char>(n));

        REQUIRE(TestMatrixClone<int64_t>(n*n*n));
        REQUIRE(TestMatrixClone<int32_t>(n*n));
        REQUIRE(TestMatrixClone<int16_t>(n));
        REQUIRE(TestMatrixClone<int8_t>(n));

        REQUIRE(TestMatrixClone<uint64_t>(n*n*n));
        REQUIRE(TestMatrixClone<uint32_t>(n*n));
        REQUIRE(TestMatrixClone<uint16_t>(n));
        REQUIRE(TestMatrixClone<uint8_t>(n));
    }
}
