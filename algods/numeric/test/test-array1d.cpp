//
// test-array1d.cpp
//
// Copyright (c) 2020 Carlos Braga
// This program is free software; you can redistribute it and/or modify it
// under the terms of the MIT License. See accompanying LICENSE.md or
// https://opensource.org/licenses/MIT.
//

#include "catch2/catch.hpp"
#include "test-array1d.h"

TEST_CASE("Vector") {
    SECTION("Test Vector API") {
        const size_t n = 24;            // minimum vector size

        // Test vector read/write
        REQUIRE(TestVectorReadWrite<long double>(n*n*n));
        REQUIRE(TestVectorReadWrite<double>(n*n*n));
        REQUIRE(TestVectorReadWrite<float>(n*n));

        REQUIRE(TestVectorReadWrite<long long>(n*n*n));
        REQUIRE(TestVectorReadWrite<long>(n*n*n));
        REQUIRE(TestVectorReadWrite<int>(n*n));
        REQUIRE(TestVectorReadWrite<short>(n));

        REQUIRE(TestVectorReadWrite<unsigned long long>(n*n*n));
        REQUIRE(TestVectorReadWrite<unsigned long>(n*n*n));
        REQUIRE(TestVectorReadWrite<unsigned int>(n*n));
        REQUIRE(TestVectorReadWrite<unsigned short>(n));

        REQUIRE(TestVectorReadWrite<int64_t>(n*n*n));
        REQUIRE(TestVectorReadWrite<int32_t>(n*n));
        REQUIRE(TestVectorReadWrite<int16_t>(n));

        REQUIRE(TestVectorReadWrite<uint64_t>(n*n*n));
        REQUIRE(TestVectorReadWrite<uint32_t>(n*n));
        REQUIRE(TestVectorReadWrite<uint16_t>(n));

        // Test vector copy/assign
        REQUIRE(TestVectorClone<long double>(n*n*n));
        REQUIRE(TestVectorClone<double>(n*n*n));
        REQUIRE(TestVectorClone<float>(n*n));

        REQUIRE(TestVectorClone<long long>(n*n*n));
        REQUIRE(TestVectorClone<long>(n*n*n));
        REQUIRE(TestVectorClone<int>(n*n));
        REQUIRE(TestVectorClone<short>(n));
        REQUIRE(TestVectorClone<char>(n));

        REQUIRE(TestVectorClone<unsigned long long>(n*n*n));
        REQUIRE(TestVectorClone<unsigned long>(n*n*n));
        REQUIRE(TestVectorClone<unsigned int>(n*n));
        REQUIRE(TestVectorClone<unsigned short>(n));
        REQUIRE(TestVectorClone<unsigned char>(n));

        REQUIRE(TestVectorClone<int64_t>(n*n*n));
        REQUIRE(TestVectorClone<int32_t>(n*n));
        REQUIRE(TestVectorClone<int16_t>(n));
        REQUIRE(TestVectorClone<int8_t>(n));

        REQUIRE(TestVectorClone<uint64_t>(n*n*n));
        REQUIRE(TestVectorClone<uint32_t>(n*n));
        REQUIRE(TestVectorClone<uint16_t>(n));
        REQUIRE(TestVectorClone<uint8_t>(n));
    }
}
