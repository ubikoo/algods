//
// test-array3d.cpp
//
// Copyright (c) 2020 Carlos Braga
// This program is free software; you can redistribute it and/or modify it
// under the terms of the MIT License. See accompanying LICENSE.md or
// https://opensource.org/licenses/MIT.
//

#include "catch2/catch.hpp"
#include "test-array3d.h"

TEST_CASE("Tensor") {
    SECTION("Test Tensor API") {
        const size_t n = 4;             // minimum tensor size

        // Test tensor read/write
        REQUIRE(TestTensorReadWrite<long double>(n*n*n));
        REQUIRE(TestTensorReadWrite<double>(n*n*n));
        REQUIRE(TestTensorReadWrite<float>(n*n));

        REQUIRE(TestTensorReadWrite<long long>(n*n*n));
        REQUIRE(TestTensorReadWrite<long>(n*n*n));
        REQUIRE(TestTensorReadWrite<int>(n*n));
        REQUIRE(TestTensorReadWrite<short>(n));

        REQUIRE(TestTensorReadWrite<unsigned long long>(n*n*n));
        REQUIRE(TestTensorReadWrite<unsigned long>(n*n*n));
        REQUIRE(TestTensorReadWrite<unsigned int>(n*n));
        REQUIRE(TestTensorReadWrite<unsigned short>(n));

        REQUIRE(TestTensorReadWrite<int64_t>(n*n*n));
        REQUIRE(TestTensorReadWrite<int32_t>(n*n));
        REQUIRE(TestTensorReadWrite<int16_t>(n));

        REQUIRE(TestTensorReadWrite<uint64_t>(n*n*n));
        REQUIRE(TestTensorReadWrite<uint32_t>(n*n));
        REQUIRE(TestTensorReadWrite<uint16_t>(n));

        // Test tensor copy/assign
        REQUIRE(TestTensorClone<long double>(n*n*n));
        REQUIRE(TestTensorClone<double>(n*n*n));
        REQUIRE(TestTensorClone<float>(n*n));

        REQUIRE(TestTensorClone<long long>(n*n*n));
        REQUIRE(TestTensorClone<long>(n*n*n));
        REQUIRE(TestTensorClone<int>(n*n));
        REQUIRE(TestTensorClone<short>(n));
        REQUIRE(TestTensorClone<char>(n));

        REQUIRE(TestTensorClone<unsigned long long>(n*n*n));
        REQUIRE(TestTensorClone<unsigned long>(n*n*n));
        REQUIRE(TestTensorClone<unsigned int>(n*n));
        REQUIRE(TestTensorClone<unsigned short>(n));
        REQUIRE(TestTensorClone<unsigned char>(n));

        REQUIRE(TestTensorClone<int64_t>(n*n*n));
        REQUIRE(TestTensorClone<int32_t>(n*n));
        REQUIRE(TestTensorClone<int16_t>(n));
        REQUIRE(TestTensorClone<int8_t>(n));

        REQUIRE(TestTensorClone<uint64_t>(n*n*n));
        REQUIRE(TestTensorClone<uint32_t>(n*n));
        REQUIRE(TestTensorClone<uint16_t>(n));
        REQUIRE(TestTensorClone<uint8_t>(n));
    }
}
