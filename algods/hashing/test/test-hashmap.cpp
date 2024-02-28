//
// test-hashmap.cpp
//
// Copyright (c) 2020 Carlos Braga
// This program is free software; you can redistribute it and/or modify it
// under the terms of the MIT License. See accompanying LICENSE.md or
// https://opensource.org/licenses/MIT.
//

#include "catch2/catch.hpp"
#include "algods/hashing/hashmap.h"
#include <stddef.h> // size_t

/// ---- Test hashmap data structure -------------------------------------------
struct Point {
    // Member variables/accessors
    int x_ = 0;
    int y_ = 0;
    int z_ = 0;

    // Constructor/destructor
    Point() = default;
    Point(int x, int y, int z) : x_(x), y_(y), z_(z) {}
    ~Point() = default;

    bool operator== (const Point &other) {
        if (x_ != other.x_)
            return false;
        if (y_ != other.y_)
            return false;
        if (z_ != other.z_)
            return false;
        return true;
    }
};

///
/// Key hash functor.
///
template<typename Key>
struct Hasher {
    size_t operator()(const Key &key) {
        constexpr uint64_t mul = (uint64_t) 0x517cc1b727220bff;
        uint64_t h = static_cast<uint64_t>(key);
        return h * mul;
    }
};

///
/// Point hash template specialization.
///
template<>
struct Hasher<Point> {
    size_t operator()(const Point &pt) {
        return 541*pt.x_ + 79*pt.y_ + 31*pt.z_;
    }
};


TEST_CASE("Hashmap") {
    const size_t num_keys = 8000;
    const size_t dim_x = 20;
    const size_t dim_y = 20;
    const size_t dim_z = 20;

    // ----------------------------------------------------------------------
    SECTION("Test API with default hasher") {
        //
        // Generate a vector of value integers
        //
        std::vector<size_t> numbers(num_keys);
        std::iota(numbers.begin(), numbers.end(), 1);
        std::random_device rd;
        std::mt19937 rng(rd());
        std::shuffle(numbers.begin(), numbers.end(), rng);

        //
        // Generate a hashmap from the vector keying each
        // integer value by its index.
        //
        Hashmap<size_t, size_t> hashmap;
        REQUIRE(hashmap.is_empty());
        for (size_t i = 0; i < numbers.size(); ++i) {
            hashmap.put(i, numbers[i]);   // insert the(i,numbers[i]) pair
            REQUIRE(hashmap.contains(i));
            REQUIRE(hashmap.get(i) == numbers[i]);
            REQUIRE(hashmap.size() == i+1);
        }
        REQUIRE(!hashmap.is_empty());

        //
        // Print the hashmap
        //
        hashmap.print_chainloads();
#if 0
        for (auto &item : hashmap.keys()) {
            std::cout << "key " << item.first
                      << " value " << item.second
                      << std::endl;
        }
#endif

        //
        // Remove the keys from the hash map
        //
        for (size_t i = 0; i < numbers.size(); ++i) {
            hashmap.erase(i);
            REQUIRE(!hashmap.contains(i));
            REQUIRE(hashmap.size() == numbers.size()-(i+1));
        }
        REQUIRE(hashmap.is_empty());
    }


    // ----------------------------------------------------------------------
    SECTION("Test API with custom hasher") {
        //
        // Generate a vector of value integers
        //
        std::vector<size_t> numbers(num_keys);
        std::iota(numbers.begin(), numbers.end(), 1);
        std::random_device rd;
        std::mt19937 rng(rd());
        std::shuffle(numbers.begin(), numbers.end(), rng);

        //
        // Generate a hashmap from the vector keying each
        // integer value by its index.
        //
        Hashmap<size_t, size_t, Hasher<size_t>> hashmap;
        REQUIRE(hashmap.is_empty());
        for (size_t i = 0; i < numbers.size(); ++i) {
            hashmap.put(i, numbers[i]);   // insert the(i,numbers[i]) pair
            REQUIRE(hashmap.contains(i));
            REQUIRE(hashmap.get(i) == numbers[i]);
            REQUIRE(hashmap.size() == i+1);
        }
        REQUIRE(!hashmap.is_empty());

        //
        // Print the hashmap
        //
        hashmap.print_chainloads();
#if 0
        for (auto &item : hashmap.keys()) {
            std::cout << "key " << item.first
                      << " value " << item.second
                      << std::endl;
        }
#endif

        //
        // Remove the keys from the hash map
        //
        for (size_t i = 0; i < numbers.size(); ++i) {
            hashmap.erase(i);
            REQUIRE(!hashmap.contains(i));
            REQUIRE(hashmap.size() == numbers.size()-(i+1));
        }
        REQUIRE(hashmap.is_empty());
    }

    // ----------------------------------------------------------------------
    SECTION("Test API with a point class") {
        //
        // Generate a vector of value integers
        //
        std::vector<Point> numbers;
        for (size_t x = 0; x < dim_x; ++x) {
            for (size_t y = 0; y < dim_y; ++y) {
                for (size_t z = 0; z < dim_z; ++z) {
                    numbers.emplace_back(x, y, z);
                }
            }
        }

        //
        // Generate a hashmap from the vector keying each
        // integer value by its corresponding point.
        //
        Hashmap<Point, size_t, Hasher<Point>> hashmap;
        REQUIRE(hashmap.is_empty());
        for (size_t i = 0; i < numbers.size(); ++i) {
            hashmap.put(numbers[i], i);   // insert the(numbers[i], i) pair
            REQUIRE(hashmap.contains(numbers[i]));
            REQUIRE(hashmap.get(numbers[i]) == i);
            REQUIRE(hashmap.size() == i+1);
        }
        REQUIRE(!hashmap.is_empty());

        //
        // Print the hashmap
        //
        hashmap.print_chainloads();
#if 0
        for (auto &item : hashmap.keys()) {
            std::cout << "key " << item.first.x_
                      <<    " " << item.first.y_
                      <<    " " << item.first.z_
                      << " value " << item.second
                      << std::endl;
        }
#endif

        //
        // Remove the keys from the hash map
        //
        for (size_t i = 0; i < numbers.size(); ++i) {
            hashmap.erase(numbers[i]);
            REQUIRE(!hashmap.contains(numbers[i]));
            REQUIRE(hashmap.size() == numbers.size()-(i+1));
        }
        REQUIRE(hashmap.is_empty());
    }
}
