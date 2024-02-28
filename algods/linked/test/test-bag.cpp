//
// test-bag.cpp
//
// Copyright (c) 2020 Carlos Braga
// This program is free software; you can redistribute it and/or modify it
// under the terms of the MIT License. See accompanying LICENSE.md or
// https://opensource.org/licenses/MIT.
//

#include "catch2/catch.hpp"
#include "algods/linked/bag.h"

/// ---- Test Bag Data Structure -----------------------------------------------
///
TEST_CASE("Bag") {
    // Vector of strings
    std::vector<std::string> text {
        "Lorem",
        "ipsum",
        "dolor",
        "sit",
        "amet,",
        "consectetur",
        "adipiscing",
        "elit.",
        "Morbi",
        "sodales",
        "eleifend",
        "magna,",
        "in",
        "laoreet",
        "sapien."};

    SECTION("Test Bag API") {
        {
            //
            // Check initial state
            //
            Bag<std::string> bag_a;
            REQUIRE(bag_a.is_empty() == true);
            REQUIRE(bag_a.size() == 0);

            //
            // Check add()
            //
            for (auto it = text.begin(); it != text.end(); ++it) {
                bag_a.add(*it);
            }
            REQUIRE(bag_a.size() == text.size());

            //
            // Check bag iterator
            //
            auto it = bag_a.begin();
            auto jt = text.rbegin();
            while (it != bag_a.end() && jt != text.rend()) {
                REQUIRE(*it == *jt);
                ++it;
                ++jt;
            }

            //
            // Check bag clear()
            //
            bag_a.clear();
            REQUIRE(bag_a.is_empty() == true);
            REQUIRE(bag_a.size() == 0);
        }

        {
            //
            // Check initial state
            //
            Bag<std::string> bag_a;
            REQUIRE(bag_a.is_empty() == true);
            REQUIRE(bag_a.size() == 0);

            //
            // Check Bag copy constructors
            //
            for (auto it = text.begin(); it != text.end(); ++it) {
                bag_a.add(*it);
            }

            {
                Bag<std::string> bag_b(bag_a);
                REQUIRE(bag_b.size() == bag_a.size());

                auto ia = bag_a.begin();
                auto ib = bag_b.begin();

                while (ia != bag_a.end() && ib != bag_b.end()) {
                    REQUIRE(*ia == *ib);
                    ++ia;
                    ++ib;
                }
            }

            {
                Bag<std::string> bag_b = bag_a;
                REQUIRE(bag_b.size() == bag_a.size());

                auto ia = bag_a.begin();
                auto ib = bag_b.begin();

                while (ia != bag_a.end() && ib != bag_b.end()) {
                    REQUIRE(*ia == *ib);
                    ++ia;
                    ++ib;
                }
            }

            {
                Bag<std::string> bag_b;
                bag_b = bag_a;
                REQUIRE(bag_b.size() == bag_a.size());

                auto ia = bag_a.begin();
                auto ib = bag_b.begin();

                while (ia != bag_a.end() && ib != bag_b.end()) {
                    REQUIRE(*ia == *ib);
                    ++ia;
                    ++ib;
                }
            }
        }
    }
}
