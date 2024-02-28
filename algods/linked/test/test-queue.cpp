//
// test-queue.cpp
//
// Copyright (c) 2020 Carlos Braga
// This program is free software; you can redistribute it and/or modify it
// under the terms of the MIT License. See accompanying LICENSE.md or
// https://opensource.org/licenses/MIT.
//

#include "catch2/catch.hpp"
#include "algods/linked/queue.h"

/// ---- Test Queue Data Structure ---------------------------------------------
///
TEST_CASE("Queue") {
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

    SECTION("Test Queue API") {
        {
            //
            // Check initial state
            //
            Queue<std::string> queue_a;
            REQUIRE(queue_a.is_empty() == true);
            REQUIRE(queue_a.size() == 0);

            //
            // Check add()
            //
            for (auto it = text.begin(); it != text.end(); ++it) {
                queue_a.enqueue(*it);
            }
            REQUIRE(queue_a.size() == text.size());

            //
            // Check queue iterator
            //
            auto it = queue_a.begin();
            auto jt = text.begin();
            while (it != queue_a.end() && jt != text.end()) {
                REQUIRE(*it == *jt);
                ++it;
                ++jt;
            }

            //
            // Check queue clear()
            //
            queue_a.clear();
            REQUIRE(queue_a.is_empty() == true);
            REQUIRE(queue_a.size() == 0);
        }

        {
            //
            // Check initial state
            //
            Queue<std::string> queue_a;
            REQUIRE(queue_a.is_empty() == true);
            REQUIRE(queue_a.size() == 0);

            //
            // Check Queue copy constructors
            //
            for (auto it = text.begin(); it != text.end(); ++it) {
                queue_a.enqueue(*it);
            }

            {
                Queue<std::string> queue_b(queue_a);
                REQUIRE(queue_b.size() == queue_a.size());
                auto it = queue_a.begin();
                auto jt = queue_b.begin();
                while (it != queue_a.end() && jt != queue_b.end()) {
                    REQUIRE(*it == *jt);
                    ++it;
                    ++jt;
                }
            }

            {
                Queue<std::string> queue_b = queue_a;
                REQUIRE(queue_b.size() == queue_a.size());
                auto it = queue_a.begin();
                auto jt = queue_b.begin();
                while (it != queue_a.end() && jt != queue_b.end()) {
                    REQUIRE(*it == *jt);
                    ++it;
                    ++jt;
                }
            }

            {
                Queue<std::string> queue_b;
                queue_b = queue_a;
                REQUIRE(queue_b.size() == queue_a.size());
                auto it = queue_a.begin();
                auto jt = queue_b.begin();
                while (it != queue_a.end() && jt != queue_b.end()) {
                    REQUIRE(*it == *jt);
                    ++it;
                    ++jt;
                }
            }
        }
    }
}
