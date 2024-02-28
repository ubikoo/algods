//
// test-list.cpp
//
// Copyright (c) 2020 Carlos Braga
// This program is free software; you can redistribute it and/or modify it
// under the terms of the MIT License. See accompanying LICENSE.md or
// https://opensource.org/licenses/MIT.
//

#include "catch2/catch.hpp"
#include "algods/linked/list.h"

/// ---- Test Doubly Linked List -----------------------------------------------
///
TEST_CASE("Doubly Linked List") {
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

    SECTION("Check List API") {
        {
            //
            // Check initial state
            //
            List<std::string> list_a;
            REQUIRE(list_a.is_empty() == true);
            REQUIRE(list_a.size() == 0);

            //
            // Check push_back() and back()
            //
            for (auto it = text.begin(); it != text.end(); ++it) {
                list_a.push_back(*it);
                REQUIRE(list_a.back() == *it);
            }
            REQUIRE(list_a.size() == text.size());

            //
            // Check list iterator
            //
            auto it = list_a.begin();
            auto jt = text.begin();
            while (it != list_a.end() && jt != text.end()) {
                REQUIRE(*it == *jt);
                ++it;
                ++jt;
            }

            //
            // Check list clear()
            //
            list_a.clear();
            REQUIRE(list_a.is_empty() == true);
            REQUIRE(list_a.size() == 0);
        }

        {
            //
            // Check initial state
            //
            List<std::string> list_a;
            REQUIRE(list_a.is_empty() == true);
            REQUIRE(list_a.size() == 0);

            //
            // Check push_front() and front()
            //
            for (auto it = text.begin(); it != text.end(); ++it) {
                list_a.push_front(*it);
                REQUIRE(list_a.front() == *it);
            }
            REQUIRE(list_a.size() == text.size());

            //
            // Check list forward iterator
            //
            auto it1 = list_a.begin();
            auto jt1 = text.rbegin();    // reverse order
            while (it1 != list_a.end() && jt1 != text.rend()) {
                REQUIRE(*it1 == *jt1);
                ++it1;
                ++jt1;
            }

            //
            // Check list reverse iterator
            //
            auto it2 = list_a.rbegin();  // reverse order
            auto jt2 = text.begin();
            while (it2 != list_a.rend() && jt2 != text.end()) {
                REQUIRE(*it2 == *jt2);
                ++it2;
                ++jt2;
            }

            //
            // Check list clear()
            //
            list_a.clear();
            REQUIRE(list_a.is_empty() == true);
            REQUIRE(list_a.size() == 0);
        }

        {
            //
            // Check initial state
            //
            List<std::string> list_a;
            REQUIRE(list_a.is_empty() == true);
            REQUIRE(list_a.size() == 0);

            //
            // Check push_back() and pop_front()
            //
            for (auto it = text.begin(); it != text.end(); ++it) {
                list_a.push_back(*it);
            }

            size_t n = 0;
            while (!list_a.is_empty()) {
                REQUIRE(list_a.front() == text[n]);
                REQUIRE(list_a.size() == text.size()-n);
                list_a.pop_front();
                n++;
            }
            REQUIRE(list_a.is_empty() == true);
            REQUIRE(list_a.size() == 0);
        }

        {
            //
            // Check initial state
            //
            List<std::string> list_a;
            REQUIRE(list_a.is_empty() == true);
            REQUIRE(list_a.size() == 0);

            //
            // Check push_front() and pop_back()
            //
            for (auto it = text.begin(); it != text.end(); ++it) {
                list_a.push_front(*it);
            }

            size_t n = 0;
            while (!list_a.is_empty()) {
                REQUIRE(list_a.back() == text[n]);
                REQUIRE(list_a.size() == text.size()-n);
                list_a.pop_back();
                n++;
            }
            REQUIRE(list_a.is_empty() == true);
            REQUIRE(list_a.size() == 0);
        }

        {
            //
            // Check initial state
            //
            List<std::string> list_a;
            REQUIRE(list_a.is_empty() == true);
            REQUIRE(list_a.size() == 0);

            //
            // Check push_back() and push()
            //
            for (auto it = text.begin(); it != text.end(); ++it) {
                list_a.push_back(*it);
            }

            List<std::string> list_b;
            List<std::string> list_c;
            for (auto &it : list_a) {
                list_b.push_front(it);
                list_c.push_back(it);
            }
            REQUIRE(list_b.size() == list_a.size());
            REQUIRE(list_c.size() == list_a.size());

            auto ia = list_a.begin();
            auto ib = list_b.rbegin();
            auto ic = list_c.begin();
            while (ia != list_a.end() &&
                   ib != list_b.rend() &&
                   ic != list_c.end()) {
                REQUIRE(*ib == *ia);
                REQUIRE(*ic == *ia);
                ++ia;
                ++ib;
                ++ic;
            }
        }

        {
            //
            // Check initial state
            //
            List<std::string> list_a;
            REQUIRE(list_a.is_empty() == true);
            REQUIRE(list_a.size() == 0);

            //
            // Check List copy constructors
            //
            for (auto it = text.begin(); it != text.end(); ++it) {
                list_a.push_front(*it);
            }

            List<std::string> list_b(list_a);
            REQUIRE(list_b.size() == list_a.size());
            auto it = list_a.begin();
            auto jt = list_b.begin();
            while (it != list_a.end() && jt != list_b.end()) {
                REQUIRE(*it == *jt);
                ++it;
                ++jt;
            }

            List<std::string> list_c = list_a;
            REQUIRE(list_c.size() == list_a.size());
            it = list_a.begin();
            jt = list_c.begin();
            while (it != list_a.end() && jt != list_c.end()) {
                REQUIRE(*it == *jt);
                ++it;
                ++jt;
            }

            List<std::string> list_d;
            list_d = list_b;
            REQUIRE(list_d.size() == list_a.size());
            it = list_a.begin();
            jt = list_d.begin();
            while (it != list_a.end() && jt != list_d.end()) {
                REQUIRE(*it == *jt);
                ++it;
                ++jt;
            }
        }
    }
}
