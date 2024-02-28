//
// test-array-templ.cpp
//
// Copyright (c) 2020 Carlos Braga
// This program is free software; you can redistribute it and/or modify it
// under the terms of the MIT License. See accompanying LICENSE.md or
// https://opensource.org/licenses/MIT.
///

#include "catch2/catch.hpp"
#include "algods/arrays/array.h"

/// ---- Test Array of strings ------------------------------------------------
TEST_CASE("Array of strings") {
    // Vector of strings //
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

    SECTION("Test Array API") {
        {
            //
            // Check initial state
            //
            Array<std::string> array_a;
            REQUIRE(array_a.is_empty() == true);
            REQUIRE(array_a.size() == 0);

            //
            // Check push_back()
            //
            for (auto it = text.begin(); it != text.end(); ++it) {
                array_a.push_back(*it);
            }
            REQUIRE(array_a.size() == text.size());

            //
            // Check array iterator
            //
            auto it = array_a.begin();
            auto jt = text.begin();
            while (it != array_a.end() && jt != text.end()) {
                REQUIRE(*it == *jt);
                ++it;
                ++jt;
            }

            //
            // Check array clear()
            //
            array_a.clear();
            REQUIRE(array_a.is_empty() == true);
            REQUIRE(array_a.size() == 0);
        }

        {
            //
            // Check initial state
            //
            Array<std::string> array_a;
            REQUIRE(array_a.is_empty() == true);
            REQUIRE(array_a.size() == 0);

            //
            // Check push_back() / pop_back
            //
            for (auto it = text.begin(); it != text.end(); ++it) {
                array_a.push_back(*it);
            }
            REQUIRE(array_a.size() == text.size());

            auto jt = text.rbegin();
            while (!array_a.is_empty() && jt != text.rend()) {
                REQUIRE(array_a.back() == *jt);
                array_a.pop_back();
                ++jt;
            }
            REQUIRE(array_a.is_empty() == true);
            REQUIRE(array_a.size() == 0);
        }

        {
            //
            // Check initial state
            //
            Array<std::string> array_a;
            REQUIRE(array_a.is_empty() == true);
            REQUIRE(array_a.size() == 0);

            //
            // Check Array copy constructors
            //
            for (auto it = text.begin(); it != text.end(); ++it) {
                array_a.push_back(*it);
            }

            {
                Array<std::string> array_b(array_a);
                REQUIRE(array_b.size() == array_a.size());

                auto ia = array_a.begin();
                auto ib = array_b.begin();

                while (ia != array_a.end() && ib != array_b.end()) {
                    REQUIRE(*ia == *ib);
                    ++ia;
                    ++ib;
                }
            }

            {
                Array<std::string> array_b = array_a;
                REQUIRE(array_b.size() == array_a.size());

                auto ia = array_a.begin();
                auto ib = array_b.begin();

                while (ia != array_a.end() && ib != array_b.end()) {
                    REQUIRE(*ia == *ib);
                    ++ia;
                    ++ib;
                }
            }

            {
                Array<std::string> array_b;
                array_b = array_a;
                REQUIRE(array_b.size() == array_a.size());

                auto ia = array_a.begin();
                auto ib = array_b.begin();

                while (ia != array_a.end() && ib != array_b.end()) {
                    REQUIRE(*ia == *ib);
                    ++ia;
                    ++ib;
                }
            }
        }
    }
}

/// ---- Test Array of Items --------------------------------------------------
struct Item {
    size_t v_;
    const size_t &v() { return v_; }

    Item() : v_(0) {}
    Item(size_t v) : v_(v) {}
    ~Item() = default;

    Item(const Item &other) { v_ = other.v_; }
    Item &operator=(const Item &other) {
        if (this == &other)
            return *this;
        v_ = other.v_;
        return *this;
    }

    Item(Item &&other) { v_ = std::move(other.v_); }
    Item &operator=(Item &&other) {
        if (this == &other)
            return *this;
        v_ = std::move(other.v_);
        return *this;
    }
};

TEST_CASE("Array of items") {
    // Vector of strings //
    std::vector<Item> vec_items;
    const size_t num_items = 2048;
    for (size_t i = 0; i < num_items; ++i) {
        vec_items.push_back(Item(i));
    }

    SECTION("Test Array API") {
        {
            //
            // Check initial state
            //
            Array<Item> array_a;
            REQUIRE(array_a.is_empty() == true);
            REQUIRE(array_a.size() == 0);

            //
            // Check push_back()
            //
            for (auto it = vec_items.begin(); it != vec_items.end(); ++it) {
                array_a.push_back(*it);
            }
            REQUIRE(array_a.size() == vec_items.size());

            //
            // Check array iterator
            //
            auto it = array_a.begin();
            auto jt = vec_items.begin();
            while (it != array_a.end() && jt != vec_items.end()) {
                REQUIRE(it->v_ == jt->v_);
                ++it;
                ++jt;
            }

            //
            // Check array clear()
            //
            array_a.clear();
            REQUIRE(array_a.is_empty() == true);
            REQUIRE(array_a.size() == 0);
        }

        {
            //
            // Check initial state
            //
            Array<Item> array_a;
            REQUIRE(array_a.is_empty() == true);
            REQUIRE(array_a.size() == 0);

            //
            // Check push_back()
            //
            for (auto it = vec_items.begin(); it != vec_items.end(); ++it) {
                array_a.push_back(*it);
            }
            REQUIRE(array_a.size() == vec_items.size());

            //
            // Check array iterator
            //
            auto it = array_a.begin();
            auto jt = array_a.rbegin();
            while (it != array_a.end() && jt != array_a.rend()) {
                REQUIRE(it->v_ + jt->v_ == array_a.size()-1);
                ++it;
                ++jt;
            }
        }

        {
            //
            // Check initial state
            //
            Array<Item> array_a;
            REQUIRE(array_a.is_empty() == true);
            REQUIRE(array_a.size() == 0);

            //
            // Check push_back()
            //
            array_a.clear();
            for (auto it = vec_items.begin(); it != vec_items.end(); ++it) {
                array_a.push_back(*it);
            }
            REQUIRE(array_a.size() == vec_items.size());

            //
            // Check array operator[]
            //
            for (size_t i = 0; i < array_a.size(); ++i) {
                REQUIRE(array_a[i].v_ == vec_items[i].v_);
            }
        }

        {
            //
            // Check size()/resize()/capacity
            //
            Array<Item> array_b;
            REQUIRE(array_b.data() != nullptr);
            REQUIRE(array_b.capacity() == 1);
            REQUIRE(array_b.size() == 0);

            array_b.push_back(Item(1));
            REQUIRE(array_b.capacity() == 1);
            REQUIRE(array_b.size() == 1);

            array_b.push_back(Item(2));
            REQUIRE(array_b.capacity() == 2);
            REQUIRE(array_b.size() == 2);

            array_b.push_back(Item(3));
            REQUIRE(array_b.capacity() == 4);
            REQUIRE(array_b.size() == 3);

            array_b.push_back(Item(4));
            REQUIRE(array_b.capacity() == 4);
            REQUIRE(array_b.size() == 4);

            array_b.push_back(Item(5));
            REQUIRE(array_b.capacity() == 8);
            REQUIRE(array_b.size() == 5);

            array_b.push_back(Item(6));
            array_b.push_back(Item(7));
            array_b.push_back(Item(8));

            array_b.resize(12);
            REQUIRE(array_b.capacity() == 12);
            REQUIRE(array_b.size() == 8);

            array_b.push_back(Item(9));
            array_b.push_back(Item(10));
            array_b.push_back(Item(11));
            array_b.push_back(Item(12));
            array_b.push_back(Item(13));
            array_b.push_back(Item(14));

            REQUIRE(array_b.capacity() == 24);
            REQUIRE(array_b.size() == 14);

            array_b.resize(14);
            REQUIRE(array_b.capacity() == 14);
            REQUIRE(array_b.size() == 14);

            array_b.push_back(Item(15));
            REQUIRE(array_b.capacity() == 28);
            REQUIRE(array_b.size() == 15);
        }
    }
}
