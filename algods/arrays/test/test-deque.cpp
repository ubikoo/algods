//
// test-deque-templ.cpp
//
// Copyright (c) 2020 Carlos Braga
// This program is free software; you can redistribute it and/or modify it
// under the terms of the MIT License. See accompanying LICENSE.md or
// https://opensource.org/licenses/MIT.
///

#include "catch2/catch.hpp"
#include "algods/arrays/deque.h"

/// ---- Test Deque Data Structure --------------------------------------------
TEST_CASE("Deque of strings") {
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

    SECTION("Test Deque API") {
        {
            //
            // Check initial state
            //
            Deque<std::string> deque_a;
            REQUIRE(deque_a.is_empty() == true);
            REQUIRE(deque_a.size() == 0);

            //
            // Check push_back / forward iterator
            //
            for (auto it = text.begin(); it != text.end(); ++it) {
                deque_a.push_back(*it);
            }
            REQUIRE(deque_a.size() == text.size());

            auto it1 = deque_a.begin();
            auto jt1 = text.begin();
            while (it1 != deque_a.end() && jt1 != text.end()) {
                REQUIRE(*it1 == *jt1);
                ++it1;
                ++jt1;
            }

            //
            // Check deque clear
            //
            deque_a.clear();
            REQUIRE(deque_a.is_empty() == true);
            REQUIRE(deque_a.size() == 0);

            //
            // Check push_front / reverse iterator
            //
            for (auto it = text.begin(); it != text.end(); ++it) {
                deque_a.push_front(*it);
            }
            REQUIRE(deque_a.size() == text.size());

            auto it2 = deque_a.rbegin();
            auto jt2 = text.begin();
            while (it2 != deque_a.rend() && jt2 != text.end()) {
                REQUIRE(*it2 == *jt2);
                ++it2;
                ++jt2;
            }
        }

        {
            //
            // Check initial state
            //
            Deque<std::string> deque_a;
            REQUIRE(deque_a.is_empty() == true);
            REQUIRE(deque_a.size() == 0);

            //
            // Check push_back / pop_front
            //
            for (auto it = text.begin(); it != text.end(); ++it) {
                deque_a.push_back(*it);
            }
            REQUIRE(deque_a.size() == text.size());

            auto it = deque_a.begin();
            auto jt = text.begin();
            while (it != deque_a.end() && jt != text.end()) {
                REQUIRE(*it == *jt);
                ++it;
                ++jt;
            }

            auto kt = text.begin();
            while (!deque_a.is_empty() && kt != text.end()) {
                REQUIRE(deque_a.front() == *kt);
                deque_a.pop_front();
                ++kt;
            }
            REQUIRE(deque_a.is_empty() == true);
            REQUIRE(deque_a.size() == 0);
        }

        {
            //
            // Check initial state
            //
            Deque<std::string> deque_a;
            REQUIRE(deque_a.is_empty() == true);
            REQUIRE(deque_a.size() == 0);

            //
            // Check push_back()
            //
            for (auto it = text.begin(); it != text.end(); ++it) {
                deque_a.push_back(*it);
            }
            REQUIRE(deque_a.size() == text.size());

            //
            // Check deque iterator
            //
            auto it = deque_a.begin();
            auto jt = text.begin();
            while (it != deque_a.end() && jt != text.end()) {
                REQUIRE(*it == *jt);
                ++it;
                ++jt;
            }

            //
            // Check pop_front
            //
            auto kt = text.rbegin();
            while (!deque_a.is_empty() && kt != text.rend()) {
                REQUIRE(deque_a.back() == *kt);
                deque_a.pop_back();
                ++kt;
            }
            REQUIRE(deque_a.is_empty() == true);
            REQUIRE(deque_a.size() == 0);
        }

        {
            //
            // Check initial state
            //
            Deque<std::string> deque_a;
            REQUIRE(deque_a.is_empty() == true);
            REQUIRE(deque_a.size() == 0);

            //
            // Check push_back()
            //
            for (auto it = text.begin(); it != text.end(); ++it) {
                deque_a.push_back(*it);
            }
            REQUIRE(deque_a.size() == text.size());

            //
            // Check deque iterator
            //
            auto it = deque_a.begin();
            auto jt = text.begin();
            while (it != deque_a.end() && jt != text.end()) {
                REQUIRE(*it == *jt);
                ++it;
                ++jt;
            }

            //
            // Check pop_front
            //
            auto kt = text.rbegin();
            while (!deque_a.is_empty() && kt != text.rend()) {
                REQUIRE(deque_a.back() == *kt);
                deque_a.pop_back();
                ++kt;
            }
            REQUIRE(deque_a.is_empty() == true);
            REQUIRE(deque_a.size() == 0);
        }

        {
            //
            // Check initial state
            //
            Deque<std::string> deque_a;
            REQUIRE(deque_a.is_empty() == true);
            REQUIRE(deque_a.size() == 0);

            //
            // Check Deque copy constructors
            //
            for (auto it = text.begin(); it != text.end(); ++it) {
                deque_a.push_back(*it);
            }

            {
                Deque<std::string> deque_b(deque_a);
                REQUIRE(deque_b.size() == deque_a.size());

                auto ia = deque_a.begin();
                auto ib = deque_b.begin();

                while (ia != deque_a.end() && ib != deque_b.end()) {
                    REQUIRE(*ia == *ib);
                    ++ia;
                    ++ib;
                }
            }

            {
                Deque<std::string> deque_b = deque_a;
                REQUIRE(deque_b.size() == deque_a.size());

                auto ia = deque_a.begin();
                auto ib = deque_b.begin();

                while (ia != deque_a.end() && ib != deque_b.end()) {
                    REQUIRE(*ia == *ib);
                    ++ia;
                    ++ib;
                }
            }

            {
                Deque<std::string> deque_b;
                deque_b = deque_a;
                REQUIRE(deque_b.size() == deque_a.size());

                auto ia = deque_a.begin();
                auto ib = deque_b.begin();

                while (ia != deque_a.end() && ib != deque_b.end()) {
                    REQUIRE(*ia == *ib);
                    ++ia;
                    ++ib;
                }
            }
        }
    }
}

/// ---- Test Deque Data Structure --------------------------------------------
struct Item {
    // Member variables/accessors
    size_t v_;
    const size_t &v() { return v_; }

    // Constructor/destructor
    Item() : v_(0) {}
    Item(size_t v) : v_(v) {}
    ~Item() = default;

    // Copy constructor/assignment
    Item(const Item &other) { v_ = other.v_; }
    Item &operator=(const Item &other) {
        if (this == &other)
            return *this;
        v_ = other.v_;
        return *this;
    }

    // Move constructor/assignment
    Item(Item &&other) { v_ = std::move(other.v_); }
    Item &operator=(Item &&other) {
        if (this == &other)
            return *this;
        v_ = std::move(other.v_);
        return *this;
    }
};


TEST_CASE("Deque of items") {
    // Vector of strings
    std::vector<Item> vec_items;
    const size_t num_items = 4096;
    for (size_t i = 0; i < num_items; ++i) {
        vec_items.push_back(Item(i));
    }

    SECTION("Test Deque API") {
        {
            //
            // Check push_back()
            //
            Deque<Item> deque_a;
            REQUIRE(deque_a.is_empty() == true);
            REQUIRE(deque_a.capacity() == 1);
            REQUIRE(deque_a.size() == 0);
            for (auto it = vec_items.begin(); it != vec_items.end(); ++it) {
                deque_a.push_back(*it);
            }
            REQUIRE(deque_a.size() == vec_items.size());

            //
            // Check deque forward iterator
            //
            auto it1 = deque_a.begin();
            auto jt1 = vec_items.begin();
            while (it1 != deque_a.end() && jt1 != vec_items.end()) {
                REQUIRE(it1->v_ == jt1->v_);
                ++it1;
                ++jt1;
            }

            //
            // Check deque reverse iterator
            //
            auto it2 = deque_a.rbegin();
            auto jt2 = vec_items.rbegin();
            while (it2 != deque_a.rend() && jt2 != vec_items.rend()) {
                REQUIRE(it2->v_ == jt2->v_);
                ++it2;
                ++jt2;
            }

            //
            // Check deque operator[]
            //
            for (size_t i = 0; i < deque_a.size(); ++i) {
                REQUIRE(deque_a[i].v_ == vec_items[i].v_);
            }

            //
            // Check deque clear()
            //
            deque_a.clear();
            REQUIRE(deque_a.is_empty() == true);
            REQUIRE(deque_a.size() == 0);
        }

        {
            //
            // Check transfer
            //
            Deque<Item> deque_a;
            REQUIRE(deque_a.is_empty() == true);
            REQUIRE(deque_a.capacity() == 1);
            REQUIRE(deque_a.size() == 0);
            for (auto it = vec_items.begin(); it != vec_items.end(); ++it) {
                deque_a.push_front(*it);
            }
            REQUIRE(deque_a.size() == vec_items.size());

            Deque<Item> deque_b;
            Deque<Item> deque_c;
            for (auto item : deque_a) {
                deque_b.push_back(item);
                deque_c.push_front(item);
            }

            auto ia = deque_a.begin();
            auto ib = deque_b.begin();
            auto ic = deque_c.rbegin();
            while (ia != deque_a.end() &&
                   ib != deque_b.end() &&
                   ic != deque_c.rend()) {
                REQUIRE(ia->v_ == ib->v_);
                REQUIRE(ia->v_ == ic->v_);
            }

            auto ja = deque_a.rbegin();
            auto jb = deque_b.rbegin();
            auto jc = deque_c.begin();
            while (ja != deque_a.rend() &&
                   jb != deque_b.rend() &&
                   jc != deque_c.end()) {
                REQUIRE(ja->v_ == jb->v_);
                REQUIRE(ja->v_ == jc->v_);
            }
        }

        {
            //
            // Check transfer
            //
            Deque<Item> deque_a;
            REQUIRE(deque_a.is_empty() == true);
            REQUIRE(deque_a.capacity() == 1);
            REQUIRE(deque_a.size() == 0);

            deque_a.push_back(Item(1));
            REQUIRE(deque_a.capacity() == 1);
            REQUIRE(deque_a.size() == 1);

            deque_a.push_front(Item(2));
            REQUIRE(deque_a.capacity() == 2);
            REQUIRE(deque_a.size() == 2);

            REQUIRE(deque_a.back().v_ == 1);
            REQUIRE(deque_a.front().v_ == 2);

            deque_a.push_front(Item(3));
            REQUIRE(deque_a.capacity() == 4);
            REQUIRE(deque_a.size() == 3);

            deque_a.push_front(Item(4));
            REQUIRE(deque_a.capacity() == 4);
            REQUIRE(deque_a.size() == 4);

            REQUIRE(deque_a[0].v_ == 4);
            REQUIRE(deque_a[1].v_ == 3);
            REQUIRE(deque_a[2].v_ == 2);
            REQUIRE(deque_a[3].v_ == 1);

            deque_a.push_back(Item(5));
            REQUIRE(deque_a.capacity() == 8);
            REQUIRE(deque_a.size() == 5);

            REQUIRE(deque_a[0].v_ == 4);
            REQUIRE(deque_a[1].v_ == 3);
            REQUIRE(deque_a[2].v_ == 2);
            REQUIRE(deque_a[3].v_ == 1);
            REQUIRE(deque_a[4].v_ == 5);

            REQUIRE(deque_a.front().v_ == 4);
            REQUIRE(deque_a.back().v_ == 5);

            deque_a.resize(5);
            REQUIRE(deque_a.capacity() == 5);
            REQUIRE(deque_a.size() == 5);

            deque_a.push_back(Item(6));
            REQUIRE(deque_a.capacity() == 10);
            REQUIRE(deque_a.size() == 6);

            REQUIRE(deque_a.back().v_ == 6);
            deque_a.pop_back();
            REQUIRE(deque_a.back().v_ == 5);
            deque_a.pop_back();
            REQUIRE(deque_a.back().v_ == 1);
            deque_a.pop_back();
            REQUIRE(deque_a.back().v_ == 2);
            deque_a.pop_back();

            REQUIRE(deque_a.back().v_ == 3);
            deque_a.pop_back();

            REQUIRE(deque_a.capacity() == 2);
            REQUIRE(deque_a.size() == 1);

            REQUIRE(deque_a.back().v_ == 4);
            deque_a.pop_back();

            REQUIRE(deque_a.capacity() == 2);
            REQUIRE(deque_a.size() == 0);
        }
    }
}
