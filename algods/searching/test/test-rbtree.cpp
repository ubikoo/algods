//
// test-rbtree.cpp
//
// Copyright (c) 2020 Carlos Braga
// This program is free software; you can redistribute it and/or modify it
// under the terms of the MIT License. See accompanying LICENSE.md or
// https://opensource.org/licenses/MIT.
//

#include "catch2/catch.hpp"
#include "algods/searching/rbtree.h"
#include <map>
#include <vector>
#include <random>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <algorithm>

/// ---- Test red black tree data structure ------------------------------------
///
TEST_CASE("Red black tree") {
    SECTION("Test API") {
        //
        // Generate tree
        //
        std::vector<char> letters = {
            'S','E','A','R','C','H','E','X','A','M','P','L','E'};
        std::vector<char> letters_min = {
            'S','E','A','A','A','A','A','A','A','A','A','A','A'};
        std::vector<char> letters_max = {
            'S','S','S','S','S','S','S','X','X','X','X','X','X'};

        RBTree<char, size_t> tree;
        std::map<char, size_t> tree_map;

        REQUIRE(tree.is_empty() == true);
        REQUIRE(tree.size() == 0);
        for (size_t i = 0; i < letters.size(); ++i) {
            tree.put(letters[i], i);
            tree_map[letters[i]] = i;

            REQUIRE(tree.is_empty() == false);
            REQUIRE(tree.size() == tree_map.size());

            REQUIRE(tree.contains(letters[i]) == true);
            REQUIRE(tree.get(letters[i]) == i);

            REQUIRE(tree.min() == letters_min[i]);
            REQUIRE(tree.max() == letters_max[i]);
        }

        //
        //  Check tree
        //
        //  A 8
        //  C 4
        //  E 12
        //  H 5
        //  L 11
        //  M 9
        //  P 10
        //  R 3
        //  S 0
        //  X 7
        //
        for (auto k : tree.keys()) {
            REQUIRE(tree.contains(k) == true);
            REQUIRE(tree.get(k) == tree_map[k]);
        }

        //
        // Check tree floor, ceil, select and rank
        //
        REQUIRE(tree.floor('B') == 'A');
        REQUIRE(tree.ceil('B') == 'C');

        REQUIRE(tree.floor('K') == 'H');
        REQUIRE(tree.ceil('K') == 'L');

        REQUIRE(tree.floor('O') == 'M');
        REQUIRE(tree.ceil('O') == 'P');

        for (size_t k = 0; k < tree.size(); ++k) {
            REQUIRE(tree.rank(tree.select(k)) == k);
        }

        //
        // Print tree
        //
        std::cout << "\n-----------------------------------------\n";
        std::cout << tree.graph();
        std::cout << "\n-----------------------------------------\n";
        std::cout << tree.graphviz();
        std::cout << "\n-----------------------------------------\n";
    }

    SECTION("Test vector") {
        const size_t numiters = 10;

        //
        // Generate tree
        //
        std::vector<size_t> numbers(262144);
        std::iota(numbers.begin(), numbers.end(), 1);
        std::random_device rd;
        std::mt19937 rng(rd());

        RBTree<size_t, size_t> tree;
        std::map<size_t, size_t> tree_map;

        for (size_t iter = 0; iter < numiters; ++iter) {
            std::cout << "Red black tree iteration " << iter << "\n";

            std::shuffle(numbers.begin(), numbers.end(), rng);
            tree.clear();
            tree_map.clear();
            for (auto n : numbers) {
                tree.put(n, 2*n);
                tree_map[n] = 2*n;
                REQUIRE(tree.size() == tree_map.size());
            }

            //
            // Test tree
            //
            for (auto k : tree.keys()) {
                REQUIRE(tree.contains(k) == true);
                REQUIRE(tree.get(k) == tree_map[k]);
            }

            std::vector<size_t> keys = tree.keys();
            for (size_t i = 0; i < keys.size() - 1; ++i) {
                REQUIRE(keys[i] < keys[i+1]);
            }

            for (size_t k = 0; k < tree.size(); ++k) {
                REQUIRE(tree.rank(tree.select(k)) == k);
            }
        }
    }
}
