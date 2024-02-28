//
// test-priority-queue.cpp
//
// Copyright (c) 2020 Carlos Braga
// This program is free software; you can redistribute it and/or modify it
// under the terms of the MIT License. See accompanying LICENSE.md or
// https://opensource.org/licenses/MIT.
//

#include "catch2/catch.hpp"
#include <iostream>
#include "test-priority-queue.h"

/// ---- Test PriorityQueue ----------------------------------------------------
TEST_CASE("PriorityQueue") {
    const size_t num = 100000;

    // Maximum priority queue
    SECTION("Test MaxPriorityQueue") {
        std::cout << "Test MaxPriorityQueue\n";

        // Create a min priority queue of those values
        PriorityQueue<int, std::greater<int>> queue;
        REQUIRE(queue.size() == 0);
        REQUIRE(queue.is_empty());

        queue.push(1);
        print_queue("push 1:", queue);
        REQUIRE(queue.size() > 0);
        REQUIRE(!queue.is_empty());
        REQUIRE(queue.top() == 1);

        queue.push(2);
        print_queue("push 2:", queue);
        REQUIRE(queue.top() == 2);

        queue.push(3);
        print_queue("push 3:", queue);
        REQUIRE(queue.top() == 3);

        queue.pop();
        print_queue("pop:", queue);
        REQUIRE(queue.top() == 2);

        queue.push(4);
        print_queue("push 4:", queue);
        REQUIRE(queue.top() == 4);
    }

    // Minimum priority queue
    SECTION("Test MinPriorityQueue") {
        std::cout << "Test MinPriorityQueue\n";

        // Create a min priority queue of those values
        PriorityQueue<int, std::less<int>> queue;
        REQUIRE(queue.size() == 0);
        REQUIRE(queue.is_empty());

        queue.push(4);
        print_queue("push 4:", queue);
        REQUIRE(queue.size() > 0);
        REQUIRE(!queue.is_empty());
        REQUIRE(queue.top() == 4);

        queue.push(3);
        print_queue("push 3:", queue);
        REQUIRE(queue.top() == 3);

        queue.push(2);
        print_queue("push 2:", queue);
        REQUIRE(queue.top() == 2);

        queue.pop();
        print_queue("pop:", queue);
        REQUIRE(queue.top() == 3);

        queue.push(1);
        print_queue("push 1:", queue);
        REQUIRE(queue.top() == 1);
    }

    // Minimum priority queue sort
    SECTION("Test MinPriorityQueue Sort") {
        std::cout << "Test MinPriorityQueue Sort\n";

        std::vector<unsigned long> values(num);
        std::iota(values.begin(), values.end(), 0);

        std::random_device rd;
        std::mt19937 rng(rd());
        std::shuffle(values.begin(), values.end(), rng);

        //
        // Create a max priority queue of those values
        //
        PriorityQueue<unsigned long> queue;
        for (auto &v : values) {
            queue.push(v);
        }

        size_t prev = queue.top();
        queue.pop();
        while (!queue.is_empty()) {
            REQUIRE(prev < queue.top());
            prev = queue.top();
            queue.pop();
        }
    }


    // Maximum priority queue sort
    SECTION("Test MaxPriorityQueue Sort") {
        std::cout << "Test MaxPriorityQueue Sort\n";

        std::vector<unsigned long> values(num);
        std::iota(values.begin(), values.end(), 0);

        std::random_device rd;
        std::mt19937 rng(rd());
        std::shuffle(values.begin(), values.end(), rng);

        // Create a max priority queue of those values
        PriorityQueue<unsigned long, std::greater<unsigned long>> queue;
        for (auto &v : values) {
            queue.push(v);
        }

        size_t prev = queue.top();
        queue.pop();
        while (!queue.is_empty()) {
            REQUIRE(prev > queue.top());
            prev = queue.top();
            queue.pop();
        }
    }

    // ---- Test IndexPriorityQueue PriorityNode ------------------------------
    // Create a random vector of nodes
    //
    std::vector<PriorityNode> nodes(num);
    for (size_t i = 0; i < num; ++i) {
        nodes[i].v(i);
    }

    std::random_device rd;
    std::mt19937 rng(rd());
    std::shuffle(nodes.begin(), nodes.end(), rng);

    // Create a minimum priority queue of nodes
    SECTION("Test MinPriorityQueue PriorityNode") {
        std::cout << "Test MinPriorityQueue PriorityNode\n";

        auto node_less =
            [&] (const PriorityNode &a, const PriorityNode &b) -> bool {
                return (a.v() < b.v());
            };

        PriorityQueue<
            PriorityNode,
            decltype(node_less)> queue(node_less);

        for (auto n : nodes) {
            queue.push(n);
        }

        auto prev = queue.top();
        queue.pop();
        while (!queue.is_empty()) {
            REQUIRE(prev.v() < queue.top().v());
            prev = queue.top();
            queue.pop();
        }
    }

    // Create a maximum priority queue of nodes
    SECTION("Test MaxPriorityQueue PriorityNode") {
        std::cout << "Test MaxPriorityQueue PriorityNode\n";

        auto node_greater =
            [&] (const PriorityNode &a, const PriorityNode &b) -> bool {
                return (a.v() > b.v());
            };
        PriorityQueue<
            PriorityNode,
            decltype(node_greater)> queue(node_greater);

        for (auto n : nodes) {
            queue.push(n);
        }

        auto prev = queue.top();
        queue.pop();
        while (!queue.is_empty()) {
            REQUIRE(prev.v() > queue.top().v());
            prev = queue.top();
            queue.pop();
        }
    }
}

/// ---- Test IndexPriorityQueue -----------------------------------------------
///
TEST_CASE("IndexPriorityQueue") {
    const size_t num = 100000;

    // Maximum priority queue
    SECTION("Test MaxIndexPriorityQueue") {
        std::cout << "Test MaxIndexPriorityQueue\n";

        // Create a min priority queue of those values
        IndexPriorityQueue<int, std::greater<int>> queue;
        REQUIRE(queue.size() == 0);
        REQUIRE(queue.is_empty());

        size_t ix1 = queue.push(1);
        print_queue("push 1:", queue);
        REQUIRE(queue.size() > 0);
        REQUIRE(!queue.is_empty());
        REQUIRE(queue.top() == 1);
        REQUIRE(queue.contains(ix1));
        REQUIRE(queue.key(ix1) == 1);

        size_t ix2 = queue.push(2);
        print_queue("push 2:", queue);
        REQUIRE(queue.top() == 2);
        REQUIRE(queue.contains(ix2));
        REQUIRE(queue.key(ix2) == 2);

        queue.modify(ix1, 3);
        print_queue("modify 1 -> 3:", queue);
        REQUIRE(queue.top() == 3);
        REQUIRE(queue.contains(ix1));
        REQUIRE(queue.key(ix1) == 3);

        size_t ix3 = queue.push(4);
        print_queue("push 4:", queue);
        REQUIRE(queue.top() == 4);
        REQUIRE(queue.contains(ix3));
        REQUIRE(queue.key(ix3) == 4);

        size_t ix4 = queue.push(5);
        print_queue("push 5:", queue);
        REQUIRE(queue.top() == 5);
        REQUIRE(queue.contains(ix4));
        REQUIRE(queue.key(ix4) == 5);

        queue.remove(ix1);
        print_queue("remove 3:", queue);
        REQUIRE(!queue.contains(ix1));

        queue.remove(ix3);
        print_queue("remove 4:", queue);
        REQUIRE(!queue.contains(ix3));
    }

    // Minimum priority queue
    SECTION("Test MinIndexPriorityQueue") {
        std::cout << "Test MinIndexPriorityQueue\n";

        // Create a min priority queue of those values
        IndexPriorityQueue<int, std::less<int>> queue;
        REQUIRE(queue.size() == 0);
        REQUIRE(queue.is_empty());

        size_t ix1 = queue.push(5);
        print_queue("push 5:", queue);
        REQUIRE(queue.size() > 0);
        REQUIRE(!queue.is_empty());
        REQUIRE(queue.top() == 5);
        REQUIRE(queue.contains(ix1));
        REQUIRE(queue.key(ix1) == 5);

        size_t ix2 = queue.push(4);
        print_queue("push 4:", queue);
        REQUIRE(queue.top() == 4);
        REQUIRE(queue.contains(ix2));
        REQUIRE(queue.key(ix2) == 4);

        queue.modify(ix1, 3);
        print_queue("modify 5 -> 3:", queue);
        REQUIRE(queue.top() == 3);
        REQUIRE(queue.contains(ix1));
        REQUIRE(queue.key(ix1) == 3);

        size_t ix3 = queue.push(2);
        print_queue("push 2:", queue);
        REQUIRE(queue.top() == 2);
        REQUIRE(queue.contains(ix3));
        REQUIRE(queue.key(ix3) == 2);

        size_t ix4 = queue.push(1);
        print_queue("push 1:", queue);
        REQUIRE(queue.top() == 1);
        REQUIRE(queue.contains(ix4));
        REQUIRE(queue.key(ix4) == 1);

        queue.remove(ix1);
        print_queue("remove 3:", queue);
        REQUIRE(!queue.contains(ix1));

        queue.remove(ix3);
        print_queue("remove 2:", queue);
        REQUIRE(!queue.contains(ix3));
    }

    // Minimum priority queue sort
    SECTION("Test MinIndexPriorityQueue Sort") {
        std::cout << "Test MinIndexPriorityQueue Sort\n";

        std::vector<unsigned long> values(num);
        std::iota(values.begin(), values.end(), 0);

        std::random_device rd;
        std::mt19937 rng(rd());
        std::shuffle(values.begin(), values.end(), rng);

        // Create a max priority queue of those values
        IndexPriorityQueue<unsigned long> queue;
        std::vector<size_t> queue_ix;
        for (auto &v : values) {
            queue_ix.push_back(queue.push(v));
        }

        size_t prev = queue.top();
        queue.pop();
        while (!queue.is_empty()) {
            REQUIRE(prev < queue.top());
            prev = queue.top();
            queue.pop();
        }
    }

    // Maximum priority queue sort
    SECTION("Test MaxIndexPriorityQueue Sort") {
        std::cout << "Test MaxIndexPriorityQueue Sort\n";

        std::vector<unsigned long> values(num);
        std::iota(values.begin(), values.end(), 0);

        std::random_device rd;
        std::mt19937 rng(rd());
        std::shuffle(values.begin(), values.end(), rng);

        // Create a max priority queue of those values
        IndexPriorityQueue<
            unsigned long,
            std::greater<unsigned long>> queue;
        std::vector<size_t> queue_ix;
        for (auto &v : values) {
            queue_ix.push_back(queue.push(v));
        }

        size_t prev = queue.top();
        queue.pop();
        while (!queue.is_empty()) {
            REQUIRE(prev > queue.top());
            prev = queue.top();
            queue.pop();
        }
    }


    // ---- Test IndexPriorityQueue PriorityNode ------------------------------
    // Create a random vector of nodes
    //
    std::vector<PriorityNode> nodes(num);
    for (size_t i = 0; i < num; ++i) {
        nodes[i].v(i);
    }

    std::random_device rd;
    std::mt19937 rng(rd());
    std::shuffle(nodes.begin(), nodes.end(), rng);

    // Create a minimum priority queue of nodes
    SECTION("Test MinIndexPriorityQueue PriorityNode") {
        std::cout << "Test MinIndexPriorityQueue PriorityNode\n";

        auto node_less =
            [&] (const PriorityNode &a, const PriorityNode &b) -> bool {
                return (a.v() < b.v());
            };

        IndexPriorityQueue<
            PriorityNode,
            decltype(node_less)> queue(node_less);

        for (auto n : nodes) {
            queue.push(n);
        }

        auto prev = queue.top();
        queue.pop();
        while (!queue.is_empty()) {
            REQUIRE(prev.v() < queue.top().v());
            prev = queue.top();
            queue.pop();
        }
    }

    // Create a maximum priority queue of nodes
    SECTION("Test MaxIndexPriorityQueue PriorityNode") {
        std::cout << "Test MaxIndexPriorityQueue PriorityNode\n";

        auto node_greater =
            [&] (const PriorityNode &a, const PriorityNode &b) -> bool {
                return (a.v() > b.v());
            };
        IndexPriorityQueue<
            PriorityNode,
            decltype(node_greater)> queue(node_greater);

        for (auto n : nodes) {
            queue.push(n);
        }

        auto prev = queue.top();
        queue.pop();
        while (!queue.is_empty()) {
            REQUIRE(prev.v() > queue.top().v());
            prev = queue.top();
            queue.pop();
        }
    }
}
