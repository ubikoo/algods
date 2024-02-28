//
// test-priority-queue.h
//
// Copyright (c) 2020 Carlos Braga
// This program is free software; you can redistribute it and/or modify it
// under the terms of the MIT License. See accompanying LICENSE.md or
// https://opensource.org/licenses/MIT.
//

#ifndef TEST_PRIORITY_ARGQUEUE_H_
#define TEST_PRIORITY_ARGQUEUE_H_

#include <iostream>
#include "algods/sorting/priority-queue.h"

/// ---- Priority node struct --------------------------------------------------
/// @brief Test the priority queue using a node structure.
///
struct PriorityNode {
    size_t v_ = 0;
    void v(size_t val) { v_ = val; }
    const size_t v() const { return v_; }

    PriorityNode() = default;
    PriorityNode(size_t v) : v_(v) {}
    ~PriorityNode() = default;

    PriorityNode(const PriorityNode &other) : v_(other.v_) {}
    PriorityNode &operator=(const PriorityNode &other) {
        if (this == &other)
            return *this;
        v_ = other.v_;
        return *this;
    }

    PriorityNode(PriorityNode &&other) : v_(std::move(other.v_)) {}
    PriorityNode &operator=(PriorityNode &&other) {
        if (this == &other)
            return *this;
        v_ = std::move(other.v_);
        return *this;
    }
};

/// ---- IndexPriorityQueue print function -------------------------------------
template<typename Key, typename Compare>
void print_queue(const std::string &name,
                  const PriorityQueue<Key, Compare> &queue)
{
    PriorityQueue<Key, Compare> copy(queue);
    std::cout << name << std::endl;
    while (!copy.is_empty()) {
        std::cout << copy.top() << " ";
        copy.pop();
    }
    std::cout << std::endl;
}

/// ---- IndexPriorityQueue print function -------------------------------------
template<typename Key, typename Compare>
void print_queue(const std::string &name,
                  const IndexPriorityQueue<Key, Compare> &queue)
{
    IndexPriorityQueue<Key, Compare> copy(queue);
    std::cout << name << std::endl;
    while (!copy.is_empty()) {
        std::cout << copy.top() << " ";
        copy.pop();
    }
    std::cout << std::endl;
}

#endif // TEST_PRIORITY_ARGQUEUE_H_
