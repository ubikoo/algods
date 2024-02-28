//
// graph.h
//
// Copyright (c) 2020 Carlos Braga
// This program is free software; you can redistribute it and/or modify it
// under the terms of the MIT License. See accompanying LICENSE.md or
// https://opensource.org/licenses/MIT.
//

#ifndef ALGO_GRAPHS_GRAPH_H_
#define ALGO_GRAPHS_GRAPH_H_

#include "atto/core/core.h"

/// ---- Graph implementation -------------------------------------------------
/// class Graph
/// @brief Graph implements a graph in adjacency list representation.
///
/// API Graph:
///
/// @see Algorithms, 4th Edition by Robert Sedgewick, Kevin Wayne.
///
class Graph {
public:
    Graph(size_t capacity)
        : parent_(capacity)
        , size_(capacity)
        , count_(capacity)
        , capacity_(capacity) {
        for (size_t i = 0; i < capacity; ++i) {
           parent_[i] = i;      // each {key} is initially its parent set
           size_[i] = 1;        // with a single element
        }
    }
    ~Graph() = default;

private:

};

#endif // ALGO_GRAPHS_GRAPH_H_
