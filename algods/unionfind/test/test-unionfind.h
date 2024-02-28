//
// test-union-find.h
//
// Copyright (c) 2020 Carlos Braga
// This program is free software; you can redistribute it and/or modify it
// under the terms of the MIT License. See accompanying LICENSE.md or
// https://opensource.org/licenses/MIT.
//

#ifndef TEST_ALGO_GRAPHS_UNION_FIND_H_
#define TEST_ALGO_GRAPHS_UNION_FIND_H_

#include <cassert>
#include <vector>
#include <cmath>
#include <numeric>
#include <random>
#include <algorithm>
#include "algods/unionfind/unionfind.h"
#include "algods/unionfind/unionfind-templ.h"

/// ---- UnionFind percolation structure ---------------------------------------
struct Percolation {
    // Member variables.
    size_t ndim_;                   // grid linear dimension
    size_t top_;                    // top boundary root
    size_t bottom_;                 // bottom boundary root
    std::vector<uint8_t> site_;     // lattice sites
    size_t max_cluster_size_;       // max cluster size
    size_t num_open_sites_;         // num open sites
    UnionFind uf_;                  // union find solver

    // Member functions.
    size_t index(size_t i, size_t j) const { return i*ndim_ + j; }

    uint8_t &site(size_t i, size_t j) { return site_[index(i, j)]; }
    const uint8_t &site(size_t i, size_t j) const {
        return site_[index(i, j)];
    }

    bool is_open(size_t i, size_t j) const { return site(i, j) == 1; }
    bool is_valid(size_t i, size_t j) const {
        return (i >= 0 && i < ndim_ && j >= 0 && j < ndim_);
    }

    double compute();
    void open(size_t i, size_t j);
    bool percolates();

    // Percolation constructor/destructor
    Percolation(size_t ndim)
        : ndim_(ndim)
        , top_(ndim*ndim)
        , bottom_(ndim*ndim+1)
        , site_(ndim*ndim)
        , max_cluster_size_(0)
        , num_open_sites_(0)
        , uf_(ndim*ndim+2) {}
    ~Percolation() = default;

    // Disable copy constructor/assignment.
    Percolation(const Percolation &other) = delete;
    Percolation &operator=(const Percolation &other) = delete;
};

///
/// @brief Run a percolation simulation. Randomly open sites on an initially
/// closed lattice until a percolating set spanning the system from the top to
/// the bottom is generated. Return the total number of open sites.
///
double Percolation::compute()
{
    // Shuffle site indices.
    std::random_device rd("/dev/urandom");  // rng device
    std::mt19937 rng(rd());                 // rng engine seeded with rd

    std::vector<size_t> site_idx(site_.size()); // lattice site indices
    std::iota(site_idx.begin(), site_idx.end(), 0);
    std::shuffle(site_idx.begin(), site_idx.end(), rng);

    // Open sites until it percolates.
    std::fill(site_.begin(), site_.end(), 0);

    size_t ix = 0;
    while (!percolates()) {
        size_t i = site_idx[ix] / ndim_;
        size_t j = site_idx[ix] - i*ndim_;
        open(i, j);
        ix++;
    }

    // Return fraction of open sites.
    return static_cast<double>(num_open_sites_) /
           static_cast<double>(site_.size());
}

///
/// @brief Connect site i and site j.
///
void Percolation::open(size_t i, size_t j)
{
    // Nothing to do if site is already open. Otherwise, open site.
    if (is_open(i, j)) {
        return;
    }
    site(i, j) = 1;

    // If any neighbour site is valid and open, connect both.
    size_t ii, jj;

    // top neighbour
    ii = i - 1;
    jj = j;
    if (is_valid(ii, jj) && is_open(ii, jj)) {
        uf_.join(index(i, j), index(ii, jj));
    } else if (!is_valid(ii, jj)) {
        uf_.join(index(i, j), top_);
    }

    // bottom neighbour
    ii = i + 1;
    jj = j;
    if (is_valid(ii, jj) && is_open(ii, jj)) {
        uf_.join(index(i, j), index(ii, jj));
    } else if (!is_valid(ii, jj)) {
        uf_.join(index(i, j), bottom_);
    }

    // left neighbour
    ii = i;
    jj = j - 1;
    if (is_valid(ii, jj) && is_open(ii, jj)) {
        uf_.join(index(i, j), index(ii, jj));
    }

    // right neighbour
    ii = i;
    jj = j + 1;
    if (is_valid(ii, jj) && is_open(ii, jj)) {
        uf_.join(index(i, j), index(ii, jj));
    }

    // Update max cluster size.
    max_cluster_size_ = std::max(max_cluster_size_, uf_.size(index(i, j)));
    num_open_sites_++;
}

///
/// @brief Does the system percolate along the i-axis?
///
bool Percolation::percolates()
{
    return (uf_.find(top_) == uf_.find(bottom_));
}

/// ---- UnionFindTempl percolation structure ---------------------------------
struct PercolationTempl {
    //
    // Member variables.
    //
    size_t ndim_;                       // grid linear dimension
    size_t top_;                        // top boundary root
    size_t bottom_;                     // bottom boundary root
    std::vector<uint8_t> site_;         // lattice sites
    size_t max_cluster_size_;           // max cluster size
    size_t num_open_sites_;             // num open sites
    UnionFindTempl<size_t> uf_;         // union find solver

    //
    // Member functions.
    //
    size_t index(size_t i, size_t j) const { return i*ndim_ + j; }

    uint8_t &site(size_t i, size_t j) { return site_[index(i, j)]; }
    const uint8_t &site(size_t i, size_t j) const {
        return site_[index(i, j)];
    }

    bool is_open(size_t i, size_t j) const { return site(i, j) == 1; }
    bool is_valid(size_t i, size_t j) const {
        return (i >= 0 && i < ndim_ && j >= 0 && j < ndim_);
    }

    double compute();
    void open(size_t i, size_t j);
    bool percolates();

    // Constructor/destructor
    PercolationTempl(size_t ndim)
        : ndim_(ndim)
        , top_(ndim*ndim)
        , bottom_(ndim*ndim+1)
        , site_(ndim*ndim)
        , max_cluster_size_(0)
        , num_open_sites_(0) {
        // make top and bottom sets
        uf_.insert(top_);
        uf_.insert(bottom_);
    }
    ~PercolationTempl() = default;

    // Copy constructor/assignment
    PercolationTempl(const PercolationTempl &other) = delete;
    PercolationTempl &operator=(const PercolationTempl &other) = delete;
};

///
/// @brief Run a percolation simulation. Randomly open sites on an initially
/// closed lattice until a percolating set spanning the system from the top to
/// the bottom is generated. Return the total number of open sites.
///
double PercolationTempl::compute()
{
    // Shuffle site indices.
    std::random_device rd("/dev/urandom");  // rng device
    std::mt19937 rng(rd());                 // rng engine seeded with rd

    std::vector<size_t> site_idx(site_.size()); // lattice site indices
    std::iota(site_idx.begin(), site_idx.end(), 0);
    std::shuffle(site_idx.begin(), site_idx.end(), rng);

    // Open sites until it percolates.
    std::fill(site_.begin(), site_.end(), 0);

    size_t ix = 0;
    while (!percolates()) {
        size_t i = site_idx[ix] / ndim_;
        size_t j = site_idx[ix] - i*ndim_;
        open(i, j);
        ix++;
    }

    // Return fraction of open sites.
    return static_cast<double>(num_open_sites_) /
           static_cast<double>(site_.size());
}

///
/// @brief Connect site i and site j.
///
void PercolationTempl::open(size_t i, size_t j)
{
    // Nothing to do if site is already open. Otherwise, open site.
    if (is_open(i, j)) {
        return;
    }
    site(i, j) = 1;
    if (!uf_.contains(index(i, j))) {
        uf_.insert(index(i, j));
    }

    // If any neighbour site is valid and open, connect both.
    size_t ii, jj;

    // top neighbour
    ii = i - 1;
    jj = j;
    if (is_valid(ii, jj) && is_open(ii, jj)) {
        uf_.join(index(i, j), index(ii, jj));
    } else if (!is_valid(ii, jj)) {
        uf_.join(index(i, j), top_);
    }

    // bottom neighbour
    ii = i + 1;
    jj = j;
    if (is_valid(ii, jj) && is_open(ii, jj)) {
        uf_.join(index(i, j), index(ii, jj));
    } else if (!is_valid(ii, jj)) {
        uf_.join(index(i, j), bottom_);
    }

    // left neighbour
    ii = i;
    jj = j - 1;
    if (is_valid(ii, jj) && is_open(ii, jj)) {
        uf_.join(index(i, j), index(ii, jj));
    }

    // right neighbour
    ii = i;
    jj = j + 1;
    if (is_valid(ii, jj) && is_open(ii, jj)) {
        uf_.join(index(i, j), index(ii, jj));
    }

    // Update max cluster size.
    max_cluster_size_ = std::max(max_cluster_size_, uf_.size(index(i, j)));
    num_open_sites_++;
}

///
/// @brief Return true if the system percolates along the i-axis and return false
/// otherwise.
///
bool PercolationTempl::percolates()
{
    return (uf_.find(top_) == uf_.find(bottom_));
}

#endif // TEST_ALGO_GRAPHS_UNION_FIND_H_
