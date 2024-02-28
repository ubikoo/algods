//
// union-find.h
//
// Copyright (c) 2020 Carlos Braga
// This program is free software; you can redistribute it and/or modify it
// under the terms of the MIT License. See accompanying LICENSE.md or
// https://opensource.org/licenses/MIT.
//

#ifndef UNIONFIND_H_
#define UNIONFIND_H_

#include <cassert>
#include <map>
#include <vector>

///
/// @brief UnionFind implements a disjoint set data structure. It supports union
/// and find operations on the sets, together with a count operation that returns
/// the total number of disjoint sets.
/// This implementation uses weighted-quick-union by size with path compression
/// during the find operation.
///
/// UnionFind solves a dynamic connectivity problem. The input is taken as a
/// sequence of pairs of integers, where each integer represents an object of
/// some type.
/// The relation p "is connected to" q is an equivalence relation meaning:
///     symmetric:  if p is connected to q, then q is connected to p.
///     transitive: if p is connected to q and q is connected to r,
///                 then p is connected to r.
///     reflexive:  p is connected to p.
///
/// The underlying data structure is a site indexed array which maintains the
/// invariant that p and q are connected if and only if the root parent of p is
/// equal to the root parent of q.
/// The find operation returns the canonical element of the set containing a
/// given element.
/// The join operation merges the set containing element p with the set
/// containing element q. Which parent becomes child is determined by the
/// rank (or size) of the sets.
///
/// API UnionFind:
///      UnionFind()             create an empty disjoint-set data structure.
///      UnionFind(n);           create a disjoint-set with n elements.
///      ~UnionFind()            destroy all sets in the data structure.
///
///      size_t capacity()       return the number of keys.
///      size_t count()          return the number of components.
///      size_t size(p)          return the size of the tree containing p
///      bool contains(p)        return true if p exists in the ensemble
///
///      size_t find(p);         find the canonical component of element p
///      void join(p, q)         merge set containing p with set containing q
///
///      std::map<size_t, std::vector<size_t>> sets()
///          return a key-value map of all the components, where key is the
///          component parent index and value is a vector of all the keys
///          contained in the component.
///
/// @see Algorithms, 4th Edition by Robert Sedgewick, Kevin Wayne.
///
class UnionFind {
public:
    size_t capacity() const { return capacity_; }
    size_t count() const { return count_; }
    size_t size(size_t p) { return size_[find(p)]; }
    bool contains(size_t p) const;

    size_t find(size_t p);
    void join(size_t p, size_t q);

    std::map<size_t, std::vector<size_t>> sets();

    // Constructor/destructor.
    UnionFind(size_t capacity);
    ~UnionFind() = default;

    // Copy constructor/assignment.
    UnionFind(const UnionFind &other);
    UnionFind &operator=(const UnionFind &other);

    // Move constructor/assignment.
    UnionFind(UnionFind &&other);
    UnionFind &operator=(UnionFind &&other);

private:
    std::vector<size_t> parent_;        // parent index of each set
    std::vector<size_t> size_;          // component size of each set
    size_t count_;                      // number of disjoint sets
    size_t capacity_;                   // total number of keys
};

///
/// @brief Create a UnionFind disjoint set with a specified capacity.
///
inline UnionFind::UnionFind(size_t capacity)
    : parent_(capacity)
    , size_(capacity)
    , count_(capacity)
    , capacity_(capacity)
{
    for (size_t i = 0; i < capacity; ++i) {
        parent_[i] = i;      // each {key} is its own parent set
        size_[i] = 1;        // and the set contains a single element
    }
}

///
/// @brief Copy constructor/assignment.
///
inline UnionFind::UnionFind(const UnionFind &other)
    : parent_(other.parent_)
    , size_(other.size_)
    , count_(other.count_)
    , capacity_(other.capacity_)
{}

inline UnionFind &UnionFind::operator=(const UnionFind &other)
{
    if (this == &other)
        return *this;
    parent_ = other.parent_;
    size_ = other.size_;
    count_ = other.count_;
    capacity_ = other.capacity_;
    return *this;
}

///
/// @brief Move constructor/assignment.
///
inline UnionFind::UnionFind(UnionFind &&other)
    : parent_(std::move(other.parent_))
    , size_(std::move(other.size_))
    , count_(std::move(other.count_))
    , capacity_(std::move(other.capacity_))
{}

inline UnionFind &UnionFind::operator=(UnionFind &&other)
{
    if (this == &other)
        return *this;
    parent_ = std::move(other.parent_);
    size_ = std::move(other.size_);
    count_ = std::move(other.count_);
    capacity_ = std::move(other.capacity_);
    return *this;
}

///
/// @brief Return true if there is a set in the ensemble containing the key and
/// return false otherwise.
///
inline bool UnionFind::contains(size_t p) const
{
    return (p >= 0 && p < capacity_);
}

///
/// @brief Find the component to which p-key belongs to. The find operation uses
/// path compression - every element between the the query value p and the root
/// is set to to point to the root.
/// Path compression can be implemented using two passes:
///  - first pass finds the index of the root element
///  - second pass compresses the path towards the root id.
///
inline size_t UnionFind::find(size_t p)
{
    assert(contains(p) && "non existent key");

    // Find the root element of p.
    size_t root = p;
    while (root != parent_[root]) {
        root = parent_[root];
    }

    // Compress the path from p towards its root.
    size_t next = p;
    while (next != parent_[next]) {
        size_t curr = next;
        next = parent_[next];
        parent_[curr] = root;
    }

    return root;
}

///
/// @brief Join the component of p-value and the component of q-value. If these
/// are already connected, there is nothing to do. Otherwise, merge the canonical
/// components of p and q.
///
inline void UnionFind::join(size_t p, size_t q)
{
    // Get the components of p and q.
    size_t root_p = find(p);
    size_t root_q = find(q);

    // Both p and q belong to the same component - nothing to do.
    if (root_p == root_q) {
        return;
    }

    //
    // Join the smaller set into the larger set to minimize tree depth after the
    // merge operation.
    // If the p-set is smaller than q-set, merge p-set into q-set and increment
    // q-size with p-size.
    // Otherwise, merge q-set into p-set and increment p-size.
    //
    if (size_[root_p] < size_[root_q]) {
        parent_[root_p] = root_q;
        size_[root_q] += size_[root_p];
    } else {
        parent_[root_q] = root_p;
        size_[root_p] += size_[root_q];
    }

    count_--;
}

///
/// @brief Return a key-value map of the components in the set. Keys represent
/// the parent set index of each component. Values are vectors containing the
/// ids in each component.
///
inline std::map<size_t, std::vector<size_t>> UnionFind::sets()
{
    std::map<size_t, std::vector<size_t>> components;
    for (size_t p = 0; p < capacity_; ++p) {
        components[find(p)].push_back(p);
    }
    return components;
}

#endif // UNIONFIND_H_
