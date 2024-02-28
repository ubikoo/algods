//
// unionfind-templ.h
//
// Copyright (c) 2020 Carlos Braga
// This program is free software; you can redistribute it and/or modify it
// under the terms of the MIT License. See accompanying LICENSE.md or
// https://opensource.org/licenses/MIT.
//

#ifndef UNIONFIND_TEMPL_H_
#define UNIONFIND_TEMPL_H_

#include <cassert>
#include <map>
#include <vector>

///
/// @brief UnionFindTempl maintains a templated disjoint-set data structure.
/// It supports union and find operations on the sets, together with a count
/// operation that returns the total number of component sets.
/// This implementation maintains an index map of unique indices for each key,
/// and it uses the weighted-quick-union by size with path compression during
/// the find operation. See UnionFind for more details.
///
/// API UnionFindTempl<Key>:
///      UnionFindTempl()        create an empty disjoint-set data structure
///      ~UnionFindTempl()       destroy all sets in the data structure
///
///      size_t capacity()       return the number of keys.
///      size_t count()          return the number of components.
///      size_t size(key)        return the size of the tree containing key
///      bool contains(key)      return true if key exists in the ensemble
///
///      size_t insert(key)      create a new set in the ensemble with key
///      size_t find(key)        find the canonical component of key
///      void join(p, q)         merge set containing p with set containing q
///
///      void clear()            clear all the sets in the ensemble
///      void merge(other)       merge all the sets from another ensemble
///
///      std::map<size_t, std::vector<size_t>> sets()
///          return a key-value map of all the components, where key is the
///          component parent index and value is a vector of all the keys
///          contained in the component.
///
/// @see Algorithms, 4th Edition by Robert Sedgewick, Kevin Wayne.
///
template<typename Key>
class UnionFindTempl {
public:
    size_t capacity() const { return capacity_; }
    size_t count() const { return count_; }
    size_t size(const Key &key) { return size_[find(key)]; }
    bool contains(const Key &key) const;

    size_t insert(const Key &key);
    size_t find(const Key &key);
    void join(const Key &key_p, const Key &key_q);

    void clear();
    void merge(const UnionFindTempl &other);
    std::map<size_t, std::vector<size_t>> sets();

    // Constructor/destructor.
    UnionFindTempl() : count_(0), capacity_(0) {}
    ~UnionFindTempl() = default;

    // Copy constructor/assignment.
    UnionFindTempl(const UnionFindTempl &other);
    UnionFindTempl &operator=(const UnionFindTempl &other);

    // Move constructor/assignment.
    UnionFindTempl(UnionFindTempl &&other);
    UnionFindTempl &operator=(UnionFindTempl &&other);

private:
    std::map<Key, size_t> indexmap_;    // key-set index map
    std::vector<size_t> parent_;        // parent index of each set
    std::vector<size_t> size_;          // component size of each set
    size_t count_;                      // number of disjoint sets
    size_t capacity_;                   // total number of keys
};

///
/// @brief Copy constructor/assignment.
///
template<typename Key>
inline UnionFindTempl<Key>::UnionFindTempl(
    const UnionFindTempl<Key> &other)
    : indexmap_(other.indexmap_)
    , parent_(other.parent_)
    , size_(other.size_)
    , count_(other.count_)
    , capacity_(other.capacity_)
{}

template<typename Key>
inline UnionFindTempl<Key> &UnionFindTempl<Key>::operator=(
    const UnionFindTempl<Key> &other)
{
    if (this == &other)
        return *this;
    indexmap_ = other.indexmap_;
    parent_ = other.parent_;
    size_ = other.size_;
    count_ = other.count_;
    capacity_ = other.capacity_;
    return *this;
}

///
/// @brief Move constructor/assignment.
///
template<typename Key>
inline UnionFindTempl<Key>::UnionFindTempl(UnionFindTempl<Key> &&other)
    : indexmap_(std::move(other.indexmap_))
    , parent_(std::move(other.parent_))
    , size_(std::move(other.size_))
    , count_(std::move(other.count_))
    , capacity_(std::move(other.capacity_))
{}

template<typename Key>
inline UnionFindTempl<Key> &UnionFindTempl<Key>::operator=(
    UnionFindTempl<Key> &&other)
{
    if (this == &other)
        return *this;
    indexmap_ = std::move(other.indexmap_);
    parent_ = std::move(other.parent_);
    size_ = std::move(other.size_);
    count_ = std::move(other.count_);
    capacity_ = std::move(other.capacity_);
    return *this;
}

///
/// @brief Does the disjoint set ensmble contains the specified key?
///
template<typename Key>
inline bool UnionFindTempl<Key>::contains(const Key &key) const
{
    return (indexmap_.find(key) != indexmap_.end());
}

///
/// @brief Create a new set {key} containing the key. The key object must not
/// appear in any other set in the ensemble. The parent of the new set is the
/// key itself.
///
template<typename Key>
inline size_t UnionFindTempl<Key>::insert(const Key &key)
{
    assert(!contains(key) && "duplicate key in the ensemble");

    size_t id = indexmap_.size();   // unique index for the new set
    indexmap_[key] = id;            // create a new set {key} with index id
    parent_.push_back(id);          // each {key} is initially its parent set
    size_.push_back(1);             // with a single element
    count_++;                       // increment number of components
    capacity_++;                    // increment number of keys

    return indexmap_[key];
}

///
/// @brief Find the component to which p-key belongs to. The find operation uses
/// path compression - every element between the the query key p and the root is
/// set to to point to the root.
/// Path compression can be implemented using two passes - one to find the root
/// and the second to compress the path towards the root id.
///
template<typename Key>
inline size_t UnionFindTempl<Key>::find(const Key &key)
{
    assert(contains(key) && "non existent key");

    // Find the root component of key.
    size_t root = indexmap_[key];
    while (root != parent_[root]) {
        root = parent_[root];
    }

    // Compress the path from key towards its root.
    size_t next = indexmap_[key];
    while (next != parent_[next]) {
        size_t curr = next;
        next = parent_[next];
        parent_[curr] = root;
    }

    return root;
}

///
/// @brief Join the component of p-key and the component of q-key. If these are
/// already connected, there is nothing to do. Otherwise, merge the canonical
/// components of p and q.
///
template<typename Key>
inline void UnionFindTempl<Key>::join(const Key &key_p, const Key &key_q)
{
    // Get the components of p and q.
    size_t root_p = find(key_p);
    size_t root_q = find(key_q);

    // Both p and q belong to the same component - nothing to do.
    if (root_p == root_q) {
        return;
    }

    //
    // Join the smaller set into the larger set to minimize the tree depth after
    // the merge operation. If the p-set is smaller than q-set, merge p-set into
    // q-set and increment q-size with p-size. Otherwise, merge q-set into p-set
    // and increment p-size.
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
/// @brief Clear all the union find data structures.
///
template<typename Key>
inline void UnionFindTempl<Key>::clear()
{
    indexmap_.clear();     // clear the set of key ids
    parent_.clear();    // clear the key parent ids
    size_.clear();      // clear the key set size
    count_ = 0;         // reset the component count
}

///
/// @brief Merge the sets of the rhs ensemble into the current one.
///
template<typename Key>
inline void UnionFindTempl<Key>::merge(const UnionFindTempl &rhs)
{
    size_t offset = indexmap_.size();          // offset for each new key id

    for (const auto &it : rhs.indexmap_) {
        const Key &key = it.first;          // rhs key reference
        const size_t &id = it.second;       // rhs key id

        const size_t parent = rhs.parent_[id];
        const size_t size = rhs.size_[id];

        assert(!contains(key) && "duplicate key in the ensemble");

        indexmap_[key] = id + offset;       // set {key} with offset id
        parent_.push_back(parent + offset); // parent set with offset id
        size_.push_back(size);              // with size equal to rhs size
    }

    count_ += rhs.count_;           // increment the number of components
    capacity_ += rhs.capacity_;     // increment the number of keys
}

///
/// @brief Return a key-value map of the components in the ensemble. The keys
/// represent the parent set index of each component. The values are vectors
/// containing the keys in each component.
///
template<typename Key>
inline std::map<size_t, std::vector<size_t>> UnionFindTempl<Key>::sets()
{
    std::map<size_t, std::vector<Key>> components;
    for (auto it : indexmap_) {
        const Key &key = it.first;      // get key reference
        components[find(key)].push_back(key);
    }
    return components;
}

#endif // UNIONFIND_TEMPL_H_
