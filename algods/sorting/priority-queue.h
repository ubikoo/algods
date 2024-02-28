//
// priority-queue.h
//
// Copyright (c) 2020 Carlos Braga
// This program is free software; you can redistribute it and/or modify it
// under the terms of the MIT License. See accompanying LICENSE.md or
// https://opensource.org/licenses/MIT.
//

#ifndef SORTING_PRIORITY_QUEUE_H_
#define SORTING_PRIORITY_QUEUE_H_

#include <cassert>
#include <vector>
#include <algorithm>
#include <functional>

///
/// @brief PriorityQueue represents a binary tree. A binary tree is said to
/// be heap-ordered if it maintains the invariant that for each node (a)
/// with two children (a->b) and (a->c) the following holds true:
///      a (op) b = true
///      a (op) c = true
/// where (op) is a strict weak order binary comparison operator:
///      x (op) x = false                           (irreflexivity)
///      if x (op) y = true, then y (op) x = false  (asymmetry)
///      if x (op) y = true, and
///         y (op) z = true, then
///         x (op) z = true                         (transitivity)
///
/// The heap invariant is maintained when each node compares true against
/// both of its children. The functions swim and sink restore the heap
/// invariant by swapping the node with its parent or the node with one the
/// extremum of its children.
///
///  API PriorityQueue<Key>:
///      PriorityQueue(compare) create an empty heap with comparator
///      ~PriorityQueue()       remove all items and destroy the heap
///
///      bool is_empty()          is the heap empty?
///      size_t size()           number of keys in the heap
///      void clear()            remove all keys from the heap
///
///      void push(Key)          insert a key into the heap and return its index
///      void pop()              delete the topmost key item from the heap
///      const Key &top()        const access the topmost key item in the heap
///
///      void resize(capacity)   resize the heap underlying binary tree
///      void swap(i, j)         swap keys associated with index i and index j
///      void swim(k)            move the k-key up to restore the heap order
///      void sink(k)            move the k-key down to restore the heap order
///
template<typename Key, typename Compare = std::less<Key>>
class PriorityQueue {
public:
    // Public member functions.
    bool is_empty() const;
    size_t size() const;

    void push(const Key &key);
    void pop();

    const Key &top() const;

    // Constructor/destructor.
    explicit PriorityQueue(const Compare &compare = Compare())
        : compare_(compare)
        , pq_(1)
        , capacity_(1)
        , count_(0) {}
    ~PriorityQueue() = default;

    PriorityQueue(const PriorityQueue &other)
        : compare_(other.compare_)
        , pq_(other.pq_)
        , capacity_(other.capacity_)
        , count_(other.count_) {}
    PriorityQueue &operator=(const PriorityQueue &other) {
        if (this == &other)
            return *this;
        compare_ = other.compare_;
        pq_ = other.pq_;
        capacity_ = other.capacity_;
        count_ = other.count_;
        return *this;
    }

    PriorityQueue(PriorityQueue &&other)
        : compare_(std::move(other.compare_))
        , pq_(std::move(other.pq_))
        , capacity_(std::move(other.capacity_))
        , count_(std::move(other.count_)) {}
    PriorityQueue &operator=(PriorityQueue &&other) {
        if (this == &other)
            return *this;
        compare_ = std::move(other.compare_);
        pq_ = std::move(other.pq_);
        capacity_ = std::move(other.capacity_);
        count_ = std::move(other.count_);
        return *this;
    }
private:
    // Helper functions to restore the heap invariant.
    void resize(size_t capacity);      // grow the heap to double size
    void swap(size_t i, size_t j);     // swap i-key with j-key
    void swim(size_t k);               // bottom-up reheapify
    void sink(size_t k);               // top-down reheapify

    // PriorityQueue private member variables.
    Compare compare_;               // key comparator
    std::vector<Key> pq_;           // heap-ordered complete binary tree
    size_t capacity_;               // min capacity is 1 with unused 0th key
    size_t count_;                  // number of keys in keys[1...count]
};

/// ---- PriorityQueue API functions -------------------------------------------
/// @brief Return true if the priority queue has no keys.
///
template<typename Key, typename Compare>
inline bool PriorityQueue<Key,Compare>::is_empty() const
{
    return (count_ == 0);
}

///
/// @brief Return the number of keys in the priority queue.
///
template<typename Key, typename Compare>
inline size_t PriorityQueue<Key,Compare>::size() const
{
    return count_;
}

///
/// @brief Add a new key to the priority queue.
///
template<typename Key, typename Compare>
inline void PriorityQueue<Key,Compare>::push(const Key &key)
{
    // Grow the priority queue if we reached its capacity.
    if (count_ == capacity_ - 1) {
        resize(2*capacity_);
    }

    // Add new key to the priority queue.
    count_++;
    pq_[count_] = key;
    swim(count_);
}

///
/// @brief Remove the topmost key from the priority queue.
///
template<typename Key, typename Compare>
inline void PriorityQueue<Key,Compare>::pop()
{
    // Remove the top key from the priority queue
    swap(1, count_--);      // swap top-key with bottom-key
    sink(1);                // sink bottom-key from the top

    // Shrink the priority queue if its a quarter of its capacity
    if ((count_ > 0) &&(count_ == (capacity_ - 1)/4)) {
        resize(capacity_/2);
    }
}

///
/// @brief Return a const reference to the topmost key in the priority queue.
///
template<typename Key, typename Compare>
inline const Key &PriorityQueue<Key,Compare>::top() const
{
    assert(count_ > 0 && "priority queue underflow");
    return pq_[1];
}

/// ---- PriorityQueue API helper functions ------------------------------------
/// @brief Helper function to double the size of the heap array. The heap is not
/// a dynamic array, it can only grow.
///
template<typename Key, typename Compare>
inline void PriorityQueue<Key,Compare>::resize(size_t capacity)
{
    assert(capacity > count_ && "priority queue overflow");

    // Create a copy of the heap arrays with the new capacity.
    std::vector<Key> copy(capacity);

    // Copy the items of the current heap arrays.
    capacity_ = capacity;
    for (size_t i = 1; i <= count_; ++i) {
        copy[i] = pq_[i];
    }
    pq_ = std::move(copy);
}

///
/// @brief Exchange the order of keys in the heap by swapping their indices.
///
template<typename Key, typename Compare>
inline void PriorityQueue<Key,Compare>::swap(size_t i, size_t j)
{
    std::swap(pq_[i], pq_[j]);
}

///
/// @brief Exchange node with its parent to restore the heap invariant.
///
template<typename Key, typename Compare>
inline void PriorityQueue<Key,Compare>::swim(size_t k)
{
    while (k > 1 && compare_(pq_[k], pq_[k/2])) {
        swap(k, k/2);
        k /= 2;
    }
}

///
/// @brief Exchange node with one of its children to restore the heap invariant.
///
template<typename Key, typename Compare>
inline void PriorityQueue<Key,Compare>::sink(size_t k)
{
    while (2*k <= count_) {
        size_t j = 2*k; // left child of parent k
        if (j < count_ && compare_(pq_[j+1], pq_[j])) {
            j++;        // use right child if it compares true against left
        }

        if (compare_(pq_[k], pq_[j])) {
            break;      // parent compares true against its child
        }
        swap(k, j);     // swap parent with its child
        k = j;          // go down the childs branch
    }
}

/// ---- IndexPriorityQueue implementation -------------------------------------
/// @brief IndexPriorityQueue represents an index based binary tree that supports
/// priority queue operations. It maintains a unique index associated with each
/// one of its keys.
///
/// The key index is returned at insertion point. It can then be used to refer to
/// the key and operate on it. When mutating the key value, the heap acts
/// accordingly and reorders itself to maintain its invariant.
///
/// A binary tree is said to be heap-ordered if it maintains the invariant that
/// for each node (a) with two children (a->b) and (a->c) the following holds
/// true:
///
///      a (op) b = true
///      a (op) c = true
///
/// where (op) is a strict weak order binary comparison operator:
///
///      x (op) x = false                           (irreflexivity)
///      if x (op) y = true, then y (op) x = false  (asymmetry)
///      if x (op) y = true, and
///         y (op) z = true, then
///         x (op) z = true                         (transitivity)
///
/// The heap invariant is maintained when each node compares true against both of
/// its children. The functions swim and sink restore the heap invariant by
/// swapping the node with its parent or the node with one the extremum of its
/// children.
///
///  API IndexPriorityQueue<Key>:
///      IndexPriorityQueue(compare)    create an empty heap with a comparator
///      ~IndexPriorityQueue()          remove all items and destroy the heap
///
///      bool is_empty()          is the heap empty?
///      size_t size()           number of keys in the heap
///      void clear()            remove all keys from the heap
///
///      size_t push(Key)        insert a key into the heap and return its index
///      void pop()              delete the topmost key item from the heap
///      const Key &top()        const access the topmost key item in the heap
///      const Key &key(k)       const access the k-index key item in the heap
///      bool contains(k)        is index k associated with a key?
///      void remove(k)          remove key associated with index k
///      void modify(k, Key)     modify key associated with index k
///
///      void resize(capacity)   resize the heap underlying binary tree
///      bool is_valid(k)        validate given index
///      void swap(i, j)         swap keys associated with index i and index j
///      void swim(k)            move the k-key up to restore the heap order
///      void sink(k)            move the k-key down to restore the heap order
///
template<typename Key, typename Compare = std::less<Key>>
class IndexPriorityQueue {
public:
    // Public member functions.
    bool is_empty() const;
    size_t size() const;

    size_t push(const Key &key);
    void pop();

    const Key &top() const;
    const Key &key(size_t k) const;

    bool contains(size_t k) const;
    void remove(size_t k);
    void modify(size_t k, const Key &key);

    // Constructor/destructor.
    explicit IndexPriorityQueue(const Compare &compare = Compare())
        : compare_(compare)
        , pq_(1)
        , qp_(1)
        , keys_(1)
        , capacity_(1)
        , count_(0) {}
    ~IndexPriorityQueue() = default;

    IndexPriorityQueue(const IndexPriorityQueue &other)
        : compare_(other.compare_)
        , pq_(other.pq_)
        , qp_(other.qp_)
        , keys_(other.keys_)
        , capacity_(other.capacity_)
        , count_(other.count_) {}
    IndexPriorityQueue &operator=(const IndexPriorityQueue &other) {
        if (this == &other)
            return *this;
        compare_ = other.compare_;
        pq_ = other.pq_;
        qp_ = other.qp_;
        keys_ = other.keys_;
        capacity_ = other.capacity_;
        count_ = other.count_;
        return *this;
    }

    IndexPriorityQueue(IndexPriorityQueue &&other)
        : compare_(std::move(other.compare_))
        , pq_(std::move(other.pq_))
        , qp_(std::move(other.qp_))
        , keys_(std::move(other.keys_))
        , capacity_(std::move(other.capacity_))
        , count_(std::move(other.count_)) {}
    IndexPriorityQueue &operator=(IndexPriorityQueue &&other) {
        if (this == &other)
            return *this;
        compare_ = std::move(other.compare_);
        pq_ = std::move(other.pq_);
        qp_ = std::move(other.qp_);
        keys_ = std::move(other.keys_);
        capacity_ = std::move(other.capacity_);
        count_ = std::move(other.count_);
        return *this;
    }

private:
    // Helper functions to restore the heap invariant.
    void resize(size_t capacity);      // grow the heap to double size
    bool is_valid(size_t k) const;     // validate key index
    void swap(size_t i, size_t j);     // swap i-key with j-key
    void swim(size_t k);               // bottom-up reheapify
    void sink(size_t k);               // top-down reheapify

    // IndexPriorityQueue private member variables.
    Compare compare_;           // key comparator
    std::vector<size_t> pq_;    // heap-ordered complete binary tree
    std::vector<size_t> qp_;    // inverse of pq - qp[pq[i]] = pq[qp[i]] = i
    std::vector<Key> keys_;     // array of keys with priorities
    size_t capacity_;           // min capacity is 1 with unused 0th key
    size_t count_;              // number of keys in keys[1...count]
};


/// ---- IndexPriorityQueue API functions --------------------------------------
/// @brief Return true if the priority queue has no keys.
///
template<typename Key, typename Compare>
inline bool IndexPriorityQueue<Key,Compare>::is_empty() const
{
    return (count_ == 0);
}

///
/// @brief Return the number of keys in the priority queue.
///
template<typename Key, typename Compare>
inline size_t IndexPriorityQueue<Key,Compare>::size() const
{
    return count_;
}

///
/// @brief Add a new key to the priority queue.
///
template<typename Key, typename Compare>
inline size_t IndexPriorityQueue<Key,Compare>::push(const Key &key)
{
    // Grow the priority queue if we reached its capacity.
    if (count_ == capacity_ - 1) {
        resize(2*capacity_);
    }

    // Add new key to the priority queue.
    count_++;
    qp_[pq_[count_]] = count_;
    keys_[count_] = key;
    swim(count_);
    return count_;
}

///
/// @brief Remove the topmost key from the priority queue.
///
template<typename Key, typename Compare>
inline void IndexPriorityQueue<Key,Compare>::pop()
{
    size_t top_ix = pq_[1];
    swap(1, count_--);      // swap top-key with bottom-key
    sink(1);                // sink bottom-key from the top
    qp_[top_ix] = -1;       // unset top-key at the bottom
}

///
/// @brief Return a const reference to the topmost key in the priority queue.
///
template<typename Key, typename Compare>
inline const Key &IndexPriorityQueue<Key,Compare>::top() const
{
    assert(count_ > 0 && "priority queue underflow");
    return keys_[pq_[1]];
}

///
/// @brief Return a const reference to the key associated with index k.
///
template<typename Key, typename Compare>
inline const Key &IndexPriorityQueue<Key,Compare>::key(size_t k) const
{
    assert(is_valid(k) && "invalid key index");
    assert(contains(k) && "non existent key");
    return keys_[k];
}

///
/// @brief Return true if the queue contains a key associated wint index k.
///
template<typename Key, typename Compare>
inline bool IndexPriorityQueue<Key,Compare>::contains(size_t k) const
{
    if (is_valid(k)) {
        return (qp_[k] != static_cast<size_t>(-1));
    }
    return false;
}

///
/// @brief Remove key associated with index k.
///
template<typename Key, typename Compare>
inline void IndexPriorityQueue<Key,Compare>::remove(size_t k)
{
    assert(is_valid(k) && "invalid key index");
    assert(contains(k) && "non existent key");

    size_t k_ix = qp_[k];
    swap(k_ix, count_--);   // swap k-index key with bottom
    swim(k_ix);             // restore heap order upwards
    sink(k_ix);             // restore heap order downwards
    qp_[k] = -1;            // unset key at the bottom
}

///
/// @brief Modify key associated with index k.
///
template<typename Key, typename Compare>
inline void IndexPriorityQueue<Key,Compare>::modify(size_t k, const Key &key)
{
    assert(is_valid(k) && "invalid key index");
    assert(contains(k) && "non existent key");

    keys_[k] = key;         // store the key
    swim(qp_[k]);           // restore heap order upwards
    sink(qp_[k]);           // restore heap order downwards
}


/// ---- IndexPriorityQueue API helper functions -------------------------------
/// @brief Helper function to double the size of the heap array.
/// The heap is not a dynamic array, it can only grow.
///
template<typename Key, typename Compare>
inline void IndexPriorityQueue<Key,Compare>::resize(size_t capacity)
{
    // Create a copy of the heap arrays with the new capacity.
    std::vector<size_t> copy_pq(capacity);
    std::vector<size_t> copy_qp(capacity);
    std::vector<Key> copy_keys(capacity);

    // Copy the items of the current heap arrays.
    capacity_ = capacity;
    for (size_t i = 1; i <= count_; ++i) {
        copy_pq[i] = pq_[i];
        copy_qp[i] = qp_[i];
        copy_keys[i] = keys_[i];
    }

    for (size_t i = count_+1; i < capacity_; ++i) {
        copy_pq[i] = i;
        copy_qp[i] = -1;
    }

    pq_ = std::move(copy_pq);
    qp_ = std::move(copy_qp);
    keys_ = std::move(copy_keys);
}

///
/// @brief Validate key index.
///
template<typename Key, typename Compare>
inline bool IndexPriorityQueue<Key,Compare>::is_valid(size_t k) const
{
    return (k >= 1 && k <=count_);
}

///
/// @brief Exchange the order of keys in the heap by swapping their indices.
///
template<typename Key, typename Compare>
inline void IndexPriorityQueue<Key,Compare>::swap(size_t i, size_t j)
{
    std::swap(pq_[i], pq_[j]);
    qp_[pq_[i]] = i;
    qp_[pq_[j]] = j;
}

///
/// @brief Exchange node with its parent to restore the heap invariant.
///
template<typename Key, typename Compare>
inline void IndexPriorityQueue<Key,Compare>::swim(size_t k)
{
    while (k > 1 && compare_(keys_[pq_[k]], keys_[pq_[k/2]])) {
        swap(k, k/2);
        k /= 2;
    }
}

///
/// @brief Exchange node with one of its children to restore the heap invariant.
///
template<typename Key, typename Compare>
inline void IndexPriorityQueue<Key,Compare>::sink(size_t k)
{
    while (2*k <= count_) {
        size_t j = 2*k;     // left child of parent k
        if (j < count_ && compare_(keys_[pq_[j+1]], keys_[pq_[j]])) {
            j++;            // use right child if it compares true against left
        }

        if (compare_(keys_[pq_[k]], keys_[pq_[j]])) {
            break;          // parent compares true against its child
        }
        swap(k, j);         // swap parent with its child
        k = j;              // go down the childs branch
    }
}

#endif // SORTING_PRIORITY_QUEUE_H_
