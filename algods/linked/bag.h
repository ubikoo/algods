//
// bag-templ.h
//
// Copyright (c) 2020 Carlos Braga
// This program is free software; you can redistribute it and/or modify it
// under the terms of the MIT License. See accompanying LICENSE.md or
// https://opensource.org/licenses/MIT.
//

#ifndef LINKED_BAG_H_
#define LINKED_BAG_H_

#include "singly-node.h"

/// ---- Bag implementation -----------------------------------------------------
/// @brief Bag is an abstract data structure containing an unordered collection
/// of items. An item can be inserted but not removed.
/// A Bag data structure simply holds a collection of items without any specific
/// order over which it can iterate.
///
/// The original Bag implementation in Algorithms implements the Java iterable
/// interface. The current class uses a SinglyNodeIterator as the underlying
/// iterable object.
///
/// API Bag<Item>:
///      Bag()                      create an empty bag
///      ~Bag()                     remove all items and destroy the bag
///
///      bool is_empty()             is the bag empty?
///      size_t size()               number of items in the bag
///      void clear()                remove all items from the bag
///
///      void add(Item item)        add an item to the bag
/// @see Algorithms, 4th Edition by Robert Sedgewick, Kevin Wayne.
///
template<typename Item>
class Bag {
public:
    // Bag iterators.
    typedef SinglyNodeIterator<Item> iterator;
    typedef SinglyNodeConstIterator<Item> const_iterator;

    iterator begin() { return iterator(head_); }
    const_iterator begin() const { return const_iterator(head_); }

    iterator end() { return iterator(nullptr); }
    const_iterator end() const { return const_iterator(nullptr); }

    // Public member functions.
    bool is_empty() const;
    size_t size() const;
    void clear();

    void add(const Item &item);

    // Constructor/destructor.
    explicit Bag() : head_(nullptr), count_(0) {}
    ~Bag() { clear(); }

    Bag(const Bag &other);
    Bag<Item> &operator=(const Bag<Item> &other);

    Bag(Bag &&other);
    Bag<Item> &operator=(Bag<Item> &&other);

private:
    SinglyNode<Item> *head_;    // head node pointer
    size_t count_;              // node counter

    void initialize(const_iterator first, const_iterator last);
};


/// ---- Bag constructors and assign functions ----------------------------------
/// @brief Copy the contents of the other in the same order.
///
template<typename Item>
Bag<Item>::Bag(const Bag<Item> &other) : Bag<Item>::Bag()
{
    initialize(other.begin(), other.end());
}

///
/// @brief Copy the nodes of the other in the same order. Ensure is empty before
/// copying all the nodes.
///
template<typename Item>
Bag<Item> &Bag<Item>::operator=(const Bag<Item> &other)
{
    if (this == &other) {
        return *this;
    }
    initialize(other.begin(), other.end());
    return *this;
}

///
/// @brief Move the contents of the other without copying. Afterwards, the other
/// is an unspecified object.
///
template<typename Item>
Bag<Item>::Bag(Bag<Item> &&other) : Bag<Item>::Bag()
{
    head_  = std::move(other.head_);
    count_ = std::move(other.count_);
}

///
/// @brief Clear and move the contents of the other without copying. Afterwards,
/// the other is an unspecified object.
///
template<typename Item>
Bag<Item> &Bag<Item>::operator=(Bag<Item> &&other)
{
    if (this == &other) {
        return *this;
    }

    clear();
    head_  = std::move(other.head_);
    count_ = std::move(other.count_);

    return *this;
}

/// ---- Bag private member functions -------------------------------------------
/// @brief Initialise the bag with the items in the range from first to last.
/// Called by the copy/assign constructor.
///
template<typename Item>
inline void Bag<Item>::initialize(const_iterator first, const_iterator last)
{
    clear();
    try {
        SinglyNode<Item> *to = head_;
        for (const_iterator it = first; it != last; ++it) {
            SinglyNode<Item> *node = new SinglyNode<Item>(*it);
            if (to) {
                to->next_ = node;
            }
            to = node;
            ++count_;
        }
    } catch (std::exception& e) {
        clear();
        throw std::runtime_error(e.what());
    }
}

/// ---- Bag API functions ------------------------------------------------------
/// @brief Return true if the head doesn't point to a node.
///
template<typename Item>
inline bool Bag<Item>::is_empty() const
{
    return (head_ == nullptr);
}

///
/// @brief Return the number of nodes in the bag.
///
template<typename Item>
inline size_t Bag<Item>::size() const
{
    return count_;
}

///
/// @brief Delete all nodes and reset head and count.
///
template<typename Item>
inline void Bag<Item>::clear()
{
    while (!is_empty()) {
        SinglyNode<Item> *node = head_->next_;
        delete head_;
        head_ = node;
    }
    count_ = 0;
}

/// ---- Bag API mutator functions ----------------------------------------------
/// @brief Create a new node at the head of the bag.
///
template<typename Item>
inline void Bag<Item>::add(const Item &item)
{
    SinglyNode<Item> *node = nullptr;
    node = new SinglyNode<Item>(item);
    node->next_ = head_;
    head_ = node;
    ++count_;
}

#endif // LINKED_BAG_H_
