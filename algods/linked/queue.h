//
// queue.h
//
// Copyright (c) 2020 Carlos Braga
// This program is free software; you can redistribute it and/or modify it
// under the terms of the MIT License. See accompanying LICENSE.md or
// https://opensource.org/licenses/MIT.
//

#ifndef LINKED_QUEUE_H_
#define LINKED_QUEUE_H_

#include "singly-node.h"

/// ---- Queue implementation -------------------------------------------------
/// @brief Queue is an abstract data structure containing a collection of items
/// in last-in-first-out LIFO order.
/// The original Queue implementation in Algorithms implements the Java iterable
/// interface. The current class uses a NodeIterator as the underlying iterable
/// object.
///
/// API Queue<Item>:
///      Queue()                    create an empty queue
///      ~Queue()                   remove all items and destroy the queue
///
///      bool is_empty()             is the queue empty?
///      size_t size()               number of items in the queue
///      void clear()                remove all items from the queue
///
///      void enqueue(Item item)    add an item to the queue
///      void dequeue()             remove the last item from the queue
///
///      Item &front()               access the front item
///      const Item &front()         const access the front item
///      Item &back()                access the back item
///      const Item &back()          const access the back item
///
/// @see Algorithms, 4th Edition by Robert Sedgewick, Kevin Wayne.
///
template<typename Item>
class Queue {
public:
    // Queue iterators.
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

    void enqueue(const Item &item);
    void dequeue();

    Item &front();
    const Item &front() const;

    Item &back();
    const Item &back() const;

    // Constructor/destructor.
    explicit Queue() : head_(nullptr), tail_(nullptr), count_(0) {}
    ~Queue() { clear(); }

    Queue(const Queue &other);
    Queue<Item> &operator=(const Queue<Item> &other);

    Queue(Queue &&other);
    Queue<Item> &operator=(Queue<Item> &&other);

private:
    SinglyNode<Item> *head_;    // head node pointer
    SinglyNode<Item> *tail_;    // tail node pointer
    size_t count_;              // node counter

    void initialize(const_iterator first, const_iterator last);
};

/// ---- Queue constructors and assign functions ------------------------------
/// @brief Copy the contents of the other queue in the same order.
///
template<typename Item>
Queue<Item>::Queue(const Queue<Item> &other) : Queue<Item>::Queue()
{
    initialize(other.begin(), other.end());
}

///
/// @brief Copy the contents of the other queue in the same order. Ensure queue
/// is empty before copying all nodes.
///
template<typename Item>
Queue<Item> &Queue<Item>::operator=(const Queue<Item> &other)
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
Queue<Item>::Queue(Queue<Item> &&other) : Queue<Item>::Queue()
{
    head_  = std::move(other.head_);
    tail_  = std::move(other.head_);
    count_ = std::move(other.count_);
}

///
/// @brief Clear and move the contents of the other without copying. Afterwards,
/// the other is an unspecified object.
///
template<typename Item>
Queue<Item> &Queue<Item>::operator=(Queue<Item> &&other)
{
    if (this == &other) {
        return *this;
    }
    clear();
    head_  = std::move(other.head_);
    tail_  = std::move(other.head_);
    count_ = std::move(other.count_);

    return *this;
}

/// ---- Queue private member functions ---------------------------------------
/// @brief Initialise the queue with the items in the range from first to last.
/// Called by the copy/assign constructor.
///
template<typename Item>
inline void Queue<Item>::initialize(const_iterator first, const_iterator last)
{
    clear();
    for (const_iterator it = first; it != last; ++it) {
        enqueue(*it);
    }
}

/// ---- Queue API functions --------------------------------------------------
/// @brief Return true if the head doesn't point to a node.
///
template<typename Item>
inline bool Queue<Item>::is_empty() const
{
    return (head_ == nullptr);
}

///
/// @brief Return the number of nodes in the queue.
///
template<typename Item>
inline size_t Queue<Item>::size() const
{
    return count_;
}

///
/// @brief Delete all nodes and reset head and count.
///
template<typename Item>
inline void Queue<Item>::clear()
{
    while (!is_empty()) {
        dequeue();
    }
}

/// ---- Queue API mutator functions ------------------------------------------
/// @brief Create a new node at the head of the queue.
///
template<typename Item>
inline void Queue<Item>::enqueue(const Item &item)
{
    SinglyNode<Item> *node = nullptr;
    node = new SinglyNode<Item>(item);

    //
    // If queue is empty, update both head and tail.
    // Otherwise, link the back node to the new one.
    //
    if (is_empty()) {
        head_ = node;
    } else {
        tail_->next_ = node;
    }
    tail_ = node;

    ++count_;
}

///
/// @brief Delete a node from the front of the queue.
///
template<typename Item>
inline void Queue<Item>::dequeue()
{
    assert(!is_empty() && "empty queue, out of range error");

    SinglyNode<Item> *node = head_->next_;
    delete head_;
    head_ = node;

    if (is_empty()) {
        tail_ = nullptr;    // if queue is empty, reset the tail too.

    }
    --count_;
}

/// ---- Queue API accessor functions -----------------------------------------
template<typename Item>
inline Item &Queue<Item>::front()
{
    assert(!is_empty() && "empty queue, out of range error");
    return head_->item_;
}

template<typename Item>
inline const Item &Queue<Item>::front() const
{
    assert(!is_empty() && "empty queue, out of range error");
    return head_->item_;
}

template<typename Item>
inline Item &Queue<Item>::back()
{
    assert(!is_empty() && "empty queue, out of range error");
    return tail_->item_;
}

template<typename Item>
inline const Item &Queue<Item>::back() const
{
    assert(!is_empty() && "empty queue, out of range error");
    return tail_->item_;
}

#endif // LINKED_QUEUE_H_
