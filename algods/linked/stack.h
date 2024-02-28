//
// stack.h
//
// Copyright (c) 2020 Carlos Braga
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the MIT Licensse.
// See accompanying LICENSE.md or https://opensource.org/licenses/MIT.
//

#ifndef LINKED_STACK_H_
#define LINKED_STACK_H_

#include "singly-node.h"

/// ---- Stack implementation ---------------------------------------------------
/// @brief Stack is an abstract data structure containing a collection of items
/// in last-in-first-out LIFO order.
///
/// The original Stack implementation in Algorithms implements the Java iterable
/// interface. The current class uses a SinglyNodeIterator as the underlying
/// iterable object.
///
/// API Stack<Item>:
///      Stack()                    create an empty stack
///      ~Stack()                   remove all items and destroy the stack
///
///      bool is_empty()             is the stack empty?
///      size_t size()               number of items in the stack
///      void clear()                remove all items from the stack
///
///      void push(Item item)       add an item to the stack
///      void pop()                 remove the last item from the stack
///
///      Item &top()                 access the top item
///      const Item &top()           const access the top item
///
/// @see Algorithms, 4th Edition by Robert Sedgewick, Kevin Wayne.
///
template<typename Item>
class Stack {
public:
    // Stack iterators.
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

    void push(const Item &item);
    void pop();

    Item &top();
    const Item &top() const;

    // Constructor/destructor.
    explicit Stack() : head_(nullptr), count_(0) {}
    ~Stack() { clear(); }

    Stack(const Stack &other);
    Stack<Item> &operator=(const Stack<Item> &other);

    Stack(Stack &&other);
    Stack<Item> &operator=(Stack<Item> &&other);

private:
    SinglyNode<Item> *head_;    // head node pointer
    size_t count_;              // node counter

    void initialize(const_iterator first, const_iterator last);
};

/// ---- Stack constructors and assign functions --------------------------------
/// @brief Copy the contents of the other stack in the same order.
///
template<typename Item>
Stack<Item>::Stack(const Stack<Item> &other) : Stack<Item>::Stack()
{
    initialize(other.begin(), other.end());
}

///
/// @brief Copy the contents of the other stack in the same order. Ensure stack
/// is empty before copying all nodes.
///
template<typename Item>
Stack<Item> &Stack<Item>::operator=(const Stack<Item> &other)
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
Stack<Item>::Stack(Stack<Item> &&other) : Stack<Item>::Stack()
{
    head_  = std::move(other.head_);
    count_ = std::move(other.count_);
}

///
/// @brief Clear and move the contents of the other without copying. Afterwards,
/// the other is an unspecified object.
///
template<typename Item>
Stack<Item> &Stack<Item>::operator=(Stack<Item> &&other)
{
    if (this == &other) {
        return *this;
    }
    clear();
    head_  = std::move(other.head_);
    count_ = std::move(other.count_);

    return *this;
}

/// ---- Stack private member functions -----------------------------------------
/// @brief Initialise the stack with the items in the range from first to last.
/// Called by the copy/assign constructor.
///
template<typename Item>
inline void Stack<Item>::initialize(const_iterator first, const_iterator last)
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

/// ---- Stack API functions ----------------------------------------------------
/// @brief Return true if the head doesn't point to a node.
///
template<typename Item>
inline bool Stack<Item>::is_empty() const
{
    return (head_ == nullptr);
}

///
/// @brief Return the number of nodes in the stack.
///
template<typename Item>
inline size_t Stack<Item>::size() const
{
    return count_;
}

///
/// @brief Delete all nodes and reset head and count.
///
template<typename Item>
inline void Stack<Item>::clear()
{
    while (!is_empty()) {
        SinglyNode<Item> *node = head_->next_;
        delete head_;
        head_ = node;
    }
    count_ = 0;
}

/// ---- Stack API mutator functions --------------------------------------------
/// @brief Create a new node at the head of the stack.
///
template<typename Item>
inline void Stack<Item>::push(const Item &item)
{
    SinglyNode<Item> *node = nullptr;
    node = new SinglyNode<Item>(item);
    node->next_ = head_;
    head_ = node;

    ++count_;
}

///
/// @brief Delete a node from the front of the stack.
///
template<typename Item>
inline void Stack<Item>::pop()
{
    assert(!is_empty() && "empty stack, out of range error");
    SinglyNode<Item> *node = head_->next_;
    delete head_;
    head_ = node;
    --count_;
}

/// ---- Stack API accessor functions -------------------------------------------
template<typename Item>
inline Item &Stack<Item>::top()
{
    assert(!is_empty() && "empty stack, out of range error");
    return head_->item_;
}

template<typename Item>
inline const Item &Stack<Item>::top() const
{
    assert(!is_empty() && "empty stack, out of range error");
    return head_->item_;
}

#endif // LINKED_STACK_H_
