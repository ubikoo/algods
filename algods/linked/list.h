//
// list-templ.h
//
// Copyright (c) 2020 Carlos Braga
// This program is free software; you can redistribute it and/or modify it
// under the terms of the MIT License. See accompanying LICENSE.md or
// https://opensource.org/licenses/MIT.
//

#ifndef LINKED_LIST_H_
#define LINKED_LIST_H_

#include "doubly-node.h"

/// ---- List implementation ----------------------------------------------------
/// class List<Item>
/// @brief List represents a finite list of N(>=0) elements stored in
/// double linked order over which it can iterate.
///
/// Each doubly node contains three fields - a link to the previous
/// node, a link to the next and a data item.
/// The class supports front (and back) insertion (and removal) of
/// items with constant O(1) complexity, and forward(and reverse)
/// iteration.
///
/// The list has two pointers - head and tail - signalling the first
/// and last nodes in the list respectivelly.
/// If the list has only a single node, then the node is the first
/// and also the last - head and tail point to the same node.
/// If the list is empty, head and tail pointers are nullptr.
///
///  API List<Item>:
///      List()                 create an empty list
///      ~List()                remove all items and destroy the list
///
///      bool is_empty()         is the list empty?
///      size_t size()           number of items in the list
///      void clear()            remove all items from the list
///
///      void push_front(item)   insert an item at the front of the list
///      void push_back(item)    insert an item at the back of the list
///      void pop_front()        remove an item from the front of the list
///      void pop_back()         remove an item from the back of the list
///
///      Item &front()           access the front item
///      const Item &front()     const access the front item
///      Item &back()            access the back item
///      const Item &back()      const access the back item
///
/// @see Algorithms 4th Edition by Robert Sedgewick, Kevin Wayne.
///
template<typename Item>
class List {
public:
    // List iterators.
    typedef DoublyNodeIterator<Item> iterator;
    typedef DoublyNodeConstIterator<Item> const_iterator;

    typedef DoublyNodeReverseIterator<Item> reverse_iterator;
    typedef DoublyNodeConstReverseIterator<Item> const_reverse_iterator;

    iterator begin() { return iterator(head_); }
    const_iterator begin() const { return const_iterator(head_); }

    iterator end() { return iterator(nullptr); }
    const_iterator end() const { return const_iterator(nullptr); }

    reverse_iterator rbegin() { return reverse_iterator(tail_); }
    const_reverse_iterator rbegin() const {
        return const_reverse_iterator(tail_);
    }

    reverse_iterator rend() { return reverse_iterator(nullptr); }
    const_reverse_iterator rend() const {
        return const_reverse_iterator(nullptr);
    }

    // Public member functions.
    bool is_empty() const;
    size_t size() const;
    void clear();

    void push_front(const Item &item);
    void push_back(const Item &item);

    void pop_front();
    void pop_back();

    Item &front();
    const Item &front() const;

    Item &back();
    const Item &back() const;

    // Constructor/destructor.
    explicit List() : head_(nullptr), tail_(nullptr), count_(0) {}
    ~List() { clear(); }

    List(const List &other);
    List<Item> &operator=(const List<Item> &other);

    List(List &&other);
    List<Item> &operator=(List<Item> &&other);

private:
    DoublyNode<Item> *head_;        // head node pointer
    DoublyNode<Item> *tail_;        // tail node pointer
    size_t count_;                  // node counter

    void link(DoublyNode<Item> *left, DoublyNode<Item> *right);
    void unlink(DoublyNode<Item> *curr);
    void insert_before(DoublyNode<Item> *position, DoublyNode<Item> *node);
    void insert_after(DoublyNode<Item> *position, DoublyNode<Item> *node);
    void erase(DoublyNode<Item> *node);
};

/// ---- List constructors and assign functions ---------------------------------
/// @brief Copy the contents of the other list in the same order.
///
template<typename Item>
List<Item>::List(const List<Item> &other)
    : List<Item>::List()
{
    for (const auto it : other) {
        push_back(it);  // insert at the back to preserve sequence order
    }
}

///
/// @brief Copy the contents of the other list in the same order. Ensure list is
/// empty before copying all nodes.
///
template<typename Item>
List<Item> &List<Item>::operator=(const List<Item> &other)
{
    if (this == &other) {
        return *this;
    }
    clear();
    for (const auto &it : other) {
        push_back(it);      // insert at the back to preserve sequence order
    }

    return *this;
}

///
/// @brief Move the contents of the other without copying. Afterwards, the other
/// is an unspecified object.
///
template<typename Item>
List<Item>::List(List<Item> &&other)
    : List<Item>::List()
{
    head_  = std::move(other.head_);
    tail_  = std::move(other.tail_);
    count_ = std::move(other.count_);
}

///
/// @brief Clear and move the contents of the other without copying. Afterwards,
/// the other is an unspecified object.
///
template<typename Item>
List<Item> &List<Item>::operator=(List<Item> &&other)
{
    if (this == &other) {
        return *this;
    }
    clear();
    head_  = std::move(other.head_);
    tail_  = std::move(other.tail_);
    count_ = std::move(other.count_);

    return *this;
}

/// ---- List private member functions ------------------------------------------
/// @brief Link the left node to the right node.
///
template<typename Item>
inline void List<Item>::link(
    DoublyNode<Item> *left,
    DoublyNode<Item> *right)
{
    left->next_ = right;
    right->prev_ = left;
}

///
/// @brief Unlink the node from the list./ Note that deleting the last node of
/// the list (with null curr->next and curr->prev pointers) sets both head and
/// tail pointers to null, thereby handling the removal from a one-element list.
///
template<typename Item>
inline void List<Item>::unlink(DoublyNode<Item> *curr)
{
    //
    // If we have a previous node, link it to the next node link.
    // Otherwise current node is the front, relink the head instead.
    //
    if (curr->prev_) {
        curr->prev_->next_ = curr->next_;
    } else {
        head_ = curr->next_;
    }

    //
    // If we have a next node, link it to the previous node link.
    // Otherwise current node is the back, relink the tail instead.
    //
    if (curr->next_) {
        curr->next_->prev_ = curr->prev_;
    } else {
        tail_ = curr->prev_;
    }
}

///
/// @brief Insert a node before a given position position.
///
template<typename Item>
inline void List<Item>::insert_before(
    DoublyNode<Item> *position,
    DoublyNode<Item> *node)
{
    if (position->prev_) {
       link(position->prev_, node);
    }
    link(node, position);
}

///
/// @brief Insert a node after a given position.
///
template<typename Item>
inline void List<Item>::insert_after(
    DoublyNode<Item> *position,
    DoublyNode<Item> *node)
{
    if (position->next_) {
       link(node, position->next_);
    }
    link(position, node);
}

///
/// @brief Unlink and delete a node from the list.
///
template<typename Item>
inline void List<Item>::erase(DoublyNode<Item> *node)
{
    unlink(node);
    delete node;
}

/// ---- List API functions -----------------------------------------------------
/// @brief Return true if the head doesn't point to a node.
///
template<typename Item>
inline bool List<Item>::is_empty() const
{
    return (head_ == nullptr);
}

///
/// @brief Return the number of nodes in the list
///
template<typename Item>
inline size_t List<Item>::size() const
{
    return count_;
}

///
/// @brief Delete all the nodes and reset the head and tail pointers.
///
template<typename Item>
inline void List<Item>::clear()
{
    while (!is_empty()) {
        pop_front();
    }
}

/// ---- List API mutator functions ---------------------------------------------
/// @brief Create a new node at the head of the list. If list is empty, update
/// head and tail to point to the new node. Otherwise, insert the node before the
/// head.
///
template<typename Item>
inline void List<Item>::push_front(const Item &item)
{
    DoublyNode<Item> *node = nullptr;
    node = new DoublyNode<Item>(item);

    if (is_empty()) {
        tail_ = node;
    } else {
        insert_before(head_, node);
    }
    head_ = node;

    ++count_;
}

///
/// @brief Create a new node at the tail of the list. If list is empty, update
/// head and tail to point to the new node. Otherwise, insert the node after the
/// head.
///
template<typename Item>
inline void List<Item>::push_back(const Item &item)
{
    DoublyNode<Item> *node = nullptr;
    node = new DoublyNode<Item>(item);

    if (is_empty()) {
        head_ = node;
    } else {
        insert_after(tail_, node);
    }
    tail_ = node;

    ++count_;
}

///
/// @brief Delete a node from the front of the list.
///
template<typename Item>
inline void List<Item>::pop_front()
{
    assert(!is_empty() && "empty list, out of range error");
    erase(head_);
    --count_;
}

///
/// @brief Delete a node from the front of the list.
///
template<typename Item>
inline void List<Item>::pop_back()
{
    assert(!is_empty() && "empty list, out of range error");
    erase(tail_);
    --count_;
}

/// ---- List API accessor functions --------------------------------------------
template<typename Item>
inline Item &List<Item>::front()
{
    assert(!is_empty() && "empty list, out of range error");
    return head_->item_;
}

template<typename Item>
inline const Item &List<Item>::front() const
{
    assert(!is_empty() && "empty list, out of range error");
    return head_->item_;
}

template<typename Item>
inline Item &List<Item>::back()
{
    assert(!is_empty() && "empty list, out of range error");
    return tail_->item_;
}

template<typename Item>
inline const Item &List<Item>::back() const
{
    assert(!is_empty() && "empty list, out of range error");
    return tail_->item_;
}

#endif // LINKED_LIST_H_
