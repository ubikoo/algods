//
// singly-node.h
//
// Copyright (c) 2020 Carlos Braga
// This program is free software; you can redistribute it and/or modify it
// under the terms of the MIT License. See accompanying LICENSE.md or
// https://opensource.org/licenses/MIT.
//

#ifndef LINKED_SINGLY_NODE_H_
#define LINKED_SINGLY_NODE_H_

/// -----------------------------------------------------------------------------
/// @brief SinglyNode is a recursive data type with a data item and a reference
/// to another node of the same type.
/// It is an abstraction of an associative data structure used to by more complex
/// data structures - bags, stacks and queues.
///
/// @note Every friend of SinglyNode class has access to its members, functions
/// and constructors. These are all marked as private thereby blocking everything
/// except these.
///
/// @see Algorithms 4th Edition by Robert Sedgewick, Kevin Wayne.
/// https://en.wikipedia.org/wiki/Linked_list
///
template<typename Item> class Bag;
template<typename Item> class Stack;
template<typename Item> class Queue;
template<typename Item> class SinglyNodeIterator;
template<typename Item> class SinglyNodeConstIterator;

template<typename Item>
class  SinglyNode {
public:
    // Grant access to Bag, Stack, Queue and SinglyNodeIterator classes
    friend class Bag<Item>;
    friend class Stack<Item>;
    friend class Queue<Item>;
    friend class SinglyNodeIterator<Item>;
    friend class SinglyNodeConstIterator<Item>;

private:
    SinglyNode<Item> *next_;
    Item item_;

    explicit SinglyNode()
        : next_(nullptr) {}
    explicit SinglyNode(const Item &item)
        : next_(nullptr)
        , item_(item) {}
    ~SinglyNode() {}

    SinglyNode(const SinglyNode &other)
        : next_(other.next_)
        , item_(other.item_) {}
    SinglyNode &operator=(const SinglyNode &other) {
        if (this == &other)
            return *this;
        next_ = other.next_;
        item_ = other.item_;
        return *this;
    }

    SinglyNode(SinglyNode &&other)
        : next_(std::move(other.next_))
        , item_(std::move(other.item_)) {}
    SinglyNode &operator=(SinglyNode &&other) {
        if (this == &other)
            return *this;
        next_ = std::move(other.next_);
        item_ = std::move(other.item_);
        return *this;
    }
};

/// -----------------------------------------------------------------------------
/// @brief SinglyNodeIterator is a forward iterator, accessing a list of nodes in
/// the direction from head to tail. It needs to implement the following at least:
///
///  X a;                // default-constructor and destructor.
///  X b(a);             // copy-constructor.
///  X b = a;            // copy-assign.
///
///  *a, a->m;           // dereference operators(const).
///  ++a, a++, *a++;     // increment operators(dereferenceable).
///  a == b, a != b;     // comparison operators(const).
///  swap(a,b);          // swap operator.
///
/// The implementation needs to define the following types describing
/// the iterator properties:
///
///  typedef typename Iterator::value_type            value_type;
///  typedef typename Iterator::pointer               pointer;
///  typedef typename Iterator::reference             reference;
///  typedef typename Iterator::difference_type       difference_type;
///  typedef typename Iterator::iterator_category     iterator_category;
///
/// Standard library containers provide two types of iterators,
///  iterator type - pointing to mutable data, and
///  const_iterator type - pointing to immutable data.
///
/// It is possible to support both by providing a conversion operator:
///  operator Iterator<const Item>() const {
///      return Iterator<const Item>(itr);
///  }
///
/// Or by defining explicitly a SinglyNodeConstIterator. The later is used here.
///
/// @see Algorithms, 4th Edition by Robert Sedgewick, Kevin Wayne.
/// http://www.cplusplus.com/reference/iterator
/// https://en.cppreference.com/w/cpp/iterator
/// https://www.justsoftwaresolutions.co.uk/cplusplus/generating_sequences.html
/// https://codereview.stackexchange.com/questions/74609
///
template<typename Item>
class SinglyNodeIterator {
public:
    // Iterator member types.
    typedef SinglyNodeIterator<Item>    self;

    typedef Item                        value_type;
    typedef Item*                       pointer;
    typedef Item&                       reference;
    typedef ptrdiff_t                   difference_type;
    typedef std::forward_iterator_tag   iterator_category;

    // Dereference operators(const)
    reference operator*() const {
        assert(node_ != nullptr && "invalid iterator dereference");
        return node_->item_;
    }

    pointer operator->() const {
        assert(node_ != nullptr && "invalid iterator dereference");
        return &(operator*());
    }

    // Increment operators
    self operator++() {
        assert(node_ != nullptr && "invalid iterator increment");
        node_ = node_->next_;
        return *this;
    }

    self operator++(int) {
        assert(node_ != nullptr && "invalid iterator increment");
        SinglyNodeIterator<Item> tmp(*this);
        node_ = node_->next_;
        return tmp;
    }

    // Comparison operators(const)
    bool operator==(const SinglyNodeIterator<Item> &other) const {
        return (node_ == other.node_);
    }

    bool operator!=(const SinglyNodeIterator<Item> &other) const {
        return !operator==(other);
    }

    // Swap operator
    void swap(SinglyNodeIterator &other) {
        std::swap(node_, other.node_);
    }

    // Constructor/destructor.
    SinglyNodeIterator() : node_() {}
    explicit SinglyNodeIterator(SinglyNode<Item> *node) : node_(node) {}
    ~SinglyNodeIterator() = default;

private:
    SinglyNode<Item> *node_;
};

/// -----------------------------------------------------------------------------
/// @brief SinglyNodeConstIterator is a constant forward iterator. It is the
/// constant counterpart of SinglyNodeIterator.
///
template<typename Item>
class SinglyNodeConstIterator {
public:
    // Iterator member types.
    typedef SinglyNodeConstIterator<Item>   self;

    typedef Item                            value_type;
    typedef const Item*                     pointer;
    typedef const Item&                     reference;
    typedef ptrdiff_t                       difference_type;
    typedef std::forward_iterator_tag       iterator_category;

    // Constructor/destructor.
    SinglyNodeConstIterator() : node_() {}
    explicit SinglyNodeConstIterator(SinglyNode<Item> *node) : node_(node) {}
    ~SinglyNodeConstIterator() = default;

    // Dereference operators(const)
    reference operator*() const {
        assert(node_ != nullptr && "invalid const_iterator dereference");
        return node_->item_;
    }

    pointer operator->() const {
        assert(node_ != nullptr && "invalid const_iterator dereference");
        return &(operator*());
    }

    // Increment operators
    self operator++() {
        assert(node_ != nullptr && "invalid const_iterator increment");
        node_ = node_->next_;
        return *this;
    }

    self operator++(int) {
        assert(node_ != nullptr && "invalid const_iterator increment");
        SinglyNodeConstIterator<Item> tmp(*this);
        node_ = node_->next_;
        return tmp;
    }

    // Comparison operators(const)
    bool operator==(const SinglyNodeConstIterator<Item> &other) const {
        return (node_ == other.node_);
    }

    bool operator!=(const SinglyNodeConstIterator<Item> &other) const {
        return !operator==(other);
    }

    // Swap operator
    void swap(SinglyNodeConstIterator &other) {
        std::swap(node_, other.node_);
    }

private:
    const SinglyNode<Item> *node_;
};

#endif // LINKED_SINGLY_NODE_H_
