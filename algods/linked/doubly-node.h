//
// doubly-node.h
//
// Copyright (c) 2020 Carlos Braga
// This program is free software; you can redistribute it and/or modify it
// under the terms of the MIT License. See accompanying LICENSE.md or
// https://opensource.org/licenses/MIT.
//

#ifndef LINKED_DOUBLY_NODE_H_
#define LINKED_DOUBLY_NODE_H_

/// -----------------------------------------------------------------------------
/// @brief DoublyNode is a recursive data type with a data item and a reference
/// to a next and previous nodes of the same type.
/// It is an abstraction of an associative data structure to be used by doubly
/// linked lists.
///
/// @note Every friend of DoublyNode class has access to its members, functions
/// and constructors. These are all marked as private thereby blocking everything
/// except these.
///
/// @see Algorithms 4th Edition by Robert Sedgewick, Kevin Wayne.
/// https://en.wikipedia.org/wiki/Doubly_linked_list
///
template<typename Item> class List;
template<typename Item> class DoublyNodeIterator;
template<typename Item> class DoublyNodeConstIterator;
template<typename Item> class DoublyNodeReverseIterator;
template<typename Item> class DoublyNodeConstReverseIterator;

template<typename Item>
class DoublyNode {
public:
    // Grant access to List and List Iterator classes
    friend class List<Item>;
    friend class DoublyNodeIterator<Item>;
    friend class DoublyNodeConstIterator<Item>;
    friend class DoublyNodeReverseIterator<Item>;
    friend class DoublyNodeConstReverseIterator<Item>;

private:
    // DoublyNode members
    DoublyNode<Item> *next_;
    DoublyNode<Item> *prev_;
    Item item_;

    // Constructor/destructor.
    explicit DoublyNode() : next_(nullptr), prev_(nullptr) {}
    explicit DoublyNode(const Item &item)
        : next_(nullptr)
        , prev_(nullptr)
        , item_(item) {}
    ~DoublyNode() {}

    DoublyNode(const DoublyNode &other)
        : next_(other.next_)
        , prev_(other.prev_)
        , item_(other.item_) {}
    DoublyNode &operator=(const DoublyNode &other) {
        if (this == &other)
            return *this;
        next_ = other.next_;
        prev_ = other.prev_;
        item_ = other.item_;
        return *this;
    }

    DoublyNode(DoublyNode &&other)
        : next_(std::move(other.next_))
        , prev_(std::move(other.prev_))
        , item_(std::move(other.item_)) {}
    DoublyNode &operator=(DoublyNode &&other) {
        if (this == &other)
            return *this;
        next_ = std::move(other.next_);
        prev_ = std::move(other.prev_);
        item_ = std::move(other.item_);
        return *this;
    }
};

/// -----------------------------------------------------------------------------
/// @brief DoublyNodeIterator is a forward iterator, accessing the nodes in the
/// List in the direction from the head to the tail. It needs to implement the
/// following at least:
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
/// The implementation needs to define the following types describing the iterator
/// properties:
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
/// Or by defining explicitly a DoublyNodeConstIterator. The later is used here.
///
/// @see Algorithms, 4th Edition by Robert Sedgewick, Kevin Wayne.
/// http://www.cplusplus.com/reference/iterator
/// https://en.cppreference.com/w/cpp/iterator
/// https://www.justsoftwaresolutions.co.uk/cplusplus/generating_sequences.html
/// https://codereview.stackexchange.com/questions/74609
///
template<typename Item>
class DoublyNodeIterator {
public:
    // Iterator member types.
    typedef DoublyNodeIterator<Item>    self;

    typedef Item                        value_type;
    typedef Item*                       pointer;
    typedef Item&                       reference;
    typedef ptrdiff_t                   difference_type;
    typedef std::forward_iterator_tag   iterator_category;

    // Constructor/destructor.
    DoublyNodeIterator() : node_() {}
    explicit DoublyNodeIterator(DoublyNode<Item> *node) : node_(node) {}
    ~DoublyNodeIterator() = default;

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
        DoublyNodeIterator<Item> tmp(*this);
        node_ = node_->next_;
        return tmp;
    }

    // Comparison operators(const)
    bool operator==(const DoublyNodeIterator<Item> &other) const {
        return (node_ == other.node_);
    }

    bool operator!=(const DoublyNodeIterator<Item> &other) const {
        return !operator==(other);
    }

    // Swap operator
    void swap(DoublyNodeIterator &other) {
        std::swap(node_, other.node_);
    }

private:
    DoublyNode<Item> *node_;
};

///
/// @brief DoublyNodeConstIterator is a constant forward iterator. It is the
/// constant counterpart of DoublyNodeIterator.
///
template<typename Item>
class DoublyNodeConstIterator {
public:
    // Iterator member types.
    typedef DoublyNodeConstIterator<Item>   self;

    typedef Item                            value_type;
    typedef const Item*                     pointer;
    typedef const Item&                     reference;
    typedef ptrdiff_t                       difference_type;
    typedef std::forward_iterator_tag       iterator_category;

    // Constructor/destructor.
    DoublyNodeConstIterator() : node_() {}
    explicit DoublyNodeConstIterator(DoublyNode<Item> *node) : node_(node) {}
    ~DoublyNodeConstIterator() = default;

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
        DoublyNodeConstIterator<Item> tmp(*this);
        node_ = node_->next_;
        return tmp;
    }

    // Comparison operators(const)
    bool operator==(const DoublyNodeConstIterator<Item> &other) const {
        return (node_ == other.node_);
    }

    bool operator!=(const DoublyNodeConstIterator<Item> &other) const {
        return !operator==(other);
    }

    // Swap operator
    void swap(DoublyNodeConstIterator &other) {
        std::swap(node_, other.node_);
    }

private:
    const DoublyNode<Item> *node_;
};

/// -----------------------------------------------------------------------------
/// @brief DoublyNodeReverseIterator is a reverse iterator, accessing the nodes
/// in the List in the direction from the tail to the head. It needs to implement
/// the following at least:
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
/// Or by defining explicitly a DoublyNodeConstIterator. The later is used here.
///
/// @see Algorithms, 4th Edition by Robert Sedgewick, Kevin Wayne.
/// http://www.cplusplus.com/reference/iterator
/// https://en.cppreference.com/w/cpp/iterator
/// https://www.justsoftwaresolutions.co.uk/cplusplus/generating_sequences.html
/// https://codereview.stackexchange.com/questions/74609
///
template<typename Item>
class DoublyNodeReverseIterator {
public:
    // Iterator member types.
    typedef DoublyNodeReverseIterator<Item> self;

    typedef Item                            value_type;
    typedef Item*                           pointer;
    typedef Item&                           reference;
    typedef ptrdiff_t                       difference_type;
    typedef std::forward_iterator_tag       iterator_category;

    // Constructor/destructor.
    DoublyNodeReverseIterator() : node_() {}
    explicit DoublyNodeReverseIterator(DoublyNode<Item> *node) : node_(node) {}
    ~DoublyNodeReverseIterator() = default;

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
        node_ = node_->prev_;
        return *this;
    }

    self operator++(int) {
        assert(node_ != nullptr && "invalid iterator increment");
        DoublyNodeReverseIterator<Item> tmp(*this);
        node_ = node_->prev_;
        return tmp;
    }

    // Comparison operators(const)
    bool operator==(const DoublyNodeReverseIterator<Item> &other) const {
        return (node_ == other.node_);
    }

    bool operator!=(const DoublyNodeReverseIterator<Item> &other) const {
        return !operator==(other);
    }

    // Swap operator
    void swap(DoublyNodeReverseIterator &other) {
        std::swap(node_, other.node_);
    }

private:
    DoublyNode<Item> *node_;
};

/// -----------------------------------------------------------------------------
/// @brief DoublyNodeConstReverseIterator is a constant forward iterator.
/// It is the constant counterpart of DoublyNodeReverseIterator.
///
template<typename Item>
class DoublyNodeConstReverseIterator {
public:
    // Iterator member types.
    typedef DoublyNodeConstReverseIterator<Item>    self;

    typedef Item                                    value_type;
    typedef const Item*                             pointer;
    typedef const Item&                             reference;
    typedef ptrdiff_t                               difference_type;
    typedef std::forward_iterator_tag               iterator_category;

    // Constructor/destructor.
    DoublyNodeConstReverseIterator() : node_() {}
    explicit DoublyNodeConstReverseIterator(DoublyNode<Item> *node)
        : node_(node) {}
    ~DoublyNodeConstReverseIterator() = default;

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
        node_ = node_->prev_;
        return *this;
    }

    self operator++(int) {
        assert(node_ != nullptr && "invalid const_iterator increment");
        DoublyNodeConstReverseIterator<Item> tmp(*this);
        node_ = node_->prev_;
        return tmp;
    }

    // Comparison operators(const)
    bool operator==(const DoublyNodeConstReverseIterator<Item> &other) const {
        return (node_ == other.node_);
    }

    bool operator!=(const DoublyNodeConstReverseIterator<Item> &other) const {
        return !operator==(other);
    }

    // Swap operator
    void swap(DoublyNodeConstReverseIterator &other) {
        std::swap(node_, other.node_);
    }

private:
    const DoublyNode<Item> *node_;
};

#endif // LINKED_DOUBLY_NODE_H_
