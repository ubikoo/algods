//
// deque-templ.h
//
// Copyright (c) 2020 Carlos Braga
// This program is free software; you can redistribute it and/or modify it
// under the terms of the MIT License. See accompanying LICENSE.md or
// https://opensource.org/licenses/MIT.
//

#ifndef ARRAYS_DEQUE_H_
#define ARRAYS_DEQUE_H_

///
/// @brief Deque represents a finite deque of (N >= 0) items stored in a
/// contiguous block of memory over which it can iterate.
/// The class implements a ring buffer and supports front/back insertion/removal
/// of items with constant O(1) complexity.
///
///  API Deque<Item>:
///      Deque()                create an empty deque
///      ~Deque()               remove all items and destroy the deque
///
///      bool is_empty()             is the deque empty?
///      size_t size()               number of items in the deque
///      void clear()                remove all items from the deque
///
///      void resize(size_t)         resize the number of items in the deque
///      void capacity()             max number of items in the deque
///
///      void push_front(item)       insert an item at the back of the deque
///      void pop_front()            remove an item from the back of the deque
///      void push_back(item)        insert an item at the back of the deque
///      void pop_back()             remove an item from the back of the deque
///
///      Item &back()                access the back item
///      const Item &back()          const access the back item
///      Item &front()               access the front item
///      const Item &front()         const access the front item
///
///      Item &operator[](i)        access indexed item
///      const Item &operator[](i)  const access indexed item
///
///      Item *data()               access the data memory block
///      const Item *data()         const access the data memory block
///
/// @see Algorithms, 4th Edition by Robert Sedgewick, Kevin Wayne.
///
/// @see https://en.wikipedia.org/wiki/Circular_buffer
///      https://embeddedartistry.com/blog/2017/05/17/creating-a-circular-buffer-in-c-and-c/
///      https://www.embedded.com/ring-buffer-basics
///      https://brilliant.org/wiki/modular-arithmetic
///      https://stackoverflow.com/questions/7594508
///
template<typename Item> class DequeIterator;
template<typename Item> class DequeConstIterator;
template<typename Item> class DequeReverseIterator;
template<typename Item> class DequeConstReverseIterator;

template<typename Item>
class Deque {
public:
    // Deque iterators.
    typedef DequeIterator<Item> iterator;
    typedef DequeConstIterator<Item> const_iterator;

    typedef DequeReverseIterator<Item> reverse_iterator;
    typedef DequeConstReverseIterator<Item> const_reverse_iterator;

    iterator begin() { return iterator(head_, capacity_, &data_[0]); }
    const_iterator begin() const {
        return const_iterator(head_, capacity_, &data_[0]);
    }

    iterator end() {
        size_t index = (tail_ + 1) % capacity_;
        return iterator(index, capacity_, &data_[0]);
    }
    const_iterator end() const {
        size_t index = (tail_ + 1) % capacity_;
        return const_iterator(index, capacity_, &data_[0]);
    }

    reverse_iterator rbegin() {
        return reverse_iterator(tail_, capacity_, &data_[0]);
    }
    const_reverse_iterator rbegin() const {
        return const_reverse_iterator(tail_, capacity_, &data_[0]);
    }

    reverse_iterator rend() {
        size_t index = (head_ + capacity_ - 1) % capacity_;
        return reverse_iterator(index, capacity_, &data_[0]);
    }
    const_reverse_iterator rend() const {
        size_t index = (head_ + capacity_ - 1) % capacity_;
        return const_reverse_iterator(index, capacity_, &data_[0]);
    }

    // Public member functions.
    bool is_empty() const;
    size_t size() const;
    void clear();

    void resize(size_t capacity);
    size_t capacity() const;

    void push_front(const Item &item);
    void pop_front();
    void push_back(const Item &item);
    void pop_back();

    Item &front();
    const Item &front() const;
    Item &back();
    const Item &back() const;

    Item &operator[](size_t ix);
    const Item &operator[](size_t ix) const;

    // Constructor/destructor.
    Deque();
    Deque(size_t capacity);
    ~Deque() { clear(); }

    // Copy constructor/assignment.
    Deque(const Deque &other);
    Deque<Item> &operator=(const Deque<Item> &other);

    // Move constructor/assignment.
    Deque(Deque &&other);
    Deque<Item> &operator=(Deque<Item> &&other);

private:
    // Private member variables.
    std::unique_ptr<Item[]> data_;          // deque of data items
    size_t capacity_;                       // deque item capacity
    size_t count_;                          // number of data items
    size_t head_;                           // index of the first item
    size_t tail_;                           // index of the last item

    // Private member functions.
    size_t advance(size_t &index);
    size_t retreat(size_t &index);
};

/// ---- Constructor and assignment operators ---------------------------------
/// @brief Copy the contents of the other array in the same order.
///
template<typename Item>
inline Deque<Item>::Deque()
    : data_(std::make_unique<Item[]>(1))
    , capacity_(1)
    , count_(0)
    , head_(0)
    , tail_(0)
{}

template<typename Item>
inline Deque<Item>::Deque(size_t capacity)
    : data_(std::make_unique<Item[]>(capacity))
    , capacity_(capacity)
    , count_(0)
    , head_(0)
    , tail_(0)
{}

template<typename Item>
inline Deque<Item>::Deque(const Deque<Item> &other)
    : Deque<Item>::Deque()
{
    for (const auto it : other) {
        push_back(it);
    }
}

///
/// @brief Clear the current items and copy the contents from the other array.
///
template<typename Item>
inline Deque<Item> &Deque<Item>::operator=(
    const Deque<Item> &other)
{
    if (this == &other) {
        return *this;
    }

    clear();
    for (const auto &it : other) {
        push_back(it);
    }

    return *this;
}

///
/// @brief Move the contents of the other object. Afterwards, the other is a
/// valid but unspecified object.
///
template<typename Item>
inline Deque<Item>::Deque(Deque<Item> &&other)
    : Deque<Item>::Deque()
{
    data_ = std::move(other.data_);
    capacity_ = std::move(other.capacity_);
    count_ = std::move(other.count_);
    head_ = std::move(other.head_);
    tail_ = std::move(other.tail_);
}

///
/// Clear the current items and move the contents of the other object.
/// Afterwards, the other is a valid but unspecified object.
///
template<typename Item>
inline Deque<Item> &Deque<Item>::operator=(Deque<Item> &&other)
{
    if (this == &other) {
        return *this;
    }
    clear();
    data_ = std::move(other.data_);
    capacity_ = std::move(other.capacity_);
    count_ = std::move(other.count_);
    head_ = std::move(other.head_);
    tail_ = std::move(other.tail_);
    return *this;
}

/// ---- Private functions ----------------------------------------------------
/// @brief Increment the index to the smallest positive value congruent
/// mod(capacity).
///
template<typename Item>
inline size_t Deque<Item>::advance(size_t &index)
{
    index = (index + 1) % capacity_;
    return index;
}

///
/// @brief Decrement the index to the smallest positive value congruent
/// mod(capacity).
///
template<typename Item>
inline size_t Deque<Item>::retreat(size_t &index)
{
    index = (index + capacity_ - 1) % capacity_;
    return index;
}

/// ---- Query functions ------------------------------------------------------
/// @brief Return true if the head doesn't point to an item.
///
template<typename Item>
inline bool Deque<Item>::is_empty() const
{
    return (count_ == 0);
}

///
/// @brief Return the number of items in the deque.
///
template<typename Item>
inline size_t Deque<Item>::size() const
{
    return count_;
}

///
/// @brief Reset the deque data with a single element deque and set the counters
/// to their default values.
///
template<typename Item>
inline void Deque<Item>::clear()
{
    data_ = std::make_unique<Item[]>(1);
    capacity_ = 1;
    count_ = 0;
    head_ = 0;
    tail_ = 0;
}

///
/// @brief Resize the deque to a new size. Copy the contents of the current
/// deque using the smallest positive index congruent mod(capacity) for each
/// item.
///
template<typename Item>
inline void Deque<Item>::resize(size_t capacity)
{
    assert(capacity > 0 && "invalid size");

    std::unique_ptr<Item[]> copy = std::make_unique<Item[]>(capacity);
    count_ = count_ < capacity ? count_ : capacity;
    for (size_t i = 0; i < count_; ++i) {
        size_t j = (head_ + i) % capacity_;
       copy[i] = data_[j];
    }
    capacity_ = capacity;
    data_.swap(copy);
    head_ = 0;
    tail_ = count_-1;
}

///
/// @brief Return the deque item capacity.
///
template<typename Item>
inline size_t Deque<Item>::capacity() const
{
    return capacity_;
}

/// ---- Mutator functions ----------------------------------------------------
/// @brief Create a new item at the front of the deque. If deque is empty, head
/// and tail already point to the new item. Otherwise, insert the item and
/// retreat the head to the smallest index value congruent mod(capacity).
///
template<typename Item>
inline void Deque<Item>::push_front(const Item &item)
{
    // Expand the array if we reached its capacity.
    if (count_ == capacity_) {
        resize(2 * capacity_);
    }

    // Retreat the head and add the item to the array.
   data_[retreat(head_)] = item;
    ++count_;
}

///
/// @brief Delete an item from the tail of the deque and retreat the tail to the
/// smallest positive index congruent mod(capacity).
///
template<typename Item>
inline void Deque<Item>::pop_front()
{
    assert(!is_empty() && "empty deque, out of range error");

    // Decrement and the deque item count. If this is the last item, there is
    // nothing else to update. Otherwise, advance the head pointer.
    if (--count_ == 0) {
        return;
    }
    advance(head_);

    // Shrink the array if half its size is twice the number of items, to
    // maintain the 2:1 ratio for capacity vs size.
    if (count_ == (capacity_ / 4)) {
        resize(capacity_ / 2);
    }
}
///
/// @brief Create a new item at the back of the deque. If deque is empty, head
/// and tail already point to the new item. Otherwise, insert the item and
/// advance the tail to the smallest index congruent mod(capacity).
///
template<typename Item>
inline void Deque<Item>::push_back(const Item &item)
{
    // Expand the array if we reached its capacity.
    if (count_ == capacity_) {
        resize(2 * capacity_);
    }

    // Advance the tail and add the item to the array.
   data_[advance(tail_)] = item;
    ++count_;
}

///
/// @brief Delete an item from the tail of the deque and retreat the
/// tail to the smallest positive index congruent mod(capacity).
///
template<typename Item>
inline void Deque<Item>::pop_back()
{
    assert(!is_empty() && "empty deque, out of range error");

    // Decrement and the deque item count. If this is the last item, there is
    // nothing else to update. Otherwise, retreat the tail pointer.
    if (--count_ == 0) {
        return;
    }
    retreat(tail_);

    // Shrink the array if half its size is twice the number of items, to
    // maintain the 2:1 ratio for capacity vs size.
    if (count_ == (capacity_ / 4)) {
        resize(capacity_ / 2);
    }
}

/// ---- Acessor functions ----------------------------------------------------
template<typename Item>
inline Item &Deque<Item>::front()
{
    assert(!is_empty() && "empty deque, out of range error");
    return (*this)[0];
}

template<typename Item>
inline const Item &Deque<Item>::front() const
{
    assert(!is_empty() && "empty deque, out of range error");
    return (*this)[0];
}

template<typename Item>
inline Item &Deque<Item>::back()
{
    assert(!is_empty() && "empty deque, out of range error");
    return (*this)[count_ - 1];
}

template<typename Item>
inline const Item &Deque<Item>::back() const
{
    assert(!is_empty() && "empty deque, out of range error");
    return (*this)[count_ - 1];
}

template<typename Item>
inline Item &Deque<Item>::operator[](size_t ix)
{
    assert(!is_empty() && "empty deque, out of range error");
    assert(ix < count_ && "out of range error");

    size_t jx = (head_ + ix) % capacity_;
    return data_[jx];
}

template<typename Item>
inline const Item &Deque<Item>::operator[](size_t ix) const
{
    assert(!is_empty() && "empty deque, out of range error");
    assert(ix < count_ && "out of range error");

    size_t jx = (head_ + ix) % capacity_;
    return data_[jx];
}

/// ---- Iterators ------------------------------------------------------------
/// @brief DequeIterator is a forward iterator, accessing the items in the
/// Deque in the direction from the head to the tail. It needs to implement
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
/// Or by defining explicitly a DequeConstIterator. The later is used here.
///
/// @see http://www.cplusplus.com/reference/iterator
///      https://en.cppreference.com/w/cpp/iterator
///      https://www.justsoftwaresolutions.co.uk/cplusplus/generating_sequences.html
///      https://codereview.stackexchange.com/questions/74609
///
template<typename Item>
class DequeIterator {
public:
    // Iterator member types.
    typedef DequeIterator<Item>    self;
    typedef Item                        value_type;
    typedef Item*                       pointer;
    typedef Item&                       reference;
    typedef ptrdiff_t                   difference_type;
    typedef std::forward_iterator_tag   iterator_category;

    // Constructor/destructor
    DequeIterator() : index_(0), capacity_(0), data_(nullptr) {}
    explicit DequeIterator(
        size_t index,
        size_t capacity,
        Item *data)
        : index_(index)
        , capacity_(capacity)
        , data_(data) {}
    ~DequeIterator() = default;

    // Dereference operators(const)
    reference operator*() const {
        assert(data_ != nullptr && "invalid iterator dereference");
        return data_[index_];
    }

    pointer operator->() const {
        assert(data_ != nullptr && "invalid iterator dereference");
        return &(operator*());
    }

    // Increment operators
    self operator++() {
        assert(data_ != nullptr && "invalid iterator increment");
        index_ = (index_ + 1) % capacity_;
        return *this;
    }

    self operator++(int) {
        assert(data_ != nullptr && "invalid iterator increment");
        size_t current = index_;
        index_ = (index_ + 1) % capacity_;
        return DequeIterator<Item>(current, capacity_, data_);
    }

    // Comparison operators(const)
    bool operator==(const DequeIterator<Item> &other) const {
        return (index_ == other.index_  &&
                capacity_ == other.capacity_ &&
                data_ == other.data_);
    }

    bool operator!=(const DequeIterator<Item> &other) const {
        return !operator==(other);
    }

    // Swap operator
    void swap(DequeIterator &other) {
        std::swap(index_, other.index_);
        std::swap(capacity_, other.capacity_);
        std::swap(data_, other.data_);
    }

private:
    size_t index_;
    size_t capacity_;
    Item *data_;
};

/// ---------------------------------------------------------------------------
/// @brief DequeConstIterator is a constant forward iterator. It is the
/// constant counterpart of DequeIterator.
///
template<typename Item>
class DequeConstIterator {
public:
    // Iterator member types.
    typedef DequeConstIterator<Item>   self;
    typedef Item                            value_type;
    typedef const Item*                     pointer;
    typedef const Item&                     reference;
    typedef ptrdiff_t                       difference_type;
    typedef std::forward_iterator_tag       iterator_category;

    // Constructor/destructor
    DequeConstIterator() : index_(0), capacity_(0), data_(nullptr) {}
    explicit DequeConstIterator(
        size_t index,
        size_t capacity,
        const Item *data)
        : index_(index)
        , capacity_(capacity)
        , data_(data) {}
    ~DequeConstIterator() = default;

    // Dereference operators(const)
    reference operator*() const {
        assert(data_ != nullptr && "invalid const_iterator dereference");
        return data_[index_];
    }

    pointer operator->() const {
        assert(data_ != nullptr && "invalid const_iterator dereference");
        return &(operator*());
    }

    // Increment operators
    self operator++() {
        assert(data_ != nullptr && "invalid const_iterator increment");
        index_ = (index_ + 1) % capacity_;
        return *this;
    }

    self operator++(int) {
        assert(data_ != nullptr && "invalid const_iterator increment");
        size_t current = index_;
        index_ = (index_ + 1) % capacity_;
        return DequeConstIterator<Item>(current, capacity_, data_);
    }

    // Comparison operators(const)
    bool operator==(const DequeConstIterator<Item> &other) const {
        return (index_ == other.index_  &&
                capacity_ == other.capacity_ &&
                data_ == other.data_);
    }

    bool operator!=(const DequeConstIterator<Item> &other) const {
        return !operator==(other);
    }

    // Swap operator
    void swap(DequeConstIterator &other) {
        std::swap(index_, other.index_);
        std::swap(capacity_, other.capacity_);
        std::swap(data_, other.data_);
    }

private:
    size_t index_;
    size_t capacity_;
    const Item *data_;
};

/// ---------------------------------------------------------------------------
/// @brief DequeReverseIterator is a reverse iterator, accessing the items
/// in the Deque in the direction from the tail to the head.
/// It needs to implement the following at least:
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
/// Or by defining explicitly a DequeConstIterator. The later is used here.
///
/// @see Algorithms, 4th Edition by Robert Sedgewick, Kevin Wayne.
///
/// http://www.cplusplus.com/reference/iterator
/// https://en.cppreference.com/w/cpp/iterator
/// https://www.justsoftwaresolutions.co.uk/cplusplus/generating_sequences.html
/// https://codereview.stackexchange.com/questions/74609
///
template<typename Item>
class DequeReverseIterator {
public:
    // Iterator member types.
    typedef DequeReverseIterator<Item> self;
    typedef Item                            value_type;
    typedef Item*                           pointer;
    typedef Item&                           reference;
    typedef ptrdiff_t                       difference_type;
    typedef std::forward_iterator_tag       iterator_category;

    // Constructor/destructor
    DequeReverseIterator() : index_(0), capacity_(0), data_(nullptr) {}
    explicit DequeReverseIterator(
        size_t index,
        size_t capacity,
        Item *data)
        : index_(index)
        , capacity_(capacity)
        , data_(data) {}
    ~DequeReverseIterator() = default;

    // Dereference operators(const)
    reference operator*() const {
        assert(data_ != nullptr && "invalid iterator dereference");
        return data_[index_];
    }

    pointer operator->() const {
        assert(data_ != nullptr && "invalid iterator dereference");
        return &(operator*());
    }

    // Increment operators
    self operator++() {
        assert(data_ != nullptr && "invalid iterator increment");
        index_ = (index_ + capacity_ - 1) % capacity_;
        return *this;
    }

    self operator++(int) {
        assert(data_ != nullptr && "invalid iterator increment");
        size_t current = index_;
        index_ = (index_ + capacity_ - 1) % capacity_;
        return DequeReverseIterator<Item>(index_--, capacity_, data_);
    }

    // Comparison operators(const)
    bool operator==(const DequeReverseIterator<Item> &other) const {
        return (index_ == other.index_  &&
                capacity_ == other.capacity_ &&
                data_ == other.data_);
    }

    bool operator!=(const DequeReverseIterator<Item> &other) const {
        return !operator==(other);
    }

    // Swap operator
    void swap(DequeReverseIterator &other) {
        std::swap(index_, other.index_);
        std::swap(capacity_, other.capacity_);
        std::swap(data_, other.data_);
    }

private:
    size_t index_;
    size_t capacity_;
    Item *data_;
};

/// ---------------------------------------------------------------------------
/// @brief DequeConstReverseIterator is a constant forward iterator.* It is
/// the constant counterpart of DequeReverseIterator.
///
template<typename Item>
class DequeConstReverseIterator {
public:
    // Iterator member types.
    typedef DequeConstReverseIterator<Item>    self;
    typedef Item                                    value_type;
    typedef const Item*                             pointer;
    typedef const Item&                             reference;
    typedef ptrdiff_t                               difference_type;
    typedef std::forward_iterator_tag               iterator_category;

    // Constructor/destructor
    DequeConstReverseIterator()
        : index_(0), capacity_(0), data_(nullptr) {}
    explicit DequeConstReverseIterator(
        size_t index,
        size_t capacity,
        const Item *data)
        : index_(index)
        , capacity_(capacity)
        , data_(data) {}
    ~DequeConstReverseIterator() = default;

    // Dereference operators(const)
    reference operator*() const {
        assert(data_ != nullptr && "invalid const_iterator dereference");
        return data_[index_];
    }

    pointer operator->() const {
        assert(data_ != nullptr && "invalid const_iterator dereference");
        return &(operator*());
    }

    // Increment operators
    self operator++() {
        assert(data_ != nullptr && "invalid const_iterator increment");
        index_ = (index_ + capacity_ - 1) % capacity_;
        return *this;
    }

    self operator++(int) {
        assert(data_ != nullptr && "invalid const_iterator increment");
        size_t current = index_;
        index_ = (index_ + capacity_ - 1) % capacity_;
        return DequeConstReverseIterator<Item>(current, capacity_, data_);
    }

    // Comparison operators(const)
    bool operator==(const DequeConstReverseIterator<Item> &other) const {
        return (index_ == other.index_  &&
                capacity_ == other.capacity_ &&
                data_ == other.data_);
    }

    bool operator!=(const DequeConstReverseIterator<Item> &other) const {
        return !operator==(other);
    }

    // Swap operator
    void swap(DequeConstReverseIterator &other) {
        std::swap(index_, other.index_);
        std::swap(capacity_, other.capacity_);
        std::swap(data_, other.data_);
    }

private:
    size_t index_;
    size_t capacity_;
    const Item *data_;
};

#endif // ARRAYS_DEQUE_H_
