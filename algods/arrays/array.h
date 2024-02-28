//
// array-templ.h
//
// Copyright (c) 2020 Carlos Braga
// This program is free software; you can redistribute it and/or modify it
// under the terms of the MIT License. See accompanying LICENSE.md or
// https://opensource.org/licenses/MIT.
//

#ifndef ARRAYS_ARRAY_H_
#define ARRAYS_ARRAY_H_

///
/// @brief Array maintains a templated array of N items stored in a
/// contiguous block of memory over which it can iterate. The class supports
/// back insertion/removal of items with constant O(1) complexity, and
/// forward/reverse iteration.
///
/// Front insertion/removal of items is not supported because it would involve
/// displacement and resizing of the array to accomodate the new item, resuting
/// in O(N) complexity.
///
///  API Array<Item>:
///      Array()                create an empty array
///      ~Array()               destroy the array and its items
///
///      bool is_empty()             is the array empty?
///      size_t size()               number of items in the array
///      void clear()                remove all items from the array
///
///      void resize(size_t)         resize the number of items in the array
///      void capacity()             max number of items in the array
///
///      void push_back(item)        insert an item at the back of the array
///      void pop_back()             remove an item from the back of the array
///
///      Item &front()               access the front item
///      const Item &front()         const access the front item
///      Item &back()                access the back item
///      const Item &back()          const access the back item
///
///      Item &operator[](i)        access indexed item
///      const Item &operator[](i)  const access indexed item
///
///      Item *data()               access the data memory block
///      const Item *data()         const access the data memory block
///
/// @see Algorithms, 4th Edition by Robert Sedgewick, Kevin Wayne.
///
template<typename Item> class ArrayIterator;
template<typename Item> class ArrayConstIterator;
template<typename Item> class ArrayReverseIterator;
template<typename Item> class ArrayConstReverseIterator;

template<typename Item>
class Array {
public:
    // Array iterators.
    typedef ArrayIterator<Item> iterator;
    typedef ArrayConstIterator<Item> const_iterator;

    typedef ArrayReverseIterator<Item> reverse_iterator;
    typedef ArrayConstReverseIterator<Item> const_reverse_iterator;

    iterator begin() { return iterator(&data_[0]); }
    const_iterator begin() const { return const_iterator(&data_[0]); }

    iterator end() { return iterator(&data_[0] + count_); }
    const_iterator end() const { return const_iterator(&data_[0] + count_); }

    reverse_iterator rbegin() { return reverse_iterator(&data_[0] + count_ - 1); }
    const_reverse_iterator rbegin() const {
        return const_reverse_iterator(&data_[0] + count_ - 1);
    }

    reverse_iterator rend() { return reverse_iterator(&data_[0] - 1); }
    const_reverse_iterator rend() const {
        return const_reverse_iterator(&data_[0] - 1);
    }

    // Public member functions.
    bool is_empty() const;
    size_t size() const;
    void clear();

    void resize(size_t capacity);
    size_t capacity() const;

    void push_back(const Item &item);
    void pop_back();

    Item &front();
    const Item &front() const;
    Item &back();
    const Item &back() const;

    Item &operator[](size_t ix);
    const Item &operator[](size_t ix) const;

    Item *data();
    const Item *data() const;

    // Array constructor/destructor.
    Array();
    Array(size_t capacity);
    ~Array() { clear(); }

    // Copy constructor/assignment.
    Array(const Array &other);
    Array<Item> &operator=(const Array<Item> &other);

    // Move constructor/assignment.
    Array(Array &&other);
    Array<Item> &operator=(Array<Item> &&other);

private:
    std::unique_ptr<Item[]> data_;          // array of data items
    size_t capacity_;                       // data items capacity
    size_t count_;                          // number of data items
};

/// ---- Constructor and assignment operators ---------------------------------
/// @brief Copy the contents of the other array in the same order.
///
template<typename Item>
inline Array<Item>::Array()
    : data_(std::make_unique<Item[]>(1))
    , capacity_(1)
    , count_(0)
{}

template<typename Item>
inline Array<Item>::Array(size_t capacity)
    : data_(std::make_unique<Item[]>(capacity))
    , capacity_(capacity)
    , count_(0)
{}

template<typename Item>
inline Array<Item>::Array(const Array<Item> &other)
    : Array<Item>::Array()
{
    for (const auto it : other) {
        push_back(it);
    }
}

///
/// @brief Clear the current items and copy the contents from the other array.
///
template<typename Item>
inline Array<Item> &Array<Item>::operator=(
    const Array<Item> &other)
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
inline Array<Item>::Array(Array<Item> &&other)
    : Array<Item>::Array()
{
    data_ = std::move(other.data_);
    capacity_ = std::move(other.capacity_);
    count_ = std::move(other.count_);
}

///
/// @brief Clear the current items and move the contents of the other object.
/// Afterwards, the other is a valid but unspecified object.
///
template<typename Item>
inline Array<Item> &Array<Item>::operator=(Array<Item> &&other)
{
    if (this == &other) {
        return *this;
    }
    clear();
    data_ = std::move(other.data_);
    capacity_ = std::move(other.capacity_);
    count_ = std::move(other.count_);
    return *this;
}

/// ----- Query functions -----------------------------------------------------
/// @brief Return true if the head doesn't point to an item.
///
template<typename Item>
inline bool Array<Item>::is_empty() const
{
    return (count_ == 0);
}

///
/// @brief Return the number of items in the array.
///
template<typename Item>
inline size_t Array<Item>::size() const
{
    return count_;
}

///
/// @brief Reset the array data with a single element array and reset counters
/// to their default values.
///
template<typename Item>
inline void Array<Item>::clear()
{
    data_ = std::make_unique<Item[]>(1);
    capacity_ = 1;
    count_ = 0;
}

///
/// @brief Resize the array to a new size. Create a copy array with the new
// capacity and copy the data items.
///
template<typename Item>
inline void Array<Item>::resize(size_t capacity)
{
    assert(capacity > 0 && "invalid size");

    std::unique_ptr<Item[]> copy = std::make_unique<Item[]>(capacity);
    count_ = count_ < capacity ? count_ : capacity;
    capacity_ = capacity;
    for (size_t i = 0; i < count_; ++i) {
       copy[i] = data_[i];
    }
    data_.swap(copy);
}

///
/// @brief Return the array item capacity.
///
template<typename Item>
inline size_t Array<Item>::capacity() const
{
    return capacity_;
}

/// ---- Mutator functions ----------------------------------------------------
/// @brief Create a new item at the tail of the array. If array is empty, update
/// head and tail to point to the new item. Otherwise, insert the item after the
/// head.
///
template<typename Item>
inline void Array<Item>::push_back(const Item &item)
{
    if (count_ == capacity_) {
        resize(2 * capacity_);
    }
    data_[count_] = item;
    count_++;
}

///
/// @brief Delete an item from the front of the array. Decrease the array item
/// count. If item count is (1/4) of the current capacity (1/4) then we can
/// halve the capacity and maintain a 2:1 ratio for capacity:size.
///
template<typename Item>
inline void Array<Item>::pop_back()
{
    assert(!is_empty() && "empty array, out of range error");
    --count_;
    if (count_ > 0 && count_ == (capacity_ / 4)) {
        resize(capacity_ / 2);
    }
}

///
/// @brief Return a reference to the first element in the array.
///
template<typename Item>
inline Item &Array<Item>::front()
{
    assert(!is_empty() && "empty array, out of range error");
    return (*this)[0];
}

template<typename Item>
inline const Item &Array<Item>::front() const
{
    assert(!is_empty() && "empty array, out of range error");
    return (*this)[0];
}

///
/// @brief Return a reference to the last element in the array.
///
template<typename Item>
inline Item &Array<Item>::back()
{
    assert(!is_empty() && "empty array, out of range error");
    return (*this)[count_ - 1];
}

template<typename Item>
inline const Item &Array<Item>::back() const
{
    assert(!is_empty() && "empty array, out of range error");
    return (*this)[count_ - 1];
}

///
/// @brief Return a reference to an item in the array.
///
template<typename Item>
inline Item &Array<Item>::operator[](size_t ix)
{
    assert(!is_empty() && "empty array, out of range error");
    assert(ix < count_ && "out of range error");
    return data_[ix];
}

template<typename Item>
inline const Item &Array<Item>::operator[](size_t ix) const
{
    assert(!is_empty() && "empty array, out of range error");
    assert(ix < count_ && "out of range error");
    return data_[ix];
}

///
/// @brief Return a pointer to the array data.
///
template<typename Item>
inline Item *Array<Item>::data()
{
    return data_.get();
}

///
/// @brief Return a pointer the array data.
///
template<typename Item>
inline const Item *Array<Item>::data() const
{
    return data_.get();
}

/// ---- Iterators ------------------------------------------------------------
/// @brief ArrayIterator as a forward iterator, accessing array items in
/// direction from head to tail. It needs to implement the following at least:
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
/// The implementation needs to define the following types describing the
/// iterator properties:
///
///  typedef typename Iterator::value_type            value_type;
///  typedef typename Iterator::pointer               pointer;
///  typedef typename Iterator::reference             reference;
///  typedef typename Iterator::difference_type       difference_type;
///  typedef typename Iterator::iterator_category     iterator_category;
///
/// Standard library containers provide two types of iterators:
///  iterator type - pointing to mutable data, and
///  const_iterator type - pointing to immutable data.
///
/// It is possible to support both by providing a conversion operator:
///  operator Iterator<const Item>() const {
///      return Iterator<const Item>(itr);
///  }
/// or by defining explicitly a ArrayConstIterator. The later is used here.
///
/// @see https://en.cppreference.com/w/cpp/iterator
///      http://www.cplusplus.com/reference/iterator
///      https://www.justsoftwaresolutions.co.uk/cplusplus/generating_sequences.html
///      https://codereview.stackexchange.com/questions/74609
///
template<typename Item>
class ArrayIterator {
public:
    // Iterator member types.
    typedef ArrayIterator<Item>    self;
    typedef Item                        value_type;
    typedef Item*                       pointer;
    typedef Item&                       reference;
    typedef ptrdiff_t                   difference_type;
    typedef std::forward_iterator_tag   iterator_category;

    // Constructor/destructor
    ArrayIterator() : item_(nullptr) {}
    explicit ArrayIterator(Item *item) : item_(item) {}
    ~ArrayIterator() = default;

    // Dereference operators(const)
    reference operator*() const {
        assert(item_ != nullptr && "invalid iterator dereference");
        return *item_;
    }

    pointer operator->() const {
        assert(item_ != nullptr && "invalid iterator dereference");
        return &(operator*());
    }

    // Increment operators
    self operator++() {
        assert(item_ != nullptr && "invalid iterator increment");
        ++item_;
        return *this;
    }

    self operator++(int) {
        assert(item_ != nullptr && "invalid iterator increment");
        return ArrayIterator<Item>(item_++);
    }

    // Comparison operators(const)
    bool operator==(const ArrayIterator<Item> &other) const {
        return (item_ == other.item_);
    }

    bool operator!=(const ArrayIterator<Item> &other) const {
        return !operator==(other);
    }

    // Swap operator
    void swap(ArrayIterator &other) {
        std::swap(item_, other.item_);
    }

private:
    Item *item_;
};

/// ---------------------------------------------------------------------------
/// @brief Constant counterpart of the forward iterator ArrayIterator.
///
template<typename Item>
class ArrayConstIterator {
public:
    // Iterator member types.
    typedef ArrayConstIterator<Item>   self;
    typedef Item                            value_type;
    typedef const Item*                     pointer;
    typedef const Item&                     reference;
    typedef ptrdiff_t                       difference_type;
    typedef std::forward_iterator_tag       iterator_category;

    // Constructor/destructor
    ArrayConstIterator() : item_(nullptr) {}
    explicit ArrayConstIterator(const Item *item) : item_(item) {}
    ~ArrayConstIterator() = default;

    // Dereference operators(const)
    reference operator*() const {
        assert(item_ != nullptr && "invalid const_iterator dereference");
        return *item_;
    }

    pointer operator->() const {
        assert(item_ != nullptr && "invalid const_iterator dereference");
        return &(operator*());
    }

    // Increment operators
    self operator++() {
        assert(item_ != nullptr && "invalid const_iterator increment");
        ++item_;
        return *this;
    }

    self operator++(int) {
        assert(item_ != nullptr && "invalid const_iterator increment");
        return ArrayConstIterator<Item>(item_++);
    }

    // Comparison operators(const)
    bool operator==(const ArrayConstIterator<Item> &other) const {
        return (item_ == other.item_);
    }

    bool operator!=(const ArrayConstIterator<Item> &other) const {
        return !operator==(other);
    }

    // Swap operator
    void swap(ArrayConstIterator &other) {
        std::swap(item_, other.item_);
    }

private:
    const Item *item_;
};

/// ---------------------------------------------------------------------------
/// @brief ArrayReverseIterator is a reverse iterator, accessing items in
/// Array from tail to head. It needs to implement the following at least:
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
/// The implementation needs to define the following types describing the
/// iterator properties:
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
/// Or by defining explicitly a ArrayConstIterator. The later is used here.
///
/// @see Algorithms, 4th Edition by Robert Sedgewick, Kevin Wayne.
///      http://www.cplusplus.com/reference/iterator
///      https://en.cppreference.com/w/cpp/iterator
///      https://www.justsoftwaresolutions.co.uk/cplusplus/generating_sequences.html
///      https://codereview.stackexchange.com/questions/74609
///
template<typename Item>
class ArrayReverseIterator {
public:
    // Iterator member types.
    typedef ArrayReverseIterator<Item> self;
    typedef Item                            value_type;
    typedef Item*                           pointer;
    typedef Item&                           reference;
    typedef ptrdiff_t                       difference_type;
    typedef std::forward_iterator_tag       iterator_category;

    // Constructor/destructor
    ArrayReverseIterator() : item_(nullptr) {}
    explicit ArrayReverseIterator(Item *item) : item_(item) {}
    ~ArrayReverseIterator() = default;

    // Dereference operators(const)
    reference operator*() const {
        assert(item_ != nullptr && "invalid iterator dereference");
        return *item_;
    }

    pointer operator->() const {
        assert(item_ != nullptr && "invalid iterator dereference");
        return &(operator*());
    }

    // Increment operators
    self operator++() {
        assert(item_ != nullptr && "invalid iterator increment");
        --item_;
        return *this;
    }

    self operator++(int) {
        assert(item_ != nullptr && "invalid iterator increment");
        return ArrayReverseIterator<Item>(item_--);
    }

    // Comparison operators(const)
    bool operator==(const ArrayReverseIterator<Item> &other) const {
        return (item_ == other.item_);
    }

    bool operator!=(const ArrayReverseIterator<Item> &other) const {
        return !operator==(other);
    }

    // Swap operator
    void swap(ArrayReverseIterator &other) {
        std::swap(item_, other.item_);
    }

private:
    Item *item_;
};

/// ---------------------------------------------------------------------------
/// @brief ArrayConstReverseIterator is a constant reverse iterator. It is
/// the constant counterpart of ArrayReverseIterator.
///
template<typename Item>
class ArrayConstReverseIterator {
public:
    // Iterator member types.
    typedef ArrayConstReverseIterator<Item>    self;
    typedef Item                                    value_type;
    typedef const Item*                             pointer;
    typedef const Item&                             reference;
    typedef ptrdiff_t                               difference_type;
    typedef std::forward_iterator_tag               iterator_category;

    // Constructor/destructor
    ArrayConstReverseIterator() : item_(nullptr) {}
    explicit ArrayConstReverseIterator(const Item *item)
        : item_(item) {}
    ~ArrayConstReverseIterator() = default;

    // Dereference operators(const)
    reference operator*() const {
        assert(item_ != nullptr && "invalid const_iterator dereference");
        return *item_;
    }

    pointer operator->() const {
        assert(item_ != nullptr && "invalid const_iterator dereference");
        return &(operator*());
    }

    // Increment operators
    self operator++() {
        assert(item_ != nullptr && "invalid const_iterator increment");
        --item_;
        return *this;
    }

    self operator++(int) {
        assert(item_ != nullptr && "invalid const_iterator increment");
        return ArrayConstReverseIterator<Item>(item_--);
    }

    // Comparison operators(const)
    bool operator==(const ArrayConstReverseIterator<Item> &other) const {
        return (item_ == other.item_);
    }

    bool operator!=(const ArrayConstReverseIterator<Item> &other) const {
        return !operator==(other);
    }

    // Swap operator
    void swap(ArrayConstReverseIterator &other) {
        std::swap(item_, other.item_);
    }

private:
    const Item *item_;
};

#endif // ARRAYS_ARRAY_H_
