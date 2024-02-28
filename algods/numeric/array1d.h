//
// array1d.h
//
// Copyright (c) 2020 Carlos Braga
// This program is free software; you can redistribute it and/or modify it
// under the terms of the MIT License. See accompanying LICENSE.md or
// https://opensource.org/licenses/MIT.
//

#ifndef NUMERIC_ARRAY1D_H_
#define NUMERIC_ARRAY1D_H_

#include <fstream>
#include <iostream>
#include <iomanip>
#include <memory>
#include <exception>
#include <type_traits>

///
/// @brief Array1d represents a static 1d-array of n1 numeric elements.
///  API Array1d<T>:
///     Array1d(n1)                 create a array with n1 elements
///     ~Array1d()                  destroy the array
///
///     size_t size()               number of items in the array
///     size_t n1()                 number of items in the first dimension
///
///     T &operator()(i)            access indexed item(i)
///     const T &operator()(i)      const access indexed item(i)
///
///     T *data()                   access the data memory block
///     const T *data()             const access the data memory block
///
///     clone()                     return a clone of the array
///
///     void read(ifs)              read the array from an input stream
///     void write(ofs)             write the array to an output stream
///
template<typename T>
struct Array1d {
    static_assert(std::is_arithmetic<T>::value, "non arithmetic type");

    size_t n1() const { return n1_; }
    size_t size() const { return n1_; }

    T *data() { return data_.get(); }
    const T *data() const { return data_.get(); }

    T &operator()(size_t i) { return data_[i]; }
    const T &operator()(size_t i) const { return data_[i]; }

    void read(std::ifstream &ifs);
    void write(std::ofstream &ofs);

    void clear();
    Array1d clone();

    explicit Array1d(size_t n1, T value = static_cast<T>(0));
    ~Array1d() = default;

    Array1d(Array1d &&other);
    Array1d &operator=(Array1d &&other);

private:
    size_t n1_;                     // array dimensions
    std::unique_ptr<T[]> data_;     // array memory block
};

///
/// @brief Create and 1d-array with the specified size.
///
template<typename T>
inline Array1d<T>::Array1d(size_t n1, T value)
{
    n1_ = n1;
    data_ = std::make_unique<T[]>(n1);
    for (size_t i = 0; i < size(); ++i) {
        data_[i] = value;
    }
}

///
/// @brief Array move copy/assignment.
///
template<typename T>
inline Array1d<T>::Array1d(Array1d<T> &&other)
{
    n1_ = std::move(other.n1_);
    data_ = std::move(other.data_);
}

template<typename T>
inline Array1d<T> &Array1d<T>::operator=(Array1d<T> &&other)
{
    if (this == &other) {
        return *this;
    }
    n1_ = std::move(other.n1_);
    data_ = std::move(other.data_);
    return *this;
}

///
/// @brief Clear the contents of the current array.
///
template<typename T>
inline void Array1d<T>::clear()
{
    std::memset(data(), 0, size() * sizeof(T));
}

///
/// @brief Return a clone of the current array.
///
template<typename T>
inline Array1d<T> Array1d<T>::clone()
{
    Array1d<T> array(n1_);
    std::memcpy(array.data(), data(), size() * sizeof(T));
    return std::move(array);
}

///
/// @brief Read the array elements from the input stream.
///
template<typename T>
inline void Array1d<T>::read(std::ifstream &ifs)
{
    if (!ifs) {
        throw std::runtime_error("invalid ifstream");
    }

    // Read the array data.
    ifs >> n1_;
    data_ = std::make_unique<T[]>(n1_);

    size_t ix = 0;
    while (ifs && ix < size()) {
        ifs >> data_[ix++];
    }
}

///
/// @brief Write the array elements to the output stream.
///
template<typename T>
inline void Array1d<T>::write(std::ofstream &ofs)
{
    if (!ofs) {
        throw std::runtime_error("invalid ofstream");
    }

    ofs << std::setprecision(std::numeric_limits<T>::max_digits10);
    ofs << std::scientific;
    ofs << n1_ << "\n";

    size_t ix = 0;
    while (ofs && ix < size()) {
        ofs << data_[ix++] << "\n";
    }
}

#endif // NUMERIC_ARRAY1D_H_
