//
// array3d.h
//
// Copyright (c) 2020 Carlos Braga
// This program is free software; you can redistribute it and/or modify it
// under the terms of the MIT License. See accompanying LICENSE.md or
// https://opensource.org/licenses/MIT.
//

#ifndef NUMERIC_ARRAY3D_H_
#define NUMERIC_ARRAY3D_H_

#include <fstream>
#include <iostream>
#include <iomanip>
#include <memory>
#include <exception>
#include <type_traits>

///
/// @brief Array3d represents a static 3d-array of n1*n2*n3 numeric elements.
///  API Array3d<T>:
///     Array3d(n1,n2,n3)           create a array with n1*n2*n3 elements
///     ~Array3d()                  destroy the array
///
///     size_t size()               number of items in the array
///     size_t n1()                 number of items in the first dimension
///     size_t n2()                 number of items in the second dimension
///     size_t n3()                 number of items in the third dimension
///
///     T &operator()(i,j,k)        access indexed item(i,j,k)
///     const T &operator()(i,j,k)  const access indexed item(i,j,k)
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
struct Array3d {
    static_assert(std::is_arithmetic<T>::value, "non arithmetic type");

    size_t n1() const { return n1_; }
    size_t n2() const { return n2_; }
    size_t n3() const { return n3_; }
    size_t size() const { return n1_ * n2_ * n3_; }

    T *data() { return data_.get(); }
    const T *data() const { return data_.get(); }

    T &operator()(size_t i, size_t j, size_t k) {
        return data_[i*n2_*n3_ + j*n3_ + k];
    }
    const T &operator()(size_t i, size_t j, size_t k) const {
        return data_[i*n2_*n3_ + j*n3_ + k];
    }

    void read(std::ifstream &ifs);
    void write(std::ofstream &ofs);

    void clear();
    Array3d clone();

    explicit Array3d(
        size_t n1,
        size_t n2,
        size_t n3,
        T value = static_cast<T>(0));
    ~Array3d() = default;

    Array3d(Array3d &&other);
    Array3d &operator=(Array3d &&other);

private:
    size_t n1_;                     // array dimensions
    size_t n2_;
    size_t n3_;
    std::unique_ptr<T[]> data_;     // array memory block
};

///
/// @brief Create and 1d-array with the specified size.
///
template<typename T>
inline Array3d<T>::Array3d(size_t n1, size_t n2, size_t n3, T value)
{
    n1_ = n1;
    n2_ = n2;
    n3_ = n3;
    data_ = std::make_unique<T[]>(n1 * n2 * n3);
    for (size_t i = 0; i < size(); ++i) {
        data_[i] = value;
    }
}

///
/// @brief Array move copy/assignment.
///
template<typename T>
inline Array3d<T>::Array3d(Array3d<T> &&other)
{
    n1_ = std::move(other.n1_);
    n2_ = std::move(other.n2_);
    n3_ = std::move(other.n3_);
    data_ = std::move(other.data_);
}

template<typename T>
inline Array3d<T> &Array3d<T>::operator=(Array3d<T> &&other)
{
    if (this == &other) {
        return *this;
    }
    n1_ = std::move(other.n1_);
    n2_ = std::move(other.n2_);
    n3_ = std::move(other.n3_);
    data_ = std::move(other.data_);
    return *this;
}

///
/// @brief Clear the contents of the current array.
///
template<typename T>
inline void Array3d<T>::clear()
{
    std::memset(data(), 0, size() * sizeof(T));
}

///
/// @brief Return a clone of the current array.
///
template<typename T>
inline Array3d<T> Array3d<T>::clone()
{
    Array3d<T> array(n1_, n2_, n3_);
    std::memcpy(array.data(), data(), size() * sizeof(T));
    return std::move(array);
}

///
/// @brief Read the array elements from the input stream.
///
template<typename T>
inline void Array3d<T>::read(std::ifstream &ifs)
{
    if (!ifs) {
        throw std::runtime_error("invalid ifstream");
    }

    // Read the array data.
    ifs >> n1_;
    ifs >> n2_;
    ifs >> n3_;
    data_ = std::make_unique<T[]>(n1_ * n2_ * n3_);

    size_t ix = 0;
    while (ifs && ix < size()) {
        ifs >> data_[ix++];
    }
}

///
/// @brief Write the array elements to the output stream.
///
template<typename T>
inline void Array3d<T>::write(std::ofstream &ofs)
{
    if (!ofs) {
        throw std::runtime_error("invalid ofstream");
    }

    ofs << std::setprecision(std::numeric_limits<T>::max_digits10);
    ofs << std::scientific;
    ofs << n1_ << "\n";
    ofs << n2_ << "\n";
    ofs << n3_ << "\n";

    size_t ix = 0;
    while (ofs && ix < size()) {
        ofs << data_[ix++] << "\n";
    }
}

#endif // NUMERIC_ARRAY3D_H_
