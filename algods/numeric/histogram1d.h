//
// histogram1d.h
//
// Copyright (c) 2020 Carlos Braga
// This program is free software; you can redistribute it and/or modify it
// under the terms of the MIT License. See accompanying LICENSE.md or
// https://opensource.org/licenses/MIT.
//

#ifndef NUMERIC_HISTOGRAM1D_H_
#define NUMERIC_HISTOGRAM1D_H_

#include <cstring>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <memory>
#include <exception>
#include <type_traits>
#include <algorithm>
#include "algods/numeric/array.h"

struct Histogram1d {
    size_t xbins() const { return xbins_; }
    size_t size() const { return xbins_; }
    void clear();

    void set(double xmin, double xmax, bool reset = false);
    void bin(double xval, double weight, bool midpoint = true);

    double *countdata() { return count_.data(); }
    const double *countdata() const { return count_.data(); }
    double count(size_t i) const { return count_(i); }

    double *valuedata() { return value_.data(); }
    const double *valuedata() const { return value_.data(); }
    double value(size_t i) const { return value_(i); }

    double integral() const;

    void read(std::ifstream &ifs);
    void write(std::ofstream &ofs);

    Histogram1d clone() const;

    explicit Histogram1d(size_t xbins, double xmin, double xmax);
    ~Histogram1d() = default;

    Histogram1d(Histogram1d &&other);
    Histogram1d &operator=(Histogram1d &&other);

private:
    size_t xbins_{0};
    double xmin_{0.0};
    double xmax_{0.0};
    Vector<double> count_;
    Vector<double> value_;
};

///
/// @brief Create and 1d-histogram with the specified range.
///
inline Histogram1d::Histogram1d(size_t xbins, double xmin, double xmax)
    : xbins_(xbins)
    , xmin_(xmin)
    , xmax_(xmax)
    , count_(Vector<double>(xbins))
    , value_(Vector<double>(xbins))
{}

///
/// @brief Array move copy/assignment.
///
inline Histogram1d::Histogram1d(Histogram1d &&other)
    : xbins_(std::move(other.xbins_))
    , xmin_(std::move(other.xmin_))
    , xmax_(std::move(other.xmax_))
    , count_(std::move(other.count_))
    , value_(std::move(other.value_))
{}

inline Histogram1d &Histogram1d::operator=(Histogram1d &&other)
{
    if (this == &other) {
        return *this;
    }
    xbins_ = std::move(other.xbins_);
    xmin_ = std::move(other.xmin_);
    xmax_ = std::move(other.xmax_);
    count_ = std::move(other.count_);
    value_ = std::move(other.value_);
    return *this;
}

///
/// @brief Reset the histogram count and value data.
///
inline void Histogram1d::clear()
{
    count_.clear();
    value_.clear();
}

///
/// @brief Set the histogram ranges and reset if needed.
///
inline void Histogram1d::set(double xmin, double xmax, bool reset)
{
    xmin_ = xmin;
    xmax_ = xmax;
    if (reset) {
        clear();
    }
}

///
/// @brief Bin a sample value with a given weight.
///
inline void Histogram1d::bin(double xval, double weight, bool midpoint)
{
    // Increment the bin counter and bin value.
    // Compute the bin index using the lower or midpoint rounding conventions.
    // Lower convention takes the bin of x to be:
    //      [x_k, x_k + dx], i.e., k = floor((x-x_0)/dx).
    // Midpoint convention takes the bin of x to be:
    //      [x_k - dx/2, x_k + dx/2], i.e.,
    //      k = floor((x-(x0-dx/2)/dx) = floor((x-x0)/dx + dx/2).
    double bias = midpoint ? 0.5 : 0.0;
    double ux = (double) xbins_ * (xval - xmin_) / (xmax_ - xmin_);

    const size_t zero = static_cast<size_t>(0);
    const size_t xhi = static_cast<size_t>(xbins_ - 1);

    size_t ix = std::floor(ux + bias);
    ix = std::min(std::max(ix, zero), xhi);

    count_(ix) += 1.0;
    value_(ix) += weight;
}

///
/// @brief Compute the histogram probability integral.
///
inline double Histogram1d::integral() const
{
    double count = 0.0;
    for (size_t i = 0; i < xbins_; ++i) {
        count += count_(i);
    }

    double xdel = (xmax_ - xmin_) / (double) xbins_;
    return (count * xdel);
}

///
/// @brief Return a clone of the current histogram.
///
inline Histogram1d Histogram1d::clone() const
{
    Histogram1d histogram(xbins_, xmin_, xmax_);
    std::memcpy(histogram.countdata(), countdata(), size() * sizeof(double));
    std::memcpy(histogram.valuedata(), valuedata(), size() * sizeof(double));
    return std::move(histogram);
}

///
/// @brief Read the array elements from the input stream.
///
inline void Histogram1d::read(std::ifstream &ifs)
{
    if (!ifs) {
        throw std::runtime_error("invalid ifstream");
    }

    // Read histogram data.
    ifs >> xbins_;
    ifs >> xmin_;
    ifs >> xmax_;
    count_ = Vector<double>(xbins_);
    count_.read(ifs);
    value_ = Vector<double>(xbins_);
    value_.read(ifs);
}

///
/// @brief Write the array elements to the output stream.
///
inline void Histogram1d::write(std::ofstream &ofs)
{
    if (!ofs) {
        throw std::runtime_error("invalid ofstream");
    }

    ofs << std::setprecision(std::numeric_limits<double>::max_digits10);
    ofs << std::scientific;
    ofs << xbins_ << "\n";
    ofs << xmin_ << "\n";
    ofs << xmax_ << "\n";
    count_.write(ofs);
    value_.write(ofs);
}

#endif // NUMERIC_HISTOGRAM1D_H
