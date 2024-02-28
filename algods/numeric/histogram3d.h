//
// histogram3d.h
//
// Copyright (c) 2020 Carlos Braga
// This program is free software; you can redistribute it and/or modify it
// under the terms of the MIT License. See accompanying LICENSE.md or
// https://opensource.org/licenses/MIT.
//

#ifndef NUMERIC_HISTOGRAM3D_H_
#define NUMERIC_HISTOGRAM3D_H_

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

struct Histogram3d {
    size_t xbins() const { return xbins_; }
    size_t ybins() const { return ybins_; }
    size_t zbins() const { return zbins_; }
    size_t size() const { return xbins_ * ybins_ * zbins_; }
    void clear();

    void set(
        double xmin,
        double xmax,
        double ymin,
        double ymax,
        double zmin,
        double zmax,
        bool reset = false);
    void bin(
        double xval,
        double yval,
        double zval,
        double weight,
        bool midpoint = true);

    double *countdata() { return count_.data(); }
    const double *countdata() const { return count_.data(); }
    double count(size_t i, size_t j, size_t k) const { return count_(i, j, k); }

    double *valuedata() { return value_.data(); }
    const double *valuedata() const { return value_.data(); }
    double value(size_t i, size_t j, size_t k) const { return value_(i, j, k); }

    double integral() const;

    void read(std::ifstream &ifs);
    void write(std::ofstream &ofs);

    Histogram3d clone() const;

    explicit Histogram3d(
        size_t xbins,
        size_t ybins,
        size_t zbins,
        double xmin,
        double xmax,
        double ymin,
        double ymax,
        double zmin,
        double zmax);
    ~Histogram3d() = default;

    Histogram3d(Histogram3d &&other);
    Histogram3d &operator=(Histogram3d &&other);

private:
    size_t xbins_{0};
    size_t ybins_{0};
    size_t zbins_{0};
    double xmin_{0.0};
    double xmax_{0.0};
    double ymin_{0.0};
    double ymax_{0.0};
    double zmin_{0.0};
    double zmax_{0.0};
    Tensor<double> count_;
    Tensor<double> value_;
};

///
/// @brief Create and 1d-histogram with the specified range.
///
inline Histogram3d::Histogram3d(
    size_t xbins,
    size_t ybins,
    size_t zbins,
    double xmin,
    double xmax,
    double ymin,
    double ymax,
    double zmin,
    double zmax)
    : xbins_(xbins)
    , ybins_(ybins)
    , zbins_(zbins)
    , xmin_(xmin)
    , xmax_(xmax)
    , ymin_(ymin)
    , ymax_(ymax)
    , zmin_(zmin)
    , zmax_(zmax)
    , count_(Tensor<double>(xbins, ybins, zbins))
    , value_(Tensor<double>(xbins, ybins, zbins))
{}

///
/// @brief Array move copy/assignment.
///
inline Histogram3d::Histogram3d(Histogram3d &&other)
    : xbins_(std::move(other.xbins_))
    , ybins_(std::move(other.ybins_))
    , zbins_(std::move(other.zbins_))
    , xmin_(std::move(other.xmin_))
    , xmax_(std::move(other.xmax_))
    , ymin_(std::move(other.ymin_))
    , ymax_(std::move(other.ymax_))
    , zmin_(std::move(other.zmin_))
    , zmax_(std::move(other.zmax_))
    , count_(std::move(other.count_))
    , value_(std::move(other.value_))
{}

inline Histogram3d &Histogram3d::operator=(Histogram3d &&other)
{
    if (this == &other) {
        return *this;
    }
    xbins_ = std::move(other.xbins_);
    ybins_ = std::move(other.ybins_);
    zbins_ = std::move(other.zbins_);
    xmin_ = std::move(other.xmin_);
    xmax_ = std::move(other.xmax_);
    ymin_ = std::move(other.ymin_);
    ymax_ = std::move(other.ymax_);
    zmin_ = std::move(other.zmin_);
    zmax_ = std::move(other.zmax_);
    count_ = std::move(other.count_);
    value_ = std::move(other.value_);
    return *this;
}

///
/// @brief Reset the histogram count and value data.
///
inline void Histogram3d::clear()
{
    count_.clear();
    value_.clear();
}

///
/// @brief Set the histogram ranges and reset if needed.
///
inline void Histogram3d::set(
    double xmin,
    double xmax,
    double ymin,
    double ymax,
    double zmin,
    double zmax,
    bool reset)
{
    xmin_ = xmin;
    xmax_ = xmax;
    ymin_ = ymin;
    ymax_ = ymax;
    zmin_ = zmin;
    zmax_ = zmax;
    if (reset) {
        clear();
    }
}

///
/// @brief Bin a sample value with a given weight.
///
inline void Histogram3d::bin(
    double xval,
    double yval,
    double zval,
    double weight,
    bool midpoint)
{
    // Increment the bin counter and bin value.
    // Compute the bin index using the lower or midpoint rounding conventions.
    // Lower convention takes the bin of x to be:
    //      [x_k, x_k + dx], i.e., k = floor((x-x_0)/dx).
    // Midpoint convention takes the bin of x to be:
    //      [x_k - dx/2, x_k + dx/2], i.e.,
    //      k = floor((x-(x0-dx/2)/dx) = floor((x-x0)/dx + dx/2).
    // The y- and z-directions are treated similarly.
    double bias = midpoint ? 0.5 : 0.0;

    double ux = (double) xbins_ * (xval - xmin_) / (xmax_ - xmin_);
    double uy = (double) ybins_ * (yval - ymin_) / (ymax_ - ymin_);
    double uz = (double) zbins_ * (zval - zmin_) / (zmax_ - zmin_);

    const size_t zero = static_cast<size_t>(0);
    const size_t xhi = static_cast<size_t>(xbins_ - 1);
    const size_t yhi = static_cast<size_t>(ybins_ - 1);
    const size_t zhi = static_cast<size_t>(zbins_ - 1);

    size_t ix = std::floor(ux + bias);
    size_t iy = std::floor(uy + bias);
    size_t iz = std::floor(uz + bias);
    ix = std::min(std::max(ix, zero), xhi);
    iy = std::min(std::max(iy, zero), yhi);
    iz = std::min(std::max(iz, zero), zhi);

    count_(ix, iy, iz) += 1.0;
    value_(ix, iy, iz) += weight;
}

///
/// @brief Compute the histogram probability integral.
///
inline double Histogram3d::integral() const
{
    double count = 0.0;
    for (size_t i = 0; i < xbins_; ++i) {
        for (size_t j = 0; j < ybins_; ++j) {
            for (size_t k = 0; k < zbins_; ++k) {
                count += count_(i, j, k);
            }
        }
    }

    double xdel = (xmax_ - xmin_) / (double) xbins_;
    double ydel = (ymax_ - ymin_) / (double) ybins_;
    double zdel = (zmax_ - zmin_) / (double) zbins_;
    return (count * xdel * ydel * ydel);
}

///
/// @brief Return a clone of the current histogram.
///
inline Histogram3d Histogram3d::clone() const
{
    Histogram3d histogram(xbins_, ybins_, zbins_, xmin_, xmax_, ymin_, ymax_,
        zmin_, zmax_);
    std::memcpy(histogram.countdata(), countdata(), size() * sizeof(double));
    std::memcpy(histogram.valuedata(), valuedata(), size() * sizeof(double));
    return std::move(histogram);
}

///
/// @brief Read the array elements from the input stream.
///
inline void Histogram3d::read(std::ifstream &ifs)
{
    if (!ifs) {
        throw std::runtime_error("invalid ifstream");
    }

    // Read histogram data.
    ifs >> xbins_;
    ifs >> ybins_;
    ifs >> zbins_;
    ifs >> xmin_;
    ifs >> xmax_;
    ifs >> ymin_;
    ifs >> ymax_;
    ifs >> zmin_;
    ifs >> zmax_;
    count_ = Tensor<double>(xbins_, ybins_, zbins_);
    count_.read(ifs);
    value_ = Tensor<double>(xbins_, ybins_, zbins_);
    value_.read(ifs);
}

///
/// @brief Write the array elements to the output stream.
///
inline void Histogram3d::write(std::ofstream &ofs)
{
    if (!ofs) {
        throw std::runtime_error("invalid ofstream");
    }

    ofs << std::setprecision(std::numeric_limits<double>::max_digits10);
    ofs << std::scientific;
    ofs << xbins_ << "\n";
    ofs << ybins_ << "\n";
    ofs << zbins_ << "\n";
    ofs << xmin_ << "\n";
    ofs << xmax_ << "\n";
    ofs << ymin_ << "\n";
    ofs << ymax_ << "\n";
    ofs << zmin_ << "\n";
    ofs << zmax_ << "\n";
    count_.write(ofs);
    value_.write(ofs);
}

#endif // NUMERIC_HISTOGRAM3D_H_
