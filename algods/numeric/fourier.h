//
// fourier.h
//
// Copyright (c) 2020 Carlos Braga
// This program is free software; you can redistribute it and/or modify it
// under the terms of the MIT License. See accompanying LICENSE.md or
// https://opensource.org/licenses/MIT.
//

#ifndef NUMERIC_FOURIER_H_
#define NUMERIC_FOURIER_H_

#include <omp.h>
#include "algods/numeric/array.h"

/// -----------------------------------------------------------------------------
/// @brief 1d forward dft transform of a real valued function to its complex
/// counterpart
///
///     f(n) -> (fre(n), fim(n))
///
/// and corresponding inverse transform
///
///     (fre(n), fim(n)) -> f(n)
///     n = number of waves
///     f = real function
///     fre = transform real part
///     fim = transform imaginary part
///
template<typename T = double, bool IsPar = true>
void Dft1dRow(
    size_t m1,
    Vector<T> & __restrict__ f,
    Vector<T> & __restrict__ fre,
    Vector<T> & __restrict__ fim)
{
    static_assert(std::is_floating_point<T>::value, "non floating point type");
    const size_t n = f.n1();
    double phi, k = 2.0 * M_PI / (T) n;

    fre(m1) = 0.0;
    fim(m1) = 0.0;
    for (size_t n1 = 0; n1 < n; ++n1) {
        phi = k * (T) m1 * n1;
        fre(m1) += f(n1) * std::cos(phi);
        fim(m1) -= f(n1) * std::sin(phi);
    }
}

template<typename T = double, bool IsPar = true>
void Dft1d(
    Vector<T> & __restrict__ f,
    Vector<T> & __restrict__ fre,
    Vector<T> & __restrict__ fim)
{
    static_assert(std::is_floating_point<T>::value, "non floating point type");
    assert(f.n1() == fre.n1() && f.n1() == fim.n1() && "invalid dimensions");

    const size_t n = f.n1();

    #pragma omp parallel for if(IsPar) default(none) shared(f, fre, fim) \
        schedule(static)
    for (size_t m1 = 0; m1 < n; ++m1) {
        Dft1dRow(m1, n, f, fre, fim);
    }
}

template<typename T = double, bool IsPar = true>
void InvDft1dRow(
    size_t n1,
    Vector<T> & __restrict__ fre,
    Vector<T> & __restrict__ fim,
    Vector<T> & __restrict__ f)
{
    static_assert(std::is_floating_point<T>::value, "non floating point type");
    const size_t n = f.n1();
    double phi, k = 2.0 * M_PI / (T) n;

    f(n1) = 0.0;
    for (size_t m1 = 0; m1 < n; ++m1) {
        phi = k * (T) m1 * n1;
        f(n1) += fre(m1) * std::cos(phi) - fim(m1) * std::sin(phi);
    }
    f(n1) /= (T) n;
}

template<typename T = double, bool IsPar = true>
void InvDft1d(
    Vector<T> & __restrict__ fre,
    Vector<T> & __restrict__ fim,
    Vector<T> & __restrict__ f)
{
    static_assert(std::is_floating_point<T>::value, "non floating point type");
    assert(f.n1() == fre.n1() && f.n1() == fim.n1() && "invalid dimensions");

    const size_t n = f.n1();

    #pragma omp parallel for if(IsPar) default(none) shared(f, fre, fim) \
        schedule(static)
    for (size_t n1 = 0; n1 < n; ++n1) {
        InvDft1dRow(n1, n, fre, fim, f);
    }
}

/// -----------------------------------------------------------------------------
/// @brief 2d forward dft transform of a real valued function to its complex
/// counterpart
///
///     f(n) -> (fre(n), fim(n))
///
/// and corresponding inverse transform
///
///     (fre(n), fim(n)) -> f(n)
///     n = number of waves
///     f = real function
///     fre = transform real part
///     fim = transform imaginary part
///
template<typename T = double, bool IsPar = true>
void Dft2dRow(
    size_t m1,
    Matrix<T> & __restrict__ f,
    Matrix<T> & __restrict__ fre,
    Matrix<T> & __restrict__ fim)
{
    static_assert(std::is_floating_point<T>::value, "non floating point type");
    const size_t n = f.n1();
    double phi, k = 2.0 * M_PI / (T) n;

    for (size_t m2 = 0; m2 < n; ++m2) {
        fre(m1, m2) = 0.0;
        fim(m1, m2) = 0.0;

        for (size_t n1 = 0; n1 < n; ++n1) {
            for (size_t n2 = 0; n2 < n; ++n2) {
                phi = k * (T) m1 * n1 +
                      k * (T) m2 * n2;
                fre(m1, m2) += f(n1, n2) * std::cos(phi);
                fim(m1, m2) -= f(n1, n2) * std::sin(phi);
            }
        }
    }
}

template<typename T = double, bool IsPar = true>
void Dft2d(
    Matrix<T> & __restrict__ f,
    Matrix<T> & __restrict__ fre,
    Matrix<T> & __restrict__ fim)
{
    static_assert(std::is_floating_point<T>::value, "non floating point type");
    assert(f.n1() == f.n2()   &&
           f.n1() == fre.n1() &&
           f.n2() == fre.n2() &&
           f.n1() == fim.n1() &&
           f.n2() == fim.n2() &&
           "invalid dimensions");

    const size_t n = f.n1();

    #pragma omp parallel for if(IsPar) default(none) shared(f, fre, fim) \
        schedule(static)
    for (size_t m1 = 0; m1 < n; ++m1) {
        Dft2dRow(m1, n, f, fre, fim);
    }
}

template<typename T = double, bool IsPar = true>
void InvDft2dRow(
    size_t n1,
    Matrix<T> & __restrict__ fre,
    Matrix<T> & __restrict__ fim,
    Matrix<T> & __restrict__ f)
{
    static_assert(std::is_floating_point<T>::value, "non floating point type");
    const size_t n = f.n1();
    double phi, k = 2.0 * M_PI / (T) n;

    for (size_t n2 = 0; n2 < n; ++n2) {
        f(n1, n2) = 0.0;
        for (size_t m1 = 0; m1 < n; ++m1) {
            for (size_t m2 = 0; m2 < n; ++m2) {
                phi = k * (T) m1 * n1 +
                      k * (T) m2 * n2;
                f(n1, n2) += fre(m1, m2) * std::cos(phi) -
                             fim(m1, m2) * std::sin(phi);
            }
        }
        f(n1, n2) /= (T) n * n;
    }
}

template<typename T = double, bool IsPar = true>
void InvDft2d(
    Matrix<T> & __restrict__ fre,
    Matrix<T> & __restrict__ fim,
    Matrix<T> & __restrict__ f)
{
    static_assert(std::is_floating_point<T>::value, "non floating point type");
    assert(f.n1() == f.n2()   &&
           f.n1() == fre.n1() &&
           f.n2() == fre.n2() &&
           f.n1() == fim.n1() &&
           f.n2() == fim.n2() &&
           "invalid dimensions");

    const size_t n = f.n1();

    #pragma omp parallel for if(IsPar) default(none) shared(f, fre, fim) \
        schedule(static)
    for (size_t n1 = 0; n1 < n; ++n1) {
        InvDft2dRow(n1, n, fre, fim, f);
    }
}

/// -----------------------------------------------------------------------------
/// @brief 3d forward dft transform of a real valued function to its complex
/// counterpart
///
///     f(n) -> (fre(n), fim(n))
///
/// and corresponding inverse transform
///
///     (fre(n), fim(n)) -> f(n)
///     n = number of waves
///     f = real function
///     fre = transform real part
///     fim = transform imaginary part
///
template<typename T = double, bool IsPar = true>
void Dft3dRow(
    size_t m1,
    Tensor<T> & __restrict__ f,
    Tensor<T> & __restrict__ fre,
    Tensor<T> & __restrict__ fim)
{
    static_assert(std::is_floating_point<T>::value, "non floating point type");
    const size_t n = f.n1();
    double phi, k = 2.0 * M_PI / (T) n;

    for (size_t m2 = 0; m2 < n; ++m2) {
        for (size_t m3 = 0; m3 < n; ++m3) {
            fre(m1, m2, m3) = 0.0;
            fim(m1, m2, m3) = 0.0;

            for (size_t n1 = 0; n1 < n; ++n1) {
                for (size_t n2 = 0; n2 < n; ++n2) {
                    for (size_t n3 = 0; n3 < n; ++n3) {
                        phi = k * (T) m1 * n1 +
                              k * (T) m2 * n2 +
                              k * (T) m3 * n3;
                        fre(m1, m2, m3) += f(n1, n2, n3) * std::cos(phi);
                        fim(m1, m2, m3) -= f(n1, n2, n3) * std::sin(phi);
                    }
                }
            }
        }
    }
}

template<typename T = double, bool IsPar = true>
void Dft3d(
    Tensor<T> & __restrict__ f,
    Tensor<T> & __restrict__ fre,
    Tensor<T> & __restrict__ fim)
{
    static_assert(std::is_floating_point<T>::value, "non floating point type");
    assert(f.n1() == f.n2()   &&
           f.n1() == f.n3()   &&
           f.n1() == fre.n1() &&
           f.n2() == fre.n2() &&
           f.n3() == fre.n3() &&
           f.n1() == fim.n1() &&
           f.n2() == fim.n2() &&
           f.n3() == fim.n3() &&
           "invalid dimensions");

    const size_t n = f.n1();

    #pragma omp parallel for if(IsPar) default(none) shared(f, fre, fim) \
        schedule(static)
    for (size_t m1 = 0; m1 < n; ++m1) {
        Dft3dRow(m1, n, f, fre, fim);
    }
}

template<typename T = double, bool IsPar = true>
void InvDtf3dRow(
    size_t n1,
    Tensor<T> & __restrict__ fre,
    Tensor<T> & __restrict__ fim,
    Tensor<T> & __restrict__ f)
{
    static_assert(std::is_floating_point<T>::value, "non floating point type");
    const size_t n = f.n1();
    double phi, k = 2.0 * M_PI / (T) n;

    for (size_t n2 = 0; n2 < n; ++n2) {
        for (size_t n3 = 0; n3 < n; ++n3) {
            f(n1, n2, n3) = 0.0;
            for (size_t m1 = 0; m1 < n; ++m1) {
                for (size_t m2 = 0; m2 < n; ++m2) {
                    for (size_t m3 = 0; m3 < n; ++m3) {
                        phi = k * (T) m1 * n1 +
                              k * (T) m2 * n2 +
                              k * (T) m3 * n3;
                        f(n1, n2, n3) += fre(m1, m2, m3) * std::cos(phi) -
                                         fim(m1, m2, m3) * std::sin(phi);
                    }
                }
            }
            f(n1, n2, n3) /= (T) n * n * n;
        }
    }
}

template<typename T = double, bool IsPar = true>
void InvDtf3d(
    Tensor<T> & __restrict__ fre,
    Tensor<T> & __restrict__ fim,
    Tensor<T> & __restrict__ f)
{
    static_assert(std::is_floating_point<T>::value, "non floating point type");
    assert(f.n1() == f.n2()   &&
           f.n1() == f.n3()   &&
           f.n1() == fre.n1() &&
           f.n2() == fre.n2() &&
           f.n3() == fre.n3() &&
           f.n1() == fim.n1() &&
           f.n2() == fim.n2() &&
           f.n3() == fim.n3() &&
           "invalid dimensions");

    const size_t n = f.n1();

    #pragma omp parallel for if(IsPar) default(none) shared(f, fre, fim) \
        schedule(static)
    for (size_t n1 = 0; n1 < n; ++n1) {
        InvDtf3dRow(n1, n, fre, fim, f);
    }
}

#endif // NUMERIC_FOURIER_H_
