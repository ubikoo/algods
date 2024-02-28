//
// spherical-harmonics.h
//
// Copyright (c) 2020 Carlos Braga
// This program is free software; you can redistribute it and/or modify it
// under the terms of the MIT License. See accompanying LICENSE.md or
// https://opensource.org/licenses/MIT.
//

#ifndef NUMERIC_SPHERICAL_HARMONICS_H_
#define NUMERIC_SPHERICAL_HARMONICS_H_

///
/// @brief Spherical harmonic transform implementation of a real valued function.
///
template<typename T = double>
struct SphHarm {
    static_assert(std::is_floating_point<T>::value, "non floating point type");

    typedef std::vector<std::vector<size_t>> VecIdx;
    typedef std::vector<std::vector<T>> VecVec;

    size_t lmax_;           // Legendre max order
    size_t nmax_;           // SphHarm max DoF

    VecIdx indx_;           // Legendre polynomical index
    VecVec Plm_;            // Legendre polynomial
    VecVec dxPlm_;          // Legendre polynomial derivative
    VecVec Klm_;            // spherical constants

    VecVec YlmCos_;         // Harmonic positive degree
    VecVec dthetaYlmCos_;   // Theta derivative of harmonic positive degree
    VecVec dphiYlmCos_;     // Phi derivative of harmonic positive degree

    VecVec YlmSin_;         // Harmonic negative degree
    VecVec dthetaYlmSin_;   // Theta derivative of harmonic negative degree
    VecVec dphiYlmSin_;     // Phi derivative of harmonic negative degree

    void Legendre(size_t lmax, T xval);
    void Compute(size_t lmax, T theta, T phi);

    // Constructor/destructor
    SphHarm(size_t lmax);
    ~SphHarm() = default;

    // Copy constructor/assignment.
    SphHarm(const SphHarm &other)
        : lmax_(lmax_)
        , nmax_(other.nmax_)
        , indx_(other.indx_)
        , Plm_(other.Plm_)
        , dxPlm_(other.dxPlm_)
        , Klm_(other.Klm_)
        , YlmCos_(other.YlmCos_)
        , dthetaYlmCos_(other.dthetaYlmCos_)
        , dphiYlmCos_(other.dphiYlmCos_)
        , YlmSin_(other.YlmSin_)
        , dthetaYlmSin_(other.dthetaYlmSin_)
        , dphiYlmSin_(other.dphiYlmSin_) {}
    SphHarm &operator=(const SphHarm &other) {
        if (this == &other) {
            return *this;
        }
        lmax_ = lmax_;
        nmax_ = other.nmax_;
        indx_  = other.indx_;
        Plm_   = other.Plm_;
        dxPlm_ = other.dxPlm_;
        Klm_   = other.Klm_;
        YlmCos_       = other.YlmCos_;
        dthetaYlmCos_ = other.dthetaYlmCos_;
        dphiYlmCos_   = other.dphiYlmCos_;
        YlmSin_       = other.YlmSin_;
        dthetaYlmSin_ = other.dthetaYlmSin_;
        dphiYlmSin_   = other.dphiYlmSin_;
        return *this;
    }

    // Move constructor/assignment.
    SphHarm(SphHarm &&other)
        : lmax_(std::move(lmax_))
        , nmax_(std::move(other.nmax_))
        , indx_(std::move(other.indx_))
        , Plm_(std::move(other.Plm_))
        , dxPlm_(std::move(other.dxPlm_))
        , Klm_(std::move(other.Klm_))
        , YlmCos_(std::move(other.YlmCos_))
        , dthetaYlmCos_(std::move(other.dthetaYlmCos_))
        , dphiYlmCos_(std::move(other.dphiYlmCos_))
        , YlmSin_(std::move(other.YlmSin_))
        , dthetaYlmSin_(std::move(other.dthetaYlmSin_))
        , dphiYlmSin_(std::move(other.dphiYlmSin_)) {}

    SphHarm &operator=(SphHarm &&other) {
        if (this == &other) {
            return *this;
        }
        lmax_  = std::move(lmax_);
        nmax_  = std::move(other.nmax_);
        indx_  = std::move(other.indx_);
        Plm_   = std::move(other.Plm_);
        dxPlm_ = std::move(other.dxPlm_);
        Klm_   = std::move(other.Klm_);
        YlmCos_       = std::move(other.YlmCos_);
        dthetaYlmCos_ = std::move(other.dthetaYlmCos_);
        dphiYlmCos_   = std::move(other.dphiYlmCos_);
        YlmSin_       = std::move(other.YlmSin_);
        dthetaYlmSin_ = std::move(other.dthetaYlmSin_);
        dphiYlmSin_   = std::move(other.dphiYlmSin_);
        return *this;
    }
};


///
/// @brief Generalized factorial with an arbitrary step value. Factorial and
/// double factorial are special cases with step=1 and step=2 respectivelly.
///
template<typename T>
T SphHarm::factorial(const T num, const T step = 1)
{
    static_assert(std::is_integral<T>::value, "non integer type");
    T n = num;
    T p = 1;
    while (n > 0) {
        p *= n;
        n -= step;
    }
    return p;
}

///
/// @brief Return the double factorial of a number, n!!
///
template<typename T>
T SphHarm::factorial2 (const T num)
{
    static_assert(std::is_integral<T>::value, "non integer type");
    return factorial(num, 2);
}

///
/// @brief Create a spherical harmonics object where lmax is the maximum order
/// of the Legendre Polynomial.
///
template<typename T>
SphHarm<T>::SphHarm(size_t lmax)
{
    static_assert(std::is_floating_point<T>::value, "non floating point type");
    assert(lmax > 0, "invalid spherical harmonics order");

    // Store the Legendre polynomial maximum order.
    lmax_ = lmax;

    // Create a vector of vectors of size 1 each.
    indx_  = std::move(VecIdx(lmax, 1));
    Plm_   = std::move(VecVec(lmax, 1));
    dxPlm_ = std::move(VecVec(lmax, 1));
    Klm_   = std::move(VecVec(lmax, 1));
    YlmCos_       = std::move(VecVec(lmax, 1));
    dthetaYlmCos_ = std::move(VecVec(lmax, 1));
    dphiYlmCos_   = std::move(VecVec(lmax, 1));
    YlmSin_       = std::move(VecVec(lmax, 1));
    dthetaYlmSin_ = std::move(VecVec(lmax, 1));
    dphiYlmSin_   = std::move(VecVec(lmax, 1));

    // Create the vectors for each Legendre polynomial order
    for (size_t l = 0; l < lmax; ++l) {
        indx_[l]  = std::move(Vector<size_t>(l+1));
        Plm_[l]   = std::move(Vector<T>(l+1));
        dxPlm_[l] = std::move(Vector<T>(l+1));
        Klm_[l]   = std::move(Vector<T>(l+1));
        YlmCos_[l]       = std::move(Vector<T>(l+1));
        dthetaYlmCos_[l] = std::move(Vector<T>(l+1));
        dphiYlmCos_[l]   = std::move(Vector<T>(l+1));
        YlmSin_[l]       = std::move(Vector<T>(l+1));
        dthetaYlmSin_[l] = std::move(Vector<T>(l+1));
        dphiYlmSin_[l]   = std::move(Vector<T>(l+1));
    }

    // Initialize the spherical harmonics array index
    nmax_ = 0;
    for (size_t l = 0; l < lmax; ++l) {
        for (size_t m = 0; m <= l; ++m) {
            indx_[l][m] = nmax_++;
        }
    }
}

///
/// @brief Compute the Sph Legendre polynomial
///
template<typename T>
void SphHarm<T>::Legendre(size_t lmax, T xval)
{
    assert(lmax <= lmax_, "invalid spherical harmonics order");
    assert(std::fabs(xval) < 1.0, "invalid range |x| >= 1");

    // Compute sph polynomial
    T r1 = 1.0 - xval * xval;

    // Compute P[l,l] elements
    for (size_t l = 0; l < lmax; l++) {
        T r2 = 0.5 * (T) l;
        Plm_[l][l] = std::pow(-1.0, l);
        Plm_[l][l] *= atto::math::factorial2(2 * l-1);
        Plm_[l][l] *= std::pow(r1, r2);
    }

    // Compute P[l,l-1] elements
    for (size_t l = 1; l < lmax; l++) {
        Plm_[l][l-1] = xval * (T) (2*l-1) * Plm_[l-1][l-1];
    }

    // Compute P[l+2,m] elements
    for (size_t l = 2; l < lmax; l++) {
        for (size_t m = 0; m <= l-2; m++) {
            T r1 = xval * (T) (2*l-1) / (T) (l-m);
            T r2 = (T) (l+m-1) / (T) (l-m);
            Plm_[l][m] = r1 * Plm_[l-1][m] - r2 * Plm_[l-2][m];
        }
    }

    // Compute the polynomial derivative dxPlm
    for (size_t l = 0; l < lmax; l++) {
        T r1 = xval * (T) (l) / (xval*xval - 1.0);
        dxPlm_[l][l] = r1 * Plm_[l][l];
    }

    for (size_t l = 1; l < lmax; l++) {
        for (size_t m = 0; m <= l-1; m++) {
            T r1 = xval * (T) (l) / (xval*xval - 1.0);
            T r2 = (T) (l+m) / (xval*xval - 1.0);
            dxPlm_[l][m] = r1 * Plm_[l][m] - r2 * Plm_[l-1][m];
        }
    }
}

///
/// @brief Compute the spherical harmonics.
///
template<typename T>
void SphHarm<T>::Compute(size_t lmax, T theta, T phi)
{
    assert(lmax <= lmax_, "invalid spherical harmonics order");

    constexpr T sqrt_two = std::sqrt(2.0);
    const T dcos_dtheta = -std::sin(theta);

    // Compute the Legendre polynomial Plm.
    Legendre(lmax, std::cos(theta));

    // Compute normalization constants
    constexpr double pi_4 = 4.0 * M_PI;
    for (size_t l = 0; l < lmax; l++) {
        for (size_t m = 0; m <= l; m++) {
            T r1  = (T) (2*l + 1) / pi_4;
            r1 *= Factorial(l - m) / Factorial(l + m);
            Klm_[l][m] = std::sqrt(r1);
        }
    }

    // Compute spherical harmonics positive degree(cosine based)
    for (size_t l = 0; l < lmax; l++) {
        YlmCos_[l](0) = Klm_[l](0) * Plm_[l](0);

        for (size_t m = 1; m <= l; m++) {
            T r1 = sqrt_two * Klm_[l][m] * std::cos((T) m * phi);
            YlmCos_[l][m] = r1 * Plm_[l][m];
        }
    }

    //  Compute the derivative spherical harmonics positive degree
    for (size_t l = 0; l < lmax; l++) {
        dthetaYlmCos_[l](0) = Klm_[l](0) * dxPlm_[l](0) * dcos_dtheta;

        for (size_t m = 1; m <= l; m++) {
            T r1 = sqrt_two * Klm_[l][m] * std::cos((T) m * phi) * dcos_dtheta;
            dthetaYlmCos_[l][m] = r1 * dxPlm_[l][m];
        }
    }

    for (size_t l = 0; l < lmax; l++) {
        dphiYlmCos_[l](0) = 0.0;

        for (size_t m = 1; m <= l; m++) {
            T r1 = -sqrt_two * Klm_[l][m] * (T) m * std::sin((T) m * phi);
            dphiYlmCos_[l][m] = r1 * Plm_[l][m];
        }
    }

    // Compute spherical harmonics of negative degree(sine based)
    for (size_t l = 0; l < lmax; l++) {
        YlmSin_[l](0) = Klm_[l](0) * Plm_[l](0);

        for (size_t m = 1; m <= l; m++) {
            T r1 = sqrt_two * Klm_[l][m] * std::sin((T) m * phi);
            YlmSin_[l][m] = r1 * Plm_[l][m];
        }
    }

    // Compute derivative spherical harmonics negative degree
    for (size_t l = 0; l < lmax; l++) {
        dthetaYlmSin_[l](0) = Klm_[l](0) * dxPlm_[l](0) * dcos_dtheta;

        for (size_t m = 1; m <= l; m++) {
            T r1 = sqrt_two * Klm_[l][m] * std::sin((T) m * phi) * dcos_dtheta;
            dthetaYlmSin_[l][m] = r1 * dxPlm_[l][m];
        }
    }

    for (size_t l = 0; l < lmax; l++) {
        dphiYlmSin_[l](0) = 0.0;

        for (size_t m = 1; m <= l; m++) {
            T r1 = sqrt_two * Klm_[l][m] * (T) m * std::cos((T) m * phi);
            dphiYlmSin_[l][m] = r1 * Plm_[l][m];
        }
    }
}

#endif // NUMERIC_SPHERICAL_HARMONICS_H_
