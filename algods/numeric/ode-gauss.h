//
// ode-gauss.h
//
// Copyright (c) 2020 Carlos Braga
// This program is free software; you can redistribute it and/or modify it
// under the terms of the MIT License. See accompanying LICENSE.md or
// https://opensource.org/licenses/MIT.
//

#ifndef NUMERIC_ODE_GAUSS_H_
#define NUMERIC_ODE_GAUSS_H_

#include <omp.h>
#include <limits>
#include <cmath>
#include <type_traits>
#include "algods/numeric/array.h"

///
/// @brief Solve the set of odes using implicit Gauss-Legendre algorithm:
///     dx
///     -- = f(x(t))
///     dt
///
///     x(t+dt) = x(t) + dt * f(t+dt/2, (x(t), x(t+dt)/2)
///
/// For each time step, the solution is computed using a fixed point iteration.
/// The error and the and max number of iteration steps are part of the
/// integrator state.
/// The time derivative of x at time t, f(x(t)), is computed from the derivative
/// operator Deriv. By default Deriv is a noop.
///
template<typename T = double, bool IsPar = true>
struct OdeGauss {
    static_assert(std::is_floating_point<T>::value, "non floating point type");

    size_t neq_;        // system dimension
    double maxerr_;     // max error value
    size_t maxiter_;    // max number of iterations
    Vector<T> zval_;
    Vector<T> zmid_;
    Vector<T> znew_;
    Vector<T> dzdt_;

    // Constructor/destructor
    OdeGauss(
        size_t neq,
        double maxerr = std::sqrt(std::numeric_limits<double>::epsilon()),
        size_t maxiter = 16)
        : neq_(neq)
        , maxerr_(maxerr)
        , maxiter_(maxiter)
        , zval_(Vector<T>(neq))
        , zmid_(Vector<T>(neq))
        , znew_(Vector<T>(neq))
        , dzdt_(Vector<T>(neq)) {}
    ~OdeGauss() = default;

    OdeGauss(OdeGauss &&other)
        : neq_(std::move(neq_))
        , maxerr_(std::move(other.maxerr_))
        , maxiter_(std::move(other.maxiter_))
        , zval_(std::move(other.zval_))
        , zmid_(std::move(other.zmid_))
        , znew_(std::move(other.znew_))
        , dzdt_(std::move(other.dzdt_)) {}

    OdeGauss &operator=(OdeGauss &&other) {
        if (this == &other) {
            return *this;
        }
        neq_ = std::move(neq_);
        maxerr_ = std::move(other.maxerr_);
        maxiter_ = std::move(other.maxiter_);
        zval_ = std::move(other.zval_);
        zmid_ = std::move(other.zmid_);
        znew_ = std::move(other.znew_);
        dzdt_ = std::move(other.dzdt_);
        return *this;
    }

    // Ode functions.
    template<typename Deriv>
    void Init(
        Deriv deriv,
        Vector<T> & __restrict__ x,
        Vector<T> & __restrict__ dxdt,
        T dt);

    template<typename Deriv>
    bool Step(
        Deriv deriv,
        Vector<T> & __restrict__ x,
        Vector<T> & __restrict__ dxdt,
        T dt);
};

///
/// @brief Initialize the solver internal state.
///
template<typename T, bool IsPar>
template<typename Deriv>
void OdeGauss<T,IsPar>::Init(
    Deriv deriv,
    Vector<T> & __restrict__ x,
    Vector<T> & __restrict__ dxdt,
    T dt)
{
    assert(x.n1() == dxdt.n1() && "invalid dimensions");

    // f(t) = d(x(t)) / dt
    deriv(x, dxdt);
}

///
/// @brief Integration step using explicit Gauss algorithm.
///
template<typename T, bool IsPar>
template<typename Deriv>
bool OdeGauss<T,IsPar>::Step(
    Deriv deriv,
    Vector<T> & __restrict__ x,
    Vector<T> & __restrict__ dxdt,
    T dt)
{
    assert(x.n1() == dxdt.n1() && "invalid dimensions");

    // Initialize the half step vector zval.
    #pragma omp parallel for if(IsPar) default(none) schedule(static)
    for (size_t i = 0; i < neq_; i++) {
        zval_(i) = 0.0;
    }

    // Solve the implicit integration step using a fixed point algorithm.
    double err = std::numeric_limits<T>::max();
    size_t iter = 0;
    while (err > maxerr_ && ++iter < maxiter_) {
        // Compute time derivative at zmid = (x(t+dt) + x(t)) / 2
        #pragma omp parallel for if(IsPar) default(none) shared(x) \
            schedule(static)
        for (size_t i = 0; i < neq_; i++) {
            zmid_(i) = zval_(i) + x(i);
        }
        deriv(zmid_, dzdt_);

        // Update iterate value znew = 0.5 * dt * dzdt
        #pragma omp parallel for if(IsPar) default(none) shared(dt) \
            schedule(static)
        for (size_t i = 0; i < neq_; i++) {
            znew_(i) = 0.5 * dt * dzdt_(i);
        }

        // Check convergence and loop
        err = 0.0;
        #pragma omp parallel if(IsPar) default(none) shared(err)
        {
            double thr_err = 0.0;

            #pragma omp for schedule(static)
            for (size_t i = 0; i < neq_; ++i) {
                thr_err += std::fabs(znew_(i) - zval_(i));
                zval_(i) = znew_(i);
            }

            #pragma omp critical
            {
                err += thr_err;
            } // omp critical
        } // omp parallel
    }

    // Integrate the state and update time derivatives
    //      x(t+dt) = x(t) + 2.0 * z
    //      f(t+dt) = d(x(t+dt)) / dt
    #pragma omp parallel for if(IsPar) default(none) shared(x) \
        schedule(static)
    for (size_t i = 0; i < neq_; i++) {
        x(i) += 2.0 * zval_(i);
    }
    deriv(x, dxdt);

    // Check convergence of the Gauss fixed point iteration.
    return (iter < maxiter_);
}

#endif // NUMERIC_ODE_GAUSS_H_
