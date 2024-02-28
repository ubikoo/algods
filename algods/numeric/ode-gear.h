//
// ode-gear.h
//
// Copyright (c) 2020 Carlos Braga
// This program is free software; you can redistribute it and/or modify it
// under the terms of the MIT License. See accompanying LICENSE.md or
// https://opensource.org/licenses/MIT.
//

#ifndef NUMERIC_ODE_GEAR_H_
#define NUMERIC_ODE_GEAR_H_

#include <omp.h>
#include <limits>
#include <cmath>
#include <type_traits>
#include "algods/numeric/array.h"

///
/// @brief Solve the set of odes using Gear's Predictor Corrector 4th order
/// integrator:
///     dx
///     -- = f(x(t))
///     dt
///
///     x0(t) = dt^0 * (d(x(t)) / dt^0) / (0!) = x(t)
///     x1(t) = dt^1 * (d(x(t)) / dt^1) / (1!) = dt * (dx(t) / dt)
///     x2(t) = dt^2 * (d(x(t)) / dt^2) / (2!)
///     x3(t) = dt^3 * (d(x(t)) / dt^3) / (3!)
///     x4(t) = dt^4 * (d(x(t)) / dt^4) / (4!)
///
/// Predict system state at the next time step:
///     | xp0(t+dt) |   | 1 1 1 1 1 | | x0(t) |
///     | xp1(t+dt) |   | 0 1 2 3 4 | | x1(t) |
///     | xp2(t+dt) | = | 0 0 1 3 6 | | x2(t) |
///     | xp3(t+dt) |   | 0 0 0 1 4 | | x3(t) |
///     | xp4(t+dt) |   | 0 0 0 0 1 | | x4(t) |
///
/// Use the Predicted state to compute the time derivative the next time step:
///     f(t+dt) = f(xp0(t+dt))
///
/// If the integration was exact, the Predicted time derivative would be:
///     xp1(t+dt) = dt * f(t+dt)
///
/// The Corrector step aims to fix this difference in the form:
///     | xc0(t+dt) |   | xp0(t+dt) |   | c0 |
///     | xc1(t+dt) |   | xp1(t+dt) |   | c1 |
///     | xc2(t+dt) | = | xp2(t+dt) | - | c2 | dx(t+dt)
///     | xc3(t+dt) |   | xp3(t+dt) |   | c3 |
///     | xc4(t+dt) |   | xp4(t+dt) |   | c4 |
///
/// where dx is the difference between the Predicted time derivative at t+dt and
/// the Correct time derivative
///     dx(t+dt) = xp1(t+dt) - dt * f(t+dt)
///
/// The Gear Corrector coefficients are given by
///     c0 = 251/720
///     c1 = 1
///     c2 = 11/20
///     c3 = 1/3
///     c4 = 1/24
///
/// The time derivative of x at time t, f(x(t)), is computed from the derivative
/// operator Deriv. By default Deriv is a noop.
///
template<typename T = double, bool IsPar = true>
struct OdeGear {
    static_assert(std::is_floating_point<T>::value, "non floating point type");

    size_t neq_;        // system dimension
    Vector<T> x0_;      // Gear states
    Vector<T> x1_;
    Vector<T> x2_;
    Vector<T> x3_;
    Vector<T> x4_;

    // Constructor/destructor
    OdeGear(size_t neq)
        : neq_(neq)
        , x0_(Vector<T>(neq))
        , x1_(Vector<T>(neq))
        , x2_(Vector<T>(neq))
        , x3_(Vector<T>(neq))
        , x4_(Vector<T>(neq)) {}
    ~OdeGear() = default;

    OdeGear(OdeGear &&other)
        : neq_(std::move(neq_))
        , x0_(std::move(other.x0_))
        , x1_(std::move(other.x1_))
        , x2_(std::move(other.x2_))
        , x3_(std::move(other.x3_))
        , x4_(std::move(other.x4_)) {}

    OdeGear &operator=(OdeGear &&other) {
        if (this == &other) {
            return *this;
        }
        neq_ = std::move(neq_);
        x0_ = std::move(other.x0_);
        x1_ = std::move(other.x1_);
        x2_ = std::move(other.x2_);
        x3_ = std::move(other.x3_);
        x4_ = std::move(other.x4_);
        return *this;
    }

    //
    // Ode functions.
    //
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
    void Predict(Vector<T> & __restrict__ x);
    void Correct(
        Vector<T> & __restrict__ x,
        Vector<T> & __restrict__ dxdt,
        T dt);
};

///
/// @brief Initialize the solver internal state.
///
template<typename T, bool IsPar>
template<typename Deriv>
void OdeGear<T,IsPar>::Init(
    Deriv deriv,
    Vector<T> & __restrict__ x,
    Vector<T> & __restrict__ dxdt,
    T dt)
{
    assert(x.n1() == neq_ && dxdt.n1() == neq_ && "invalid dimensions");

    // f(t) = d(x(t)) / dt
    deriv(x, dxdt);

    // x0 = x(t)
    // x1 = dt * d(x(t)) / dt
    #pragma omp parallel for if(IsPar) default(none) shared(x, dxdt, dt) \
        schedule(static)
    for (size_t i = 0; i < neq_; i++) {
        x0_(i) = x(i);
        x1_(i) = dt * dxdt(i);
    }
}

///
/// @brief Integration step using explicit Gear algorithm.
///
template<typename T, bool IsPar>
template<typename Deriv>
bool OdeGear<T,IsPar>::Step(
    Deriv deriv,
    Vector<T> & __restrict__ x,
    Vector<T> & __restrict__ dxdt,
    T dt)
{
    assert(x.n1() == neq_ && dxdt.n1() == neq_ && "invalid dimensions");

    // Predict system state at t+dt
    Predict(x);

    // Time derivative f(t+dt) = d(x(t+dt)) / dt
    deriv(x, dxdt);

    // Correct system state at t+dt
    Correct(x, dxdt, dt);
    return true;
}

///
/// @brief Gear Predictor integration step.
///
template<typename T, bool IsPar>
void OdeGear<T,IsPar>::Predict(Vector<T> & __restrict__ x)
{
    assert(x.n1() == neq_ && "invalid dimensions");

    // Predict system state at t+dt
    #pragma omp parallel if(IsPar) default(none) shared(x)
    {
        #pragma omp for schedule(static)
        for (size_t i = 0; i < neq_; i++) {
            x0_(i) = x(i);
        }

        #pragma omp for schedule(static)
        for (size_t i = 0; i < neq_; i++) {
            double x0 = x0_(i);
            double x1 = x1_(i);
            double x2 = x2_(i);
            double x3 = x3_(i);
            double x4 = x4_(i);

            x0_(i) = x0 + x1 +     x2 +     x3 +     x4;
            x1_(i) =      x1 + 2.0*x2 + 3.0*x3 + 4.0*x4;
            x2_(i) =               x2 + 3.0*x3 + 6.0*x4;
            x3_(i) =                        x3 + 4.0*x4;
            x4_(i) =                                 x4;
        }

        #pragma omp for schedule(static)
        for (size_t i = 0; i < neq_; i++) {
             x(i) = x0_(i);
        }
    } // omp parallel
}

///
/// @brief Gear Corrector integration step.
///
template<typename T, bool IsPar>
void OdeGear<T,IsPar>::Correct(
    Vector<T> & __restrict__ x,
    Vector<T> & __restrict__ dxdt,
    T dt)
{
    assert(x.n1() == neq_ && dxdt.n1() == neq_ && "invalid dimensions");

    constexpr T c0 = 251.0 / 720.0;
    constexpr T c1 = 1.0;
    constexpr T c2 = 11.0 / 12.0;
    constexpr T c3 = 1.0 / 3.0;
    constexpr T c4 = 1.0 / 24.0;

    // Correct system state at t+dt
    #pragma omp parallel if(IsPar) default(none) shared(x, dxdt, dt)
    {
        #pragma omp for schedule(static)
        for (size_t i = 0; i < neq_; i++) {
            double delx = x1_(i) - dt * dxdt(i);
            x0_(i) -= delx * c0;
            x1_(i) -= delx * c1;
            x2_(i) -= delx * c2;
            x3_(i) -= delx * c3;
            x4_(i) -= delx * c4;
        }

        #pragma omp for schedule(static)
        for (size_t i = 0; i < neq_; i++) {
            x(i) = x0_(i);
        }
    } // omp parallel
}

#endif // NUMERIC_ODE_GEAR_H_
