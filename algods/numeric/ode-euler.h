//
// ode-euler.h
//
// Copyright (c) 2020 Carlos Braga
// This program is free software; you can redistribute it and/or modify it
// under the terms of the MIT License. See accompanying LICENSE.md or
// https://opensource.org/licenses/MIT.
//

#ifndef NUMERIC_ODE_EULER_H_
#define NUMERIC_ODE_EULER_H_

#include <omp.h>
#include <limits>
#include <cmath>
#include <type_traits>
#include "algods/numeric/array.h"

///
/// @brief Solve the set of odes using Euler algorithm:
///     dx
///     -- = f(x(t))
///     dt
///
///     x(t+dt) = x(t) + dt * f(x(t))
///
/// The time derivative of x at time t, f(x(t)), is computed from the derivative
/// operator Deriv. By default Deriv is a noop.
///
template<typename T = double, bool IsPar = true>
struct OdeEuler {
    static_assert(std::is_floating_point<T>::value, "non floating point type");

    OdeEuler() {}
    ~OdeEuler() = default;

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
void OdeEuler<T,IsPar>::Init(
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
/// @brief Integration step using explicit Euler algorithm.
///
template<typename T, bool IsPar>
template<typename Deriv>
bool OdeEuler<T,IsPar>::Step(
    Deriv deriv,
    Vector<T> & __restrict__ x,
    Vector<T> & __restrict__ dxdt,
    T dt)
{
    assert(x.n1() == dxdt.n1() && "invalid dimensions");

    // x(t+dt) = x(t) + dt * f(t)
    #pragma omp parallel for if(IsPar) default(none) shared(x, dxdt, dt) \
        schedule(static)
    for (size_t i = 0; i < x.n1(); i++) {
        x(i) += dt * dxdt(i);
    }

    // f(t+dt) = d(x(t+dt)) / dt
    deriv(x, dxdt);
    return true;
}

#endif // NUMERIC_ODE_EULER_H_
