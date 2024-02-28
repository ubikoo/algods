//
// sde-ito.h
//
// Copyright (c) 2020 Carlos Braga
// This program is free software; you can redistribute it and/or modify it
// under the terms of the MIT License. See accompanying LICENSE.md or
// https://opensource.org/licenses/MIT.
//

#ifndef NUMERIC_SDE_ITO_H_
#define NUMERIC_SDE_ITO_H_

#include <omp.h>
#include <type_traits>
#include "algods/numeric/array.h"

///
/// @brief Solve the set of sdes of the type using Ito stochastic integrator:
///
///     dx(t) = f[x(t)]*dt + g[x(t)]*dw
///
/// where dt is the time step size, f(x) is the drift term and g(x) is the
/// diffusion term. The term dw is a Wiener process. It is characterized
/// by [w(t), t > 0] and has stationary independent increments dw, where
/// dw is a random Gaussian variable with mean zero and standard deviation
/// equal to sqrt(dt),
///
///     dw = sqrt(dt) * N(0,1)
///     <[dw]> = 0
///     <[dw]^2> = dt
///
/// The integration step uses the Ito interpration of stochastic differential
/// equations(Euler-Maruyama Method):
///
///     x(t+dt) = x(t) + f[x(t)]*dt + g[x(t)]*dw
///     dx = x(n+1) - x(n)
///     dt = t(n+1) - t(n)
///     dw = w(n+1) - w(n) = sqrt(dt)*N(0,1)
///
/// @see Stochastic algorithms for discontinuous multiplicative white noise,
/// Physical Review E, 81, 032104, 2010.
///
template<typename T = double,
         typename R = atto::math::rng::Kiss,
         bool IsPar = true>
struct SdeIto {
    static_assert(std::is_floating_point<T>::value, "non floating point type");

    SdeWiener<T,R,IsPar> wiener_;   // Wiener process generator
    size_t neq_;                    // number of equations
    Vector<T> dw_;                  // Wiener process

    SdeIto(size_t neq)
        : wiener_(SdeWiener<T,R,IsPar>())
        , neq_(neq)
        , dw_(Vector<T>(neq)) {}
    ~SdeIto() = default;

    SdeIto(SdeIto &&other)
        : wiener_(std::move(wiener_))
        , neq_(std::move(neq_))
        , dw_(std::move(other.dw_)) {}

    SdeIto &operator=(SdeIto &&other) {
        if (this == &other) {
            return *this;
        }
        wiener_ = std::move(wiener_);
        neq_ = std::move(neq_);
        dw_ = std::move(other.dw_);
        return *this;
    }

    template<typename Deriv>
    void Init(
        Deriv deriv,
        Vector<T> & __restrict__ x,
        Vector<T> & __restrict__ f,
        Vector<T> & __restrict__ g);

    template<typename Deriv>
    void Step(
        Deriv deriv,
        Vector<T> & __restrict__ x,
        Vector<T> & __restrict__ f,
        Vector<T> & __restrict__ g,
        T dt);
};

///
/// @brief Initialize the Sde solver. The derivative operator deriv is a functor
/// with two functions:
///     drift(x,f)  compute the drift derivative term in the sde, f[x(t)].
///     diff(x,g)   compute the diffusion derivative term in the sde, g[x(t)].
///
template<typename T,
         typename R,
         bool IsPar>
template<typename Deriv>
void SdeIto<T,R,IsPar>::Init(
    Deriv deriv,
    Vector<T> & __restrict__ x,
    Vector<T> & __restrict__ f,
    Vector<T> & __restrict__ g)
{
    assert(x.n1() == neq_ && f.n1() == neq_ && g.n1() == neq_ &&
        "invalid dimensions");

    // Drift and diffusion terms, f[x(t)] and g[x(t)]
    deriv.drift(x, f);
    deriv.diff(x, g);
}

///
/// @brief Sde integration stepper using the Ito interpration of stochastic
/// differential equations (Euler-Maruyama Method).
///
template<typename T,
         typename R,
         bool IsPar>
template<typename Deriv>
void SdeIto<T,R,IsPar>::Step(
    Deriv deriv,
    Vector<T> & __restrict__ x,
    Vector<T> & __restrict__ f,
    Vector<T> & __restrict__ g,
    T dt)
{
    assert(x.n1() == neq_ && f.n1() == neq_ && g.n1() == neq_ &&
        "invalid dimensions");

    // Wiener stochastic process, dw = sqrt(dt) * N(0,1)
    wiener_(dt, dw_);

    // Main time step, x(t+dt) = x(t) + (dt * f(t)) + (dw * g(t))
    #pragma omp parallel for if(IsPar) default(none) shared(x, f, g, dt) \
        schedule(static)
    for (size_t i = 0; i < neq_; i++) {
        x(i) += dt * f(i) + dw_(i) * g(i);
    } // omp parallel

    // Compute drift and diffusion terms, f[x(t+dt)] and g[x(t+dt)]
    deriv.drift(x, f);
    deriv.diff(x, g);
}

#endif // NUMERIC_SDE_ITO_H_
