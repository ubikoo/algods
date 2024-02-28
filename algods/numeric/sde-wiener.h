//
// sde-wiener.h
//
// Copyright (c) 2020 Carlos Braga
// This program is free software; you can redistribute it and/or modify it
// under the terms of the MIT License. See accompanying LICENSE.md or
// https://opensource.org/licenses/MIT.
//

#ifndef NUMERIC_SDE_WIENER_H_
#define NUMERIC_SDE_WIENER_H_

#include <cmath>
#include <limits>
#include <random>
#include <type_traits>
#include "algods/numeric/array.h"

///
/// @brief Compute a Wiener stochastic process
///     dw = sqrt(dt) * N(0,1)
///
/// where N(0,1) is a Gaussian distribution with zero mean and unit variance.
/// The Wiener generator is an array of random number generators that run in
/// parallel on each thread available.
///
template<typename T = double,
         bool IsPar = true>
struct SdeWiener {
    static_assert(std::is_floating_point<T>::value, "non floating point type");

    std::random_device seed_;
    std::vector<std::mt19937> engine_;
    std::vector<std::uniform_real_distribution<T>> random_;

    SdeWiener() {
        #pragma omp parallel default(none) shared(engine_, random_)
        {
            #pragma omp master
            {
                int num_threads = omp_get_num_threads();

                engine_.resize(num_threads);
                for (auto &e : engine_) {
                    e.seed(seed_());
                }

                for (size_t i = 0; i < num_threads; ++i) {
                    auto rand = std::uniform_real_distribution<T>(0.0, 1.0);
                    random_.push_back(rand);
                }
            } // omp master
        } // omp parallel
    }
    ~SdeWiener() = default;

    void operator()(T dt, Vector<T> & __restrict__ dw);
};

template<typename T, typename R, bool IsPar>
void SdeWiener<T,R,IsPar>::operator()(T dt, Vector<T> & __restrict__ dw)
{
    // dw = sqrt(dt) * N(0,1)
    #pragma omp parallel if(IsPar) default(none) shared(dt, dw)
    {
        size_t tid = omp_get_thread_num();
        T dt2 = std::sqrt(dt);

        #pragma omp for schedule(static)
        for (size_t i = 0; i < dw.n1(); i++) {
            auto &rand = random_[tid];
            auto &eng = engine_[tid];
            dw(i) = dt2 * rand(eng);
        }
    } // omp parallel
}

#endif // NUMERIC_SDE_WIENER_H_
