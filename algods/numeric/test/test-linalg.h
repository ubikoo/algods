//
// test-linalg.h
//
// Copyright (c) 2020 Carlos Braga
// This program is free software; you can redistribute it and/or modify it
// under the terms of the MIT License. See accompanying LICENSE.md or
// https://opensource.org/licenses/MIT.
//

#ifndef TEST_NUMERIC_LINALG_H_
#define TEST_NUMERIC_LINALG_H_

#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <random>
#include "algods/numeric/linalg.h"

/// ---- Test Gauss Solver ----------------------------------------------------
template<typename T>
void TestLinalgGauss(const size_t ndim)
{
    // Create random number generator.
    std::random_device seed;
    std::mt19937 rng(seed());
    std::uniform_real_distribution<double> urand(0.0,1.0);

    // Create a random matrix A, random vector B and solution vector X.
    Matrix<T> A(ndim, ndim);
    Vector<T> B(ndim);
    Vector<T> X(ndim);
    for (size_t i = 0; i < ndim; ++i) {
        // Sample the random matrix.
        for (size_t j = 0; j < ndim; ++j) {
            A(i,j) = urand(rng);
        }

        // Condition the main diagonal with a positive random number.
        A(i,i) += 1.0 + urand(rng);

        // Sample the random vector.
        B(i) = urand(rng);
    }

    // Solve the linear problem.
    Matrix<T> AA(ndim, ndim);
    Vector<T> BB(ndim);
    for (size_t i = 0; i < ndim; ++i) {
        for (size_t j = 0; j < ndim; ++j) {
            AA(i,j) = A(i,j);
        }
        BB(i) = B(i);
    }

    if (ndim < 4000) {
        GaussSolve<T>(AA, BB, X);
    } else {
        GaussSolveOmp<T>(AA, BB, X);
    }

    // Compute the error.
    Vector<T> C(ndim);
    MatmulVector<T>(A, X, C);

    double err = 0.0;
    for (size_t i = 0; i < ndim; ++i) {
        err += static_cast<double>(std::fabs(C(i) - B(i)));
    }
    err /= (double) ndim;
    std::cout << std::scientific << std::setprecision(21)
        << __FUNCTION__ << " error " << err << "\n";

    if (typeid(T) == typeid(double)) {
        REQUIRE(err < 1.0E-8);
    } else if (typeid(T) == typeid(float)) {
        REQUIRE(err < 1.0);
    }

    // Print the matrix
    //for (size_t i = 0; i < ndim; ++i) {
        //for (size_t j = 0; j < ndim; ++j) {
            //std::printf(" %g", A(i,j));
        //}
        //std::printf(" | %g %g %g\n", B(i), C(i), X(i));
    //}
}

/// ---- Test Jacobi Solver ---------------------------------------------------
template<typename T>
void TestLinalgJacobi(const size_t ndim)
{
    // Create random number generator.
    std::random_device seed;
    std::mt19937 rng(seed());
    std::uniform_real_distribution<double> urand(0.0,1.0);

    // Create a random matrix A
    Matrix<T> A(ndim, ndim);
    for (size_t i = 0; i < ndim; ++i) {
        for (size_t j = 0; j < ndim; ++j) {
            A(i,j) = urand(rng);
        }
        A(i,i) += 1.0;
    }

    // Make it symmetric by adding its transpose.
    Matrix<T> Atr(ndim, ndim);
    TransposeMatrix<T>(A, Atr);
    for (size_t i = 0; i < ndim; ++i) {
        for (size_t j = 0; j < ndim; ++j) {
            A(i,j) += Atr(i,j);
        }
    }

    // Solve the eigenvalue problem.
    Matrix<T> D(ndim, ndim);
    Matrix<T> V(ndim, ndim);

    T maxeps = std::sqrt(std::numeric_limits<T>::epsilon());
    size_t maxiter = 1000000;

    if (ndim < 160) {
        REQUIRE(EigenJacobi<T>(A, D, V, maxeps, maxiter));
    } else {
        REQUIRE(EigenJacobiOmp<T>(A, D, V, maxeps, maxiter));
    }
    std::cout << std::scientific << std::setprecision(21)
        << __FUNCTION__ << " converged to maxeps " << maxeps << "\n";

    // Sort the eigenvalues from largest to lowest.
    EigenSort<T>(D,V);

    // Compute the error.
    Matrix<T> AV(ndim, ndim);
    MatmulMatrix<T>(A, V, AV);

    Matrix<T> VD(ndim, ndim);
    MatmulMatrix<T>(V, D, VD);

    double err = 0.0;
    for (size_t i = 0; i < ndim; ++i) {
        for (size_t j = 0; j < ndim; ++j) {
            err += static_cast<double>(std::fabs(AV(i,j) - VD(i,j)));
        }
    }
    err /= (double) ndim;
    std::cout << std::scientific << std::setprecision(21)
        << __FUNCTION__ << " error " << err << "\n";

    if (typeid(T) == typeid(double)) {
        REQUIRE(err < 1.0E-4);
    } else if (typeid(T) == typeid(float)) {
        REQUIRE(err < 1.0);
    }
}

#endif // TEST_NUMERIC_LINALG_H_
