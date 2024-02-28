//
// linalg-gauss.h
//
// Copyright (c) 2020 Carlos Braga
// This program is free software; you can redistribute it and/or modify it
// under the terms of the MIT License. See accompanying LICENSE.md or
// https://opensource.org/licenses/MIT.
//

#ifndef NUMERIC_LINALG_GAUSS_H_
#define NUMERIC_LINALG_GAUSS_H_

#include <omp.h>
#include "algods/numeric/array.h"

///
/// @brief Gauss solver support function. Subtract row A(i,:) from A(j,:) and
/// B(i) from B(j):
///     j   row index below diagonal
///     i   diagonal row index
///     m   number of columns
///     A   matrix
///     B   vector
///
template<typename T = double>
void GaussSubtractRow(
    size_t j,
    size_t i,
    size_t m,
    Matrix<T> & __restrict__ A,
    Vector<T> & __restrict__ B)
{
    static_assert(std::is_floating_point<T>::value, "non floating point type");
    // Subtract row-i premultiplied by A(j,i) from row-j.
    // @note Row subtraction is computationaly heavy. Consider storing the
    // multiplier and row address to avoid any extra operations in the main loop:
    //      const double aji = A[j][i];
    //      double *aj = A[j];
    //      double *ai = A[i];
    for (size_t k = i+1; k < m; ++k) {
        A(j,k) -= A(i,k) * A(j,i);      // update submatrix
    }
    B(j) -= B(i) * A(j,i);
}

///
/// @brief Gauss solver of a linear system of equations (serial). Compute matrix
/// left division X = A\B = inv(A)*B. ie solve the linear algebra problem A*X = B
/// using Gauss Elimination with Partial Pivoting:
///     neq         number of equations
///     A(neq,neq)  coefficient matrix
///     B(neq)      nonhomogeneous vector
///     X(neq)      solution vector
///
template<typename T = double>
void GaussSolve(
    Matrix<T> & __restrict__ A,
    Vector<T> & __restrict__ B,
    Vector<T> & __restrict__ X)
{
    static_assert(std::is_floating_point<T>::value, "non floating point type");
    assert(A.n1() == A.n2() &&
           B.n1() == A.n2() &&
           X.n1() == A.n2() &&
           "invalid dimensions");

    // Number of equations
    const size_t neq = A.n1();

    // Forward Elimination with Partial Pivoting. Go over the neq-1 columns of
    // matrix A. For each column-i, find the pivot row in the submatrix
    //      A(i:neq-1, i:neq-1).
    //
    // The pivot row-k gives the index of the element in column-i
    //      |A(k,i)| = max |A(j,i)|, for (i <= j < neq).
    //
    // Store row-k by swaping it with row-i. Using the pivot row, update the
    // submatrix
    //      A(i:neq-1, i:neq-1)
    // by forward elimination of the elements A(j,i) below the diagonal, i.e.
    // (i < j < neq).
    //
    // Forward elimination step goes as follows. Store all row multipliers in the
    // zeroed entries (j,i) below the diagonal element A(i,i) before updating,
    //      A(j,i) = A(j,i) / A(i,i).
    //
    // After computing the multipliers update the rest of the matrix A, i.e. the
    // submatrix
    //      A(i+1:neq-1, i+1:neq-1)
    // by subtracting row A(i, i+1:neq-1) from rows A(k, i+1:neq-1).
    //
    for (size_t i = 0; i < neq-1; ++i) {
        // Find pivot row-k below diagonal element(i,i). For each row j below i,
        // if |A(j,i)| > |A(k,i)|, store k = j.
        size_t k = i;
        for (size_t j = i+1; j < neq; ++j) {
            if (abs(A(k,i)) < abs(A(j,i))) {
                k = j;
            }
        }

        // Check matrix condition
        assert(std::fabs(A(k,i)) > (T) 0 &&
            "ill conditioned matrix: division by zero");

        // Swap row-i with pivot row-k if row-k is below row-i
        if (i < k) {
            for (size_t j = i; j < neq; ++j) {
                std::swap(A(i,j), A(k,j));
            }
            std::swap(B(i), B(k));
        }

        // Store all row multipliers in the zeroed entries below diagonal element
        // A(i,i) before updating. Update submatrix A(i+1:neq-1, i+1:neq-1) by
        // subtracting row-i A(i, i+1:neq-1) from rows-k A(k, i+1:neq-1).
        for (size_t j = i+1; j < neq; ++j) {
            A(j,i) /= A(i,i);     // store multiplier A(j,i) / A(i,i)
        }

        for (size_t j = i+1; j < neq; ++j) {
            GaussSubtractRow<T>(j, i, neq, A, B);  // update submatrix
        }
    }

    // Backward Substitution. Normally, it would iterate (i = neq-1; i >= 0; i--).
    // Since i is an unsigned int type, when it reaches -1 at the last iteration,
    // it folds to the maximum value. The condition (i < neq) is therefore used
    // instead of(i >= 0).
    for (size_t i = neq-1; i < neq; i--) {
        T sum = (T) 0;
        for (size_t j = i+1; j < neq; ++j) {
            sum += A(i,j) * X(j);
        }
        X(i) = (B(i) - sum) / A(i,i);
    }
}

///
/// @brief Gauss solver of a linear system of equations(parallel). Compute matrix
/// left division X = A\B = inv(A)*B. i.e. solve the linear algebra problem
/// A*X = B using Gauss Elimination with Partial Pivoting:
///     neq          number of equations
///     A(neq,neq)   coefficient matrix
///     B(neq)       nonhomogeneous vector
///     X(neq)       solution vector
///
template<typename T = double>
void GaussSolveOmp(
    Matrix<T> & __restrict__ A,
    Vector<T> & __restrict__ B,
    Vector<T> & __restrict__ X)
{
    static_assert(std::is_floating_point<T>::value, "non floating point type");
    assert(A.n1() == A.n2() &&
           B.n1() == A.n2() &&
           X.n1() == A.n2() &&
           "invalid dimensions");

    const size_t neq = A.n1();

    // Forward Elimination with Partial Pivoting. Go over the neq-1 columns of
    // matrix A. For each column-i, find the pivot row in the submatrix
    //      A(i:neq-1, i:neq-1).
    // The pivot row-k gives the index of the element in column-i
    //      |A(k,i)| = max |A(j,i)|, for (i <= j < neq).
    // Store row-k by swaping it with row-i. Using the pivot row, update the
    // submatrix
    //      A(i:neq-1, i:neq-1)
    // by forward elimination of the elements A(j,i) below the diagonal, i.e.
    // (i < j < neq).
    //
    // Forward elimination step goes as follows:
    // Store all row multipliers in the zeroed entries (j,i) below the diagonal
    // element A(i,i) before updating,
    //      A(j,i) = A(j,i) / A(i,i).
    //
    // After computing the multipliers update the rest of the matrix A, i.e. the
    // submatrix
    //      A(i+1:neq-1, i+1:neq-1)
    //
    // by subtracting row A(i, i+1:neq-1) from rows A(k, i+1:neq-1).
    //
    // @note The forward elimination is sequential in row index i, therefore it
    // cannot be parallelized using a for loop construct. However, the internal
    // for loop operations (pivot finding, row swapping and row elimination)
    // can be parallized.
    // To avoid the overhead of constructing a team of threads neq times, the for
    // loop
    //      for (i=0; i < neq-1; ++i)
    //
    // is replaced by an equivalent while loop
    //      i=0;
    //      while(i < neq-1)
    //          ++i;
    //
    // where index i is a shared variable incremented by a single thread at the
    // end of the while loop step.
    size_t i = 0, k = 0;    // Initialize equation i-index and pivot k-index
    #pragma omp parallel default(none) shared(neq, i, k, A, B)
    {
        while (i < neq-1)
        {
            // Find pivot row-k below diagonal element A(i,i). For each row-j
            // below row-i, if
            //  |A[j][i]| > |A[k][i]|, then set k = j.
            // Each thread stores its private row-k index, kk. At the end,
            // k = argmax(|A[kk][i]|).
            size_t kk = i;
            #pragma omp for schedule(dynamic)
            for (size_t j = i+1; j < neq; ++j) {
                if (std::fabs(A(kk,i)) <
                    std::fabs(A(j,i))) {
                    kk = j;
                }
            }

            #pragma omp critical
            if (std::fabs(A(k,i)) < std::fabs(A(kk,i))) {
                k = kk;
            }
            // omp critical

            // Check matrix condition
            #pragma omp single
            assert(std::fabs(A(k,i)) > (T) 0 &&
                   "ill conditioned matrix: division by zero");
            // omp single


            // Swap row-i with pivot row-k if row-k is below row-i
            if (i < k) {
                #pragma omp for schedule(dynamic)
                for (size_t j = i; j < neq; ++j) {
                    std::swap(A(i,j), A(k,j));
                }

                #pragma omp single
                std::swap(B(i), B(k));
                // omp single
            }

            // Store all row multipliers in the zeroed entries below diagonal
            // element A(i,i) before updating. Update submatrix
            // A(i+1:neq-1, i+1:neq-1) by subtracting row-i A(i, i+1:neq-1)
            // from rows-k A(k, i+1:neq-1)
            #pragma omp for schedule(dynamic)
            for (size_t j = i+1; j < neq; ++j) {
                A(j,i) /= A(i,i);       // store multiplier A(j,i) / A(i,i)
            }

            #pragma omp for schedule(dynamic)
            for (size_t j = i+1; j < neq; ++j) {
                GaussSubtractRow<T>(j, i, neq, A, B); // update submatrix
            }

            // Increment equation i-index and initialise pivot row k-index.
            #pragma omp single
            k = ++i;
            // omp single
        }
    } // omp parallel

    // Backward Substitution. Normally, it would iterate (i = neq-1; i >= 0; i--)
    // but since i is an unsigned int type, when it reaches -1 at the last
    // iteration, it unfolds to the maximum positive value.
    // Assuming neq is smaller than max unsigned int, condition(i < neq) is used
    // instead of(i >= 0).
    //
    // Similar to the forward elimination step, the backward substitution step is
    // sequential in index i, i.e. X[n-1] needs to be computed before X[n-2], etc.
    //
    // The backward loop
    //      for (i = neq-1; i < neq; i--) {}
    //
    // is replaced by
    //      i = neq-1;
    //      while (i < neq)
    //          i--;
    // with the index decrement operation is computed by a single one thread.
    T sum = (T) 0;
    i = neq-1;
    #pragma omp parallel default(none) shared(neq, i, sum, A, B, X)
    {
        while (i < neq)
        {
            T thr_sum = 0.0;
            #pragma omp for schedule(dynamic)
            for (size_t j = i+1; j < neq; ++j) {
                thr_sum += A(i,j) * X(j);
            }

            #pragma omp critical
            sum += thr_sum;
            // omp critical

            // The flush directive is NOT implied in the following cases:
            //      at entry to for,
            //      at entry to or exit from master,
            //      at entry to sections,
            //      at entry to single,
            //      at exit of for, single or sections, if nowait clause is used.
            // The implicit synchronization via an implied barrier directive is
            // performed ONLY at exit from for, single and sections directives,
            // unless nowait clause is used.
            // We need a barrier here to ensure all threads add their value to the
            // sum variable before using it to compute the X[i] variable.
            // The sum variable also needs to be flushed at the entry of single
            // since this is not implicitly flushed as described above.
            #pragma omp barrier
            #pragma omp single
            {
                #pragma omp flush(sum)
                X(i) = (B(i) - sum) / A(i,i);
                sum = 0.0;
                --i;
            } // omp single
        }
    } // omp parallel
}

#endif // NUMERIC_LINALG_GAUSS_H_
