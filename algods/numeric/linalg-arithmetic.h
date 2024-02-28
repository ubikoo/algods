//
// linalg-arithmetic.h
//
// Copyright (c) 2020 Carlos Braga
// This program is free software; you can redistribute it and/or modify it
// under the terms of the MIT License. See accompanying LICENSE.md or
// https://opensource.org/licenses/MIT.
//

#ifndef NUMERIC_LINALG_ARITHMETIC_H_
#define NUMERIC_LINALG_ARITHMETIC_H_

#include <omp.h>
#include "algods/numeric/array.h"

///
/// @brief Matrix-Vector multiplication of matrix A by vector B, A * B = C,
///     C(i) = A(i,k) * B(k),
/// where C(i) is the dot product of the row-i of matrix A and vector B.
/// m is the number of rows in matrix A and in vector C
/// n is the number of columns in matrix A and rows in vector B
/// A is a pointer to matrix A
/// B is a pointer to vector B
/// C is a pointer to vector C
/// T must be a floating point data type.
///
template<typename T>
void MatmulVectorRow(
    size_t i,
    Matrix<T> & __restrict__ A,
    Vector<T> & __restrict__ B,
    Vector<T> & __restrict__ C)
{
    T sum = (T) 0;
    for (size_t j = 0; j < A.n2(); ++j) {
        sum += A(i,j) * B(j);
    }
    C(i) = sum;
}

template<typename T, bool IsPar = true>
void MatmulVector(
    Matrix<T> & __restrict__ A,
    Vector<T> & __restrict__ B,
    Vector<T> & __restrict__ C)
{
    static_assert(std::is_floating_point<T>::value, "non floating point type");
    assert(C.n1() == A.n1() && B.n1() == A.n2());

    // C(i) = A(i,j) * B(j)
    #pragma omp parallel for if(IsPar) default(none) shared(A, B, C) \
        schedule(static)
    for (size_t i = 0; i < C.n1(); ++i) {
        MatmulVectorRow<T>(i, A, B, C);
    }
}

///
/// @brief Matrix-Matrix multiplication of matrix A by matrix B, A * B = C:
///     C(i,j) = A(i,k)*B(k,j),
/// where C(i,j) is the dot product of the row i of matrix A and the vector
/// column j of matrix B.
/// m is the number of rows in matrix A and in matrix C
/// n is the number of columns in matrix A and rows in matrix B
/// p is the number of columns in matrix B and in matrix C
/// A is a pointer to matrix A
/// B is a pointer to matrix B
/// C is a pointer to result C
///
template<typename T>
void MatmulMatrixRow(
    size_t i,
    Matrix<T> & __restrict__ A,
    Matrix<T> & __restrict__ B,
    Matrix<T> & __restrict__ C)
{
    for (size_t j = 0; j < C.n2(); ++j) {
        T sum = (T) 0;
        for (size_t k = 0; k < A.n2(); ++k) {
            sum += A(i,k) * B(k,j);
        }
        C(i,j) = sum;
    }
}

template<typename T, bool IsPar = true>
void MatmulMatrix(
    Matrix<T> & __restrict__ A,
    Matrix<T> & __restrict__ B,
    Matrix<T> & __restrict__ C)
{
    static_assert(std::is_floating_point<T>::value, "non floating point type");
    assert(C.n1() == A.n1() && C.n2() == B.n2() && A.n2() == B.n1());

    // C(i,j) = A(i,k) * B(k,j)
    #pragma omp parallel for if(IsPar) default(none) shared(A, B, C) \
        schedule(static)
    for (size_t i = 0; i < C.n1(); ++i) {
        MatmulMatrixRow<T>(i, A, B, C);
    }
}

///
/// @brief Zero the vector elements.
///
template<typename T, bool IsPar = true>
void ZeroVector(Vector<T> & __restrict__ vec)
{
    static_assert(std::is_floating_point<T>::value, "non floating point type");
    assert(vec.n1() > 0);

    constexpr T zero = (T) 0;
    #pragma omp parallel for if(IsPar) default(none) shared(vec) \
        schedule(static)
    for (size_t i = 0; i < vec.n1(); ++i) {
        vec(i) = zero;
    }
}

///
/// @brief Zero the matrix elements.
///
template<typename T, bool IsPar = true>
void ZeroMatrix(Matrix<T> & __restrict__ mat)
{
    static_assert(std::is_floating_point<T>::value, "non floating point type");
    assert(mat.n1() > 0 && mat.n2() > 0);

    constexpr T zero = (T) 0;
    #pragma omp parallel for if(IsPar) default(none) shared(mat) \
        schedule(static)
    for (size_t i = 0; i < mat.n1(); ++i) {
        for (size_t j = 0; j < mat.n2(); ++j) {
            mat(i,j) = zero;
        }
    }
}

///
/// @brief Zero the tensor elements.
///
template<typename T, bool IsPar = true>
void ZeroTensor(Tensor<T> & __restrict__ tns)
{
    static_assert(std::is_floating_point<T>::value, "non floating point type");
    assert(tns.n1() > 0 && tns.n2() > 0 && tns.n3() > 0);

    constexpr T zero = (T) 0;
    #pragma omp parallel for if(IsPar) default(none) shared(tns) \
        schedule(static)
    for (size_t i = 0; i < tns.n1(); ++i) {
        for (size_t j = 0; j < tns.n2(); ++j) {
            for (size_t k = 0; k < tns.n3(); ++k) {
                tns(i,j,k) = zero;
            }
        }
    }
}

///
/// @brief Copy the vector elements.
///
template<typename T, bool IsPar = true>
void CopyVector(Vector<T> & __restrict__ src, Vector<T> & __restrict__ dst)
{
    static_assert(std::is_floating_point<T>::value, "non floating point type");
    assert(src.n1() > 0  && src.n1() == dst.n1());

    #pragma omp parallel for if(IsPar) default(none) shared(src, dst) \
        schedule(static)
    for (size_t i = 0; i < src.n1(); ++i) {
        dst(i) = src(i);
    }
}

///
/// @brief Copy the matrix elements.
///
template<typename T, bool IsPar = true>
void CopyMatrix(Matrix<T> & __restrict__ src, Matrix<T> & __restrict__ dst)
{
    static_assert(std::is_floating_point<T>::value, "non floating point type");
    assert(src.n1() > 0         &&
           src.n2() > 0         &&
           src.n1() == dst.n1() &&
           src.n2() == dst.n2());

    #pragma omp parallel for if(IsPar) default(none) shared(src, dst) \
        schedule(static)
    for (size_t i = 0; i < src.n1(); ++i) {
        for (size_t j = 0; j < src.n2(); ++j) {
            dst(i,j) = src(i,j);
        }
    }
}

///
/// @brief Copy the tensor elements.
///
template<typename T, bool IsPar = true>
void CopyTensor(Tensor<T> & __restrict__ src, Tensor<T> & __restrict__ dst)
{
    static_assert(std::is_floating_point<T>::value, "non floating point type");
    assert(src.n1() > 0         &&
           src.n2() > 0         &&
           src.n3() > 0         &&
           src.n1() == dst.n1() &&
           src.n2() == dst.n2() &&
           src.n3() == dst.n3());

    #pragma omp parallel for if(IsPar) default(none) shared(src, dst) \
        schedule(static)
    for (size_t i = 0; i < src.n1(); ++i) {
        for (size_t j = 0; j < src.n2(); ++j) {
            for (size_t k = 0; k < src.n3(); ++k) {
                dst(i,j,k) = src(i,j,k);
            }
        }
    }
}

///
/// @brief Compute diagonal matrix from vector elements.
///
template<typename T, bool IsPar = true>
void DiagMatrix(Vector<T> & __restrict__ vec, Matrix<T> & __restrict__ mat)
{
    static_assert(std::is_floating_point<T>::value, "non floating point type");
    assert(vec.n1() > 0         &&
           mat.n1() > 0         &&
           mat.n2() > 0         &&
           vec.n1() == mat.n1() &&
           vec.n1() == mat.n2());

    ZeroMatrix<T,IsPar>(mat);

    #pragma omp parallel for if(IsPar) default(none) shared(mat, vec) \
        schedule(static)
    for (size_t i = 0; i < mat.n1(); ++i) {
        mat(i,i) = vec(i);
    }
}

///
/// @brief Compute vector from matrix diagonal.
///
template<typename T, bool IsPar = true>
void DiagVector(Matrix<T> & __restrict__ mat, Vector<T> & __restrict__ vec)
{
    static_assert(std::is_floating_point<T>::value, "non floating point type");
    assert(vec.n1() > 0         &&
           mat.n1() > 0         &&
           mat.n2() > 0         &&
           vec.n1() == mat.n1() &&
           vec.n1() == mat.n2());

    ZeroVector<T,IsPar>(vec);

    #pragma omp parallel for if(IsPar) default(none) shared(mat, vec) \
        schedule(static)
    for (size_t i = 0; i < vec.n1(); ++i) {
        vec(i) = mat(i,i);
    }
}

///
/// @brief Transpose matrix.
///
template<typename T, bool IsPar = true>
void TransposeMatrix(
    Matrix<T> & __restrict__ mat,
    Matrix<T> & __restrict__ trans)
{
    static_assert(std::is_floating_point<T>::value, "non floating point type");
    assert(mat.n1() > 0             &&
           mat.n2() > 0             &&
           trans.n1() == mat.n2()   &&
           trans.n2() == mat.n1());

    #pragma omp parallel for if(IsPar) default(none) shared(mat, trans) \
        schedule(dynamic)
    for (size_t i = 0; i < mat.n1(); ++i) {
        for (size_t j = i; j < mat.n2(); ++j) {
            trans(i,j) = mat(j,i);
            trans(j,i) = mat(i,j);
        }
    }
}

///
/// @brief Identity matrix.
///
template<typename T, bool IsPar = true>
void IdentityMatrix(Matrix<T> & __restrict__ mat)
{
    static_assert(std::is_floating_point<T>::value, "non floating point type");
    assert(mat.n1() > 0         &&
           mat.n2() > 0         &&
           mat.n2() == mat.n1());

    ZeroMatrix<T,IsPar>(mat);

    constexpr T one = (T) 1;
    #pragma omp parallel for if(IsPar) default(none) shared(mat) \
        schedule(static)
    for (size_t i = 0; i < mat.n1(); ++i) {
        mat(i,i) = one;
    }
}

#endif // NUMERIC_LINALG_ARITHMETIC_H_
