//
// test-kdtree.cpp
//
// Copyright (c) 2020 Carlos Braga
// This program is free software; you can redistribute it and/or modify it
// under the terms of the MIT License. See accompanying LICENSE.md or
// https://opensource.org/licenses/MIT.
//

#include "catch2/catch.hpp"
#include "algods/searching/kdtree-static.h"
#include <random>
#include <array>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <chrono>

/// ---- KdTreeStatic time sphere_query ----------------------------------------
/// Return the elapsed time in miliseconds of a sphere_query
///
template<typename Type, std::size_t Dim>
double test_kdtree_time_sphere_query(
    KdTreeStatic<Type,Dim> &kdtree,
    const std::array<Type,Dim> &point,
    Type radius)
{
    // Start time
    auto tic = std::chrono::high_resolution_clock::now();

    // Run
    auto nearest = kdtree.sphere_query(point, radius);

    // End time
    auto toc = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double,std::ratio<1,1000>> msec = toc-tic;

    return msec.count();
}


template<typename Type, std::size_t Dim>
void test_kdtree_random_sphere_queries(
    KdTreeStatic<Type,Dim> &kdtree,
    Type point_lo,
    Type point_hi,
    Type radius_lo,
    Type radius_hi,
    size_t numqueries)
{
    std::random_device rd;      // rng device
    std::mt19937 rng(rd());     // rng engine seeded with rd
    std::uniform_real_distribution<Type> point_dist(point_lo, point_hi);
    std::uniform_real_distribution<Type> radius_dist(radius_lo, radius_hi);

    double query_time = 0.0;
    size_t query_count = 0;
    for (size_t i = 0; i < numqueries; ++i) {
        std::array<Type,Dim> point = {point_dist(rng), point_dist(rng)};
        Type radius = radius_dist(rng);
        double time = test_kdtree_time_sphere_query<Type,Dim>(
            kdtree, point, radius);

        if (time > 0.0) {
            query_time += time;
            query_count++;
        }
    }

    double avg_query_time = query_time / static_cast<double>(query_count);
    std::cout << "sphere_query avg_query_time " << avg_query_time << " msec\n";
}


/// ---- KdTreeStatic time range_query -----------------------------------------
/// Return the elapsed time in miliseconds of a sphere_query
///
template<typename Type, std::size_t Dim>
double test_kdtree_time_range_query(
    KdTreeStatic<Type,Dim> &kdtree,
    const std::array<Type,Dim> &qlo,
    const std::array<Type,Dim> &qhi)
{
    // Start time
    auto tic = std::chrono::high_resolution_clock::now();

    // Run
    auto nearest = kdtree.range_query(qlo, qhi);

    // End time
    auto toc = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double,std::ratio<1,1000>> msec = toc-tic;

    return msec.count();
}


template<typename Type, std::size_t Dim>
void test_kdtree_random_range_queries(
    KdTreeStatic<Type,Dim> &kdtree,
    Type range_lo,
    Type range_hi,
    size_t numqueries)
{
    std::random_device rd;      // rng device
    std::mt19937 rng(rd());     // rng engine seeded with rd
    std::uniform_real_distribution<Type> range_dist(range_lo, range_hi);

    double query_time = 0.0;
    size_t query_count = 0;
    for (size_t i = 0; i < numqueries; ++i) {
        std::array<Type,Dim> qlo = {range_dist(rng), range_dist(rng)};
        std::array<Type,Dim> qhi = {range_dist(rng), range_dist(rng)};
        double time = test_kdtree_time_range_query<Type,Dim>(kdtree, qlo, qhi);

        if (time > 0.0) {
            query_time += time;
            query_count++;
        }
    }

    double avg_query_time = query_time / static_cast<double>(query_count);
    std::cout << "range_query avg_query_time " << avg_query_time << " msec\n";
}


/// ---- Test KSTree data structure --------------------------------------------
///
TEST_CASE("k-d search tree") {
    //
    // Generate tree
    //
    std::random_device rd;      // rng device
    std::mt19937 rng(rd());     // rng engine seeded with rd
    std::uniform_real_distribution<double> dist(0.0, 1.0);

    const size_t numpoints = 100000;
    std::vector<std::array<double,2>> points;
    for (size_t i = 0; i < numpoints; ++i) {
        if (i == 0) {
            std::array<double,2> q = {0.0, 0.0};
            points.push_back(q);
        } else if (i == 1) {
            std::array<double,2> q = {1.0, 1.0};
            points.push_back(q);
        } else {
            std::array<double,2> q = {dist(rng), dist(rng)};
            points.push_back(q);
        }
    }

    // ----------------------------------------------------------------------
    SECTION("Test API") {
        KdTreeStatic<double, 2> kdtree;

        kdtree.insert(points[0]);
        kdtree.insert(points[1]);
        size_t count = 2;
        for (size_t i = 2; i < points.size(); ++i) {
            kdtree.insert(points[i]);
            REQUIRE(kdtree.contains(points[i]));
            REQUIRE(kdtree.size() == ++count);
            REQUIRE(!kdtree.is_empty());

            auto qmin0 = kdtree.min(0);
            REQUIRE(qmin0[0] == 0.0);
            REQUIRE(qmin0[1] == 0.0);

            auto qmin1 = kdtree.min(1);
            REQUIRE(qmin1[0] == 0.0);
            REQUIRE(qmin1[1] == 0.0);

            auto qmax0 = kdtree.max(0);
            REQUIRE(qmax0[0] == 1.0);
            REQUIRE(qmax0[1] == 1.0);

            auto qmax1 = kdtree.max(1);
            REQUIRE(qmax1[0] == 1.0);
            REQUIRE(qmax1[1] == 1.0);
        }
        kdtree.clear();
        REQUIRE(kdtree.is_empty());
    }

    // ----------------------------------------------------------------------
    SECTION("Test knn_query") {
        KdTreeStatic<double, 2> kdtree(points);

        std::array<double,2> q = {1.5, 1.5};
        auto nearest = kdtree.knn_query(q, 5);
        for (auto p : nearest) {
            std::cout << std::setw(8) << std::fixed << p[0] << "\t"
                      << std::setw(8) << std::fixed << p[1] << std::endl;
        }
    }

    // ----------------------------------------------------------------------
    SECTION("Test sphere_query") {
        KdTreeStatic<double, 2> kdtree(points);

        std::array<double,2> q = {1.5, 1.5};
        auto nearest = kdtree.sphere_query(q, 0.71 /* sqrt(2*0.5*0.5) */);
        for (auto p : nearest) {
            std::cout << std::setw(8) << std::fixed << p[0] << "\t"
                      << std::setw(8) << std::fixed << p[1] << std::endl;
        }

        KdTreeStatic<double, 2> kdtree_temp(kdtree);
        KdTreeStatic<double, 2> kdtree_copy = kdtree_temp;
        test_kdtree_random_sphere_queries<double,2>(
            kdtree_copy, -1.5, 1.5, 0.5, 1.0, 10000);
    }

    // ----------------------------------------------------------------------
    SECTION("Test range_query") {
        KdTreeStatic<double, 2> kdtree(points);

        std::array<double,2> qlo = {0.99, 0.99};
        std::array<double,2> qhi = {1.01, 1.01};
        auto nearest = kdtree.range_query(qlo, qhi);
        for (auto p : nearest) {
            std::cout << std::setw(8) << std::fixed << p[0] << "\t"
                      << std::setw(8) << std::fixed << p[1] << std::endl;
        }

        KdTreeStatic<double, 2> kdtree_temp(kdtree);
        KdTreeStatic<double, 2> kdtree_copy = kdtree_temp;
        test_kdtree_random_range_queries<double,2>(
            kdtree_copy, -1.5, 1.5, 10000);
    }
}
