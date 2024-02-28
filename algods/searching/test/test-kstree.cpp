//
// test-kstree.cpp
//
// Copyright (c) 2020 Carlos Braga
// This program is free software; you can redistribute it and/or modify it
// under the terms of the MIT License. See accompanying LICENSE.md or
// https://opensource.org/licenses/MIT.
//

#include "catch2/catch.hpp"
#include "algods/searching/kdtree-dynamic.h"
#include <array>
#include <chrono>
#include <random>
#include <limits>
#include <iomanip>
#include <iostream>

/// ---- KdTreeDynamic time sphere_query ---------------------------------------
/// Return the elapsed time in miliseconds of a sphere_query
///
template<typename Type, std::size_t Dim>
double test_kstree_time_sphere_query(
    KdTreeDynamic<Type,Dim> &kstree,
    const std::array<Type,Dim> &point,
    Type radius)
{
    // Start time
    auto tic = std::chrono::high_resolution_clock::now();

    // Run
    auto nearest = kstree.sphere_query(point, radius);

    // End time
    auto toc = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double,std::ratio<1,1000>> msec = toc-tic;

    return nearest.size() > 0 ? msec.count() : 0.0;
}


template<typename Type, std::size_t Dim>
void test_kstree_random_sphere_queries(
    KdTreeDynamic<Type,Dim> &kstree,
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
        double time = test_kstree_time_sphere_query<Type,Dim>(
            kstree, point, radius);

        if (time > 0.0) {
            query_time += time;
            query_count++;
        }
    }

    double avg_query_time = query_time / static_cast<double>(query_count);
    std::cout << "sphere_query avg_query_time " << avg_query_time << " msec\n";
}


/// ---- KdTreeDynamic time range_query ----------------------------------------
/// Return the elapsed time in miliseconds of a sphere_query
///
template<typename Type, std::size_t Dim>
double test_kstree_time_range_query(
    KdTreeDynamic<Type,Dim> &kstree,
    const std::array<Type,Dim> &qlo,
    const std::array<Type,Dim> &qhi)
{
    // Start time
    auto tic = std::chrono::high_resolution_clock::now();

    // Run
    auto nearest = kstree.range_query(qlo, qhi);

    // End time
    auto toc = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double,std::ratio<1,1000>> msec = toc-tic;

    return nearest.size() > 0 ? msec.count() : 0.0;
}


template<typename Type, std::size_t Dim>
void test_kstree_random_range_queries(
    KdTreeDynamic<Type,Dim> &kstree,
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
        double time = test_kstree_time_range_query<Type,Dim>(kstree, qlo, qhi);

        if (time > 0.0) {
            query_time += time;
            query_count++;
        }
    }

    double avg_query_time = query_time / static_cast<double>(query_count);
    std::cout << "range_query avg_query_time " << avg_query_time << " msec\n";
}


/// ---- Test KdTreeDynamic data structure -------------------------------------
///
TEST_CASE("k-d Binary search tree") {
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
        KdTreeDynamic<double, 2> kstree;

        kstree.insert(points[0]);
        kstree.insert(points[1]);
        size_t count = 2;
        for (size_t i = 2; i < points.size(); ++i) {
            kstree.insert(points[i]);
            REQUIRE(kstree.contains(points[i]));
            REQUIRE(kstree.size() == ++count);
            REQUIRE(!kstree.is_empty());

            auto qmin0 = kstree.min(0);
            REQUIRE(qmin0[0] == 0.0);
            REQUIRE(qmin0[1] == 0.0);

            auto qmin1 = kstree.min(1);
            REQUIRE(qmin1[0] == 0.0);
            REQUIRE(qmin1[1] == 0.0);

            auto qmax0 = kstree.max(0);
            REQUIRE(qmax0[0] == 1.0);
            REQUIRE(qmax0[1] == 1.0);

            auto qmax1 = kstree.max(1);
            REQUIRE(qmax1[0] == 1.0);
            REQUIRE(qmax1[1] == 1.0);
        }
        kstree.clear();
        REQUIRE(kstree.is_empty());
    }


    // ----------------------------------------------------------------------
    SECTION("Test x-min") {
        KdTreeDynamic<double, 2> kstree(points);

        std::array<double,2> qold = {std::numeric_limits<double>::lowest(),
                                     std::numeric_limits<double>::lowest()};
        while (!kstree.is_empty()) {
            auto q = kstree.min(0);
            REQUIRE(q[0] > qold[0]);
            qold = q;
            kstree.erase(q);
            REQUIRE(!kstree.contains(q));
        }
        REQUIRE(kstree.is_empty());
    }


    // ----------------------------------------------------------------------
    SECTION("Test x-max") {
        KdTreeDynamic<double, 2> kstree(points);

        std::array<double,2> qold = {std::numeric_limits<double>::max(),
                                     std::numeric_limits<double>::max()};
        while (!kstree.is_empty()) {
            auto q = kstree.max(0);
            REQUIRE(q[0] < qold[0]);
            qold = q;
            kstree.erase(q);
            REQUIRE(!kstree.contains(q));
        }
        REQUIRE(kstree.is_empty());
    }


    // ----------------------------------------------------------------------
    SECTION("Test y-min") {
        KdTreeDynamic<double, 2> kstree(points);

        std::array<double,2> qold = {std::numeric_limits<double>::lowest(),
                                     std::numeric_limits<double>::lowest()};
        while (!kstree.is_empty()) {
            auto q = kstree.min(1);
            REQUIRE(q[1] > qold[1]);
            qold = q;
            kstree.erase(q);
            REQUIRE(!kstree.contains(q));
        }
        REQUIRE(kstree.is_empty());
    }


    // ----------------------------------------------------------------------
    SECTION("Test y-max") {
        KdTreeDynamic<double, 2> kstree(points);

        std::array<double,2> qold = {std::numeric_limits<double>::max(),
                                     std::numeric_limits<double>::max()};
        while (!kstree.is_empty()) {
            auto q = kstree.max(1);
            REQUIRE(q[1] < qold[1]);
            qold = q;
            kstree.erase(q);
            REQUIRE(!kstree.contains(q));
        }
        REQUIRE(kstree.is_empty());
    }

    // ----------------------------------------------------------------------
    SECTION("Test knn_query") {
        KdTreeDynamic<double, 2> kstree(points);

        std::array<double,2> q = {1.5, 1.5};
        auto nearest = kstree.knn_query(q, 5);
        for (auto p : nearest) {
            std::cout << std::setw(8) << std::fixed << p[0] << "\t"
                      << std::setw(8) << std::fixed << p[1] << std::endl;
        }
    }

    // ----------------------------------------------------------------------
    SECTION("Test sphere_query") {
        KdTreeDynamic<double, 2> kstree(points);

        std::array<double,2> q = {1.5, 1.5};
        auto nearest = kstree.sphere_query(q, 0.71 /* sqrt(2*0.5*0.5) */);
        for (auto p : nearest) {
            std::cout << std::setw(8) << std::fixed << p[0] << "\t"
                      << std::setw(8) << std::fixed << p[1] << std::endl;
        }

        KdTreeDynamic<double, 2> kstree_temp(kstree);
        KdTreeDynamic<double, 2> kstree_copy = kstree_temp;
        test_kstree_random_sphere_queries<double,2>(
            kstree_copy, -1.5, 1.5, 0.5, 1.0, 10000);
    }

    // ----------------------------------------------------------------------
    SECTION("Test range_query") {
        KdTreeDynamic<double, 2> kstree(points);

        std::array<double,2> qlo = {0.99, 0.99};
        std::array<double,2> qhi = {1.01, 1.01};
        auto nearest = kstree.range_query(qlo, qhi);
        for (auto p : nearest) {
            std::cout << std::setw(8) << std::fixed << p[0] << "\t"
                      << std::setw(8) << std::fixed << p[1] << std::endl;
        }

        KdTreeDynamic<double, 2> kstree_temp(kstree);
        KdTreeDynamic<double, 2> kstree_copy = kstree_temp;
        test_kstree_random_range_queries<double,2>(
            kstree_copy, -1.5, 1.5, 10000);
    }
}
