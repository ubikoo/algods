//
// test-ode-kepler.cpp
//
// Copyright (c) 2020 Carlos Braga
// This program is free software; you can redistribute it and/or modify it
// under the terms of the MIT License. See accompanying LICENSE.md or
// https://opensource.org/licenses/MIT.
//

#include "catch2/catch.hpp"
#include <cmath>
#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
#include "algods/numeric/array.h"
#include "algods/numeric/ode.h"

/// ---- Kepler model ----------------------------------------------------------
struct Kepler {
    double delta_ = 0.0;

    Kepler(double delta) : delta_(delta) {}
    ~Kepler() = default;

    void operator()(Vector<double> & __restrict__ x,
                    Vector<double> & __restrict__ dxdt) {
        double r = std::sqrt(x(0)*x(0) + x(1)*x(1));
        double inv_r3 = 1.0 / (r*r*r);
        double inv_r5 = 1.0 / (r*r*r*r*r);

        dxdt(0) = x(2);
        dxdt(1) = x(3);
        dxdt(2) = -(inv_r3 + 0.5*3.0*delta_*inv_r5)*x(0);
        dxdt(3) = -(inv_r3 + 0.5*3.0*delta_*inv_r5)*x(1);
    }

    double Energy(Vector<double> & __restrict__ x) {
        double r = std::sqrt(x(0)*x(0) + x(1)*x(1));

        double e_kin = 0.5*(x(2)*x(2) + x(3)*x(3));
        double e_pot = -(1.0/r) - (0.5*delta_/r*r*r);
        return e_kin + e_pot;
    }
};

/// ---- Test Ode Gauss Integrator ---------------------------------------------
void TestOdeKeplerGauss(double dt, size_t nsteps, size_t smpfreq)
{
    // Model parameters
    const double delta = 0.01;
    Kepler kepler(delta);

    // System parameters
    const size_t ndim = 4;
    const double e = 0.6;
    Vector<double> x(ndim);          // x(t)
    x(0) = 1.0 - e;
    x(1) = 0.0;
    x(2) = 0.0;
    x(3) = std::sqrt((1.0 + e)/(1.0 - e));
    Vector<double> dxdt(ndim);       // dx(t) / dt

    // Integrate system
    double t = 0.0;
    OdeGauss<double, false> gauss(ndim);
    gauss.Init(kepler, x, dxdt, dt);

    std::ofstream ofs;
    ofs.open(std::string("/tmp/out.kepler_gauss"));
    ofs << std::setprecision(std::numeric_limits<double>::max_digits10);
    ofs << std::scientific;
    ofs << "#  t  x(0)  x(1)  x(2)  x(3)  e\n";
    for (size_t i = 0; i < nsteps; ++i) {
        // sample
        if (i%smpfreq == 0) {
            ofs << t << " "
                << x(0) << " "
                << x(1) << " "
                << x(2) << " "
                << x(3) << " "
                << kepler.Energy(x) << "\n";
        }
        // x = x + dxdt * dt
        gauss.Step(kepler, x, dxdt, dt);
        t += dt;
    }
    ofs.close();
}

/// ---- Test Ode Runge-Kutta Integrator ---------------------------------------
void TestOdeKeplerKutta(double dt, size_t nsteps, size_t smpfreq)
{
    // Model parameters
    const double delta = 0.01;
    Kepler kepler(delta);

    // System parameters
    const size_t ndim = 4;
    const double e = 0.6;
    Vector<double> x(ndim);          // x(t)
    x(0) = 1.0 - e;
    x(1) = 0.0;
    x(2) = 0.0;
    x(3) = std::sqrt((1.0 + e)/(1.0 - e));
    Vector<double> dxdt(ndim);       // dx(t) / dt

    // Integrate system
    double t = 0.0;
    OdeKutta<double, false> kutta(ndim);
    kutta.Init(kepler, x, dxdt, dt);

    std::ofstream ofs;
    ofs.open(std::string("/tmp/out.kepler_kutta"));
    ofs << std::setprecision(std::numeric_limits<double>::max_digits10);
    ofs << std::scientific;
    ofs << "#  t  x(0)  x(1)  x(2)  x(3)  e\n";
    for (size_t i = 0; i < nsteps; ++i) {
        // sample
        if (i%smpfreq == 0) {
            ofs << t << " "
                << x(0) << " "
                << x(1) << " "
                << x(2) << " "
                << x(3) << " "
                << kepler.Energy(x) << "\n";
        }
        // x = x + dxdt * dt
        kutta.Step(kepler, x, dxdt, dt);
        t += dt;
    }
    ofs.close();
}

/// ---- Test Ode Gear Predictor-Corrector Integrator --------------------------
void TestOdeKeplerGear(double dt, size_t nsteps, size_t smpfreq)
{
    // Model parameters
    const double delta = 0.01;
    Kepler kepler(delta);

    // System parameters
    const size_t ndim = 4;
    const double e = 0.6;
    Vector<double> x(ndim);          // x(t)
    x(0) = 1.0 - e;
    x(1) = 0.0;
    x(2) = 0.0;
    x(3) = std::sqrt((1.0 + e)/(1.0 - e));
    Vector<double> dxdt(ndim);       // dx(t) / dt

    // Integrate system
    double t = 0.0;
    OdeGear<double, false> gear(ndim);
    gear.Init(kepler, x, dxdt, dt);

    std::ofstream ofs;
    ofs.open(std::string("/tmp/out.kepler_gear"));
    ofs << std::setprecision(std::numeric_limits<double>::max_digits10);
    ofs << std::scientific;
    ofs << "#  t  x(0)  x(1)  x(2)  x(3)  e\n";
    for (size_t i = 0; i < nsteps; ++i) {
        // sample
        if (i%smpfreq == 0) {
            ofs << t << " "
                << x(0) << " "
                << x(1) << " "
                << x(2) << " "
                << x(3) << " "
                << kepler.Energy(x) << "\n";
        }
        // x = x + dxdt * dt
        gear.Step(kepler, x, dxdt, dt);
        t += dt;
    }
    ofs.close();
}

/// ---- Test Ode Kepler solver ------------------------------------------------
TEST_CASE("Test Ode Kepler solver") {
    SECTION("OdeGauss") {
        TestOdeKeplerGauss(0.05, 20000, 2);
    }
    SECTION("OdeKutta") {
        TestOdeKeplerKutta(0.05, 20000, 2);
    }
    SECTION("OdeGear") {
        TestOdeKeplerGear( 0.05, 20000, 2);
    }
}
