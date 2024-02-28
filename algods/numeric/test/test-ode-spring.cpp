//
// test-ode-spring.cpp
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
#include "algods/numeric/linalg.h"
#include "algods/numeric/ode.h"

/// ---- Harmonic spring model -------------------------------------------------
struct Spring {
    double kappa_ = 0.0;
    double alpha_ = 0.0;

    Spring(double kappa, double alpha)
        : kappa_(kappa)
        , alpha_(alpha) {}
    ~Spring() = default;

    void operator()(Vector<double> & __restrict__ x,
                    Vector<double> & __restrict__ dxdt) {
        dxdt(0) = x(1);
        dxdt(1) = -kappa_*x(0) - 2.0*alpha_*x(1);
    }

    double Energy(Vector<double> & __restrict__ x) {
        return (0.5*kappa_*x(0)*x(0) + 0.5*x(1)*x(1));
    }

    void Compute(double t, double v0, Vector<double> & __restrict__ x) {
        if (std::fabs(alpha_) > 0.0 && alpha_*alpha_ > kappa_) {
            // Overdamped system
            double w1 = (-std::sqrt(2.0) + 1.0)*std::sqrt(kappa_);
            double w2 = (-std::sqrt(2.0) - 1.0)*std::sqrt(kappa_);

            x(0) = v0 * (std::exp(w1*t) - std::exp(w2*t)) / (w1 - w2);
            x(1) = v0 * (w1*std::exp(w1*t) - w2*std::exp(w2*t)) / (w1 - w2);
        } else if (std::fabs(alpha_) > 0.0 && alpha_*alpha_ < kappa_) {
            // Underdamped system
            double mu = std::sqrt(kappa_ / 2.0);

            x(0) = v0 * std::exp(-mu*t) * std::sin(mu*t) / mu;
            x(1) = v0 * std::exp(-mu*t) * (std::cos(mu*t)-std::sin(mu*t));
        } else {
            // Oscillatory system
            double mu = std::sqrt(kappa_);

            x(0) = v0 * std::sin(mu*t) / mu;
            x(1) = v0 * std::cos(mu*t);
        }
    }
};

/// ---- Test Ode Euler Integrator ---------------------------------------------
void TestOdeSpringEuler(double dt, size_t nsteps, size_t smpfreq)
{
    // Model parameters
    const double kappa = 2.0*M_PI*0.5;

    // System parameters
    double t = 0.0;
    const double x0 = 0.0;
    const double v0 = 1.0;

    Vector<double> x(2);         // x(t)
    Vector<double> dxdt(2);      // dx(t) / dt
    Vector<double> xref(2);      // reference x(t)

    // Run function.
    auto Run = [&] (const std::string &filename, const double m) -> void
    {
        t = 0.0;
        x(0) = x0;                                  // x(t=0)
        x(1) = v0;                                  // v(t=0)
        CopyVector<double, false>(x, xref);       // xref(t=0)

        Spring spring(kappa, std::sqrt(m * kappa));
        OdeEuler<double, false> euler;
        euler.Init(spring, x, dxdt, dt);

        std::ofstream ofs;
        ofs.open(filename);
        ofs << std::setprecision(std::numeric_limits<double>::max_digits10);
        ofs << std::scientific;
        ofs << "#  t  x(0)  x(1)  e1  x  v  e2\n";

        for (size_t i = 0; i < nsteps; ++i) {
            // sample
            if (i%smpfreq == 0) {
                ofs << t << " "
                    << x(0) << " "
                    << x(1) << " "
                    << spring.Energy(x) << " "
                    << xref(0) << " "
                    << xref(1) << " "
                    << spring.Energy(xref) << "\n";
            }

            // x = x + dxdt * dt
            euler.Step(spring, x, dxdt, dt);
            t += dt;

            // xref
            spring.Compute(t, v0, xref);
        }

        ofs.close();
    };

    Run("/tmp/out.euler1", 2.0);    // Damped non oscillatory spring
    Run("/tmp/out.euler2", 0.5);    // Damped oscillatory spring
    Run("/tmp/out.euler3", 0.0);    // Non damped oscillatory spring
}


/// ---- Test Ode Gauss Integrator ---------------------------------------------
void TestOdeSpringGauss(double dt, size_t nsteps, size_t smpfreq)
{
    // Model parameters
    const double kappa = 2.0*M_PI*0.5;

    // System parameters
    double t = 0.0;
    const double x0 = 0.0;
    const double v0 = 1.0;

    Vector<double> x(2);         // x(t)
    Vector<double> dxdt(2);      // dx(t) / dt
    Vector<double> xref(2);      // reference x(t)

    // Run function.
    auto Run = [&] (const std::string &filename, const double m) -> void
    {
        // Solver parameters
        t = 0.0;
        x(0) = x0;                                  // x(t=0)
        x(1) = v0;                                  // v(t=0)
        CopyVector<double, false>(x, xref);       // xref(t=0)

        Spring spring(kappa, std::sqrt(m * kappa));
        OdeGauss<double, false> gauss(2, 1.0E-12, 20);
        gauss.Init(spring, x, dxdt, dt);

        std::ofstream ofs;
        ofs.open(filename);
        ofs << std::setprecision(std::numeric_limits<double>::max_digits10);
        ofs << std::scientific;
        ofs << "#  t  x(0)  x(1)  e1  x  v  e2\n";

        for (size_t i = 0; i < nsteps; ++i) {
            // sample
            if (i%smpfreq == 0) {
                ofs << t << " "
                    << x(0) << " "
                    << x(1) << " "
                    << spring.Energy(x) << " "
                    << xref(0) << " "
                    << xref(1) << " "
                    << spring.Energy(xref) << "\n";
            }

            // x = x + dxdt * dt
            gauss.Step(spring, x, dxdt, dt);
            t += dt;

            // xref
            spring.Compute(t, v0, xref);
        }

        ofs.close();
    };

    Run("/tmp/out.gauss1", 2.0);    // Damped non oscillatory spring
    Run("/tmp/out.gauss2", 0.5);    // Damped oscillatory spring
    Run("/tmp/out.gauss3", 0.0);    // Non damped oscillatory spring
}

/// ---- Test Ode Runge-Kutta Integrator ---------------------------------------
void TestOdeSpringKutta(double dt, size_t nsteps, size_t smpfreq)
{
    // Model parameters
    const double kappa = 2.0*M_PI*0.5;

    // System parameters
    double t = 0.0;
    const double x0 = 0.0;
    const double v0 = 1.0;

    Vector<double> x(2);        // x(t)
    Vector<double> dxdt(2);     // dx(t) / dt
    Vector<double> xref(2);     // reference x(t)

    // Run function.
    auto Run = [&] (const std::string &filename, const double m) -> void
    {
        // Solver parameters
        t = 0.0;
        x(0) = x0;                                  // x(t=0)
        x(1) = v0;                                  // v(t=0)
        CopyVector<double, false>(x, xref);       // xref(t=0)

        Spring spring(kappa, std::sqrt(m * kappa));
        OdeKutta<double, false> kutta(2);
        kutta.Init(spring, x, dxdt, dt);

        std::ofstream ofs;
        ofs.open(filename);
        ofs << std::setprecision(std::numeric_limits<double>::max_digits10);
        ofs << std::scientific;
        ofs << "#  t  x(0)  x(1)  e1  x  v  e2\n";

        for (size_t i = 0; i < nsteps; ++i) {
            // sample
            if (i%smpfreq == 0) {
                ofs << t << " "
                    << x(0) << " "
                    << x(1) << " "
                    << spring.Energy(x) << " "
                    << xref(0) << " "
                    << xref(1) << " "
                    << spring.Energy(xref) << "\n";
            }

            // x = x + dxdt * dt
            kutta.Step(spring, x, dxdt, dt);
            t += dt;

            // xref
            spring.Compute(t, v0, xref);
        }

        ofs.close();
    };

    Run("/tmp/out.kutta1", 2.0);    // Damped non oscillatory spring
    Run("/tmp/out.kutta2", 0.5);    // Damped oscillatory spring
    Run("/tmp/out.kutta3", 0.0);    // Non damped oscillatory spring
}

/// ---- Test Ode Gear Predictor-Corrector Integrator --------------------------
void TestOdeSpringGear(double dt, size_t nsteps, size_t smpfreq)
{
    // Model parameters
    const double kappa = 2.0*M_PI*0.5;

    // System parameters
    double t = 0.0;
    const double x0 = 0.0;
    const double v0 = 1.0;

    Vector<double> x(2);        // x(t)
    Vector<double> dxdt(2);     // dx(t) / dt
    Vector<double> xref(2);     // reference x(t)

    // Run function.
    auto Run = [&] (const std::string &filename, const double m) -> void
    {
        // Solver parameters
        t = 0.0;
        x(0) = x0;                                  // x(t=0)
        x(1) = v0;                                  // v(t=0)
        CopyVector<double, false>(x, xref);       // xref(t=0)

        Spring spring(kappa, std::sqrt(m * kappa));
        OdeGear<double, false> gear(2);
        gear.Init(spring, x, dxdt, dt);

        std::ofstream ofs;
        ofs.open(filename);
        ofs << std::setprecision(std::numeric_limits<double>::max_digits10);
        ofs << std::scientific;
        ofs << "#  t  x(0)  x(1)  e1  x  v  e2\n";

        for (size_t i = 0; i < nsteps; ++i) {
            // sample
            if (i%smpfreq == 0) {
                ofs << t << " "
                    << x(0) << " "
                    << x(1) << " "
                    << spring.Energy(x) << " "
                    << xref(0) << " "
                    << xref(1) << " "
                    << spring.Energy(xref) << "\n";
            }

            // x = x + dxdt * dt
            gear.Step(spring, x, dxdt, dt);
            t += dt;

            // xref
            spring.Compute(t, v0, xref);
        }

        ofs.close();
    };

    Run("/tmp/out.kutta1", 2.0);    // Damped non oscillatory spring
    Run("/tmp/out.kutta2", 0.5);    // Damped oscillatory spring
    Run("/tmp/out.kutta3", 0.0);    // Non damped oscillatory spring
}

/// ---- Test Ode spring solver ------------------------------------------------
TEST_CASE("Test Ode spring solver") {
    SECTION("OdeEuler") {
        TestOdeSpringEuler(0.00001, 10000000, 1000);
    }
    SECTION("OdeGauss") {
        TestOdeSpringGauss(0.001, 100000, 10);
    }
    SECTION("OdeKutta") {
        TestOdeSpringKutta(0.001, 100000, 10);
    }
    SECTION("OdeGear") {
        TestOdeSpringGear( 0.001, 100000, 10);
    }
}
