/*
 * Test_obc_flow_reflection.cpp, part of Hadrons
 *
 * Copyright (C) 2026
 *
 * Hadrons is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * Checks the openQCD open-boundary Wilson flow (flow cG = 1, force on the
 * spatial boundary links doubled; openQCD wflow.c and Luescher & Schaefer,
 * CPC 184 (2013) 519, App. B) against the periodic Wilson flow of the lattice
 * reflected across both boundaries. With Nt time slices the reflected lattice
 * has extent L = 2(Nt - 1), slice y0 >= Nt - 1 maps to x0 = L - y0, and the
 * reflected temporal links are U_t(y0) = U_t(L - y0 - 1)^dagger. The flows
 * must agree to rounding on every active link. Also checks that the unused
 * link T -> 0 is not changed, and that flow cG != 1 is detected.
 */
#include <Hadrons/Global.hpp>
#include <Hadrons/TimerArray.hpp>
#include <Hadrons/Modules/MGradientFlow/Utils_OBC.hpp>

using namespace Grid;
using namespace Hadrons;

typedef GIMPL                                          GImpl;
typedef typename GImpl::GaugeField                     GaugeField;
typedef typename GaugeField::vector_object::scalar_object SiteGauge;
typedef MGradientFlow::OBC::Evolution<MGradientFlow::OBC::WilsonAction<GImpl>,
                                      GImpl, FIMPL>    Evolution;

static constexpr int tDir = Nd - 1;

static GridCartesian *makeGrid(const int nt)
{
    // keep the SIMD layout spatial, so that any Nt is allowed
    Coordinate latt({4, 4, 4, nt});
    Coordinate simd(Nd, 1);
    int        n = vComplexD::Nsimd();

    for (int mu = 0; (mu < tDir) && (n > 1); ++mu)
    {
        simd[mu] = 2;
        n       /= 2;
    }
    if (n != 1)
    {
        std::cerr << "SIMD width too large for the 4^3 spatial volume" << std::endl;
        exit(EXIT_FAILURE);
    }

    return SpaceTimeGrid::makeFourDimGrid(latt, simd, GridDefaultMpi());
}

static void reflect(GaugeField &Ud, const GaugeField &Uo)
{
    const Coordinate dims = Ud.Grid()->GlobalDimensions();
    const int        nt   = Uo.Grid()->GlobalDimensions()[tDir];
    const int        L    = dims[tDir];
    Coordinate       y(Nd), x(Nd), xt(Nd);
    SiteGauge        so, st;

    for (int64_t idx = 0; idx < Ud.Grid()->gSites(); ++idx)
    {
        Lexicographic::CoorFromIndex(y, idx, dims);
        const int y0 = y[tDir];

        x       = y;
        x[tDir] = (y0 <= nt - 1) ? y0 : L - y0;
        peekSite(so, Uo, x);
        if (y0 >= nt - 1)
        {
            xt       = y;
            xt[tDir] = L - y0 - 1;
            peekSite(st, Uo, xt);
            so(tDir) = adj(st(tDir));
        }
        pokeSite(so, Ud, y);
    }
}

// max over active links of |U_open - U_reflected|
static RealD activeLinkDistance(const GaugeField &Uo, const GaugeField &Ud)
{
    const Coordinate dims = Uo.Grid()->GlobalDimensions();
    const int        nt   = dims[tDir];
    Coordinate       x(Nd);
    SiteGauge        so, sd;
    RealD            dist = 0.;

    for (int64_t idx = 0; idx < Uo.Grid()->gSites(); ++idx)
    {
        Lexicographic::CoorFromIndex(x, idx, dims);
        peekSite(so, Uo, x);
        peekSite(sd, Ud, x);
        for (int mu = 0; mu < Nd; ++mu)
        {
            if ((mu == tDir) && (x[tDir] == nt - 1))
            {
                continue;
            }
            dist = std::max(dist, std::sqrt(norm2(so(mu) - sd(mu))));
        }
    }

    return dist;
}

// max over the unused links T -> 0 of |U - U0|
static RealD unusedLinkDistance(const GaugeField &U, const GaugeField &U0)
{
    const Coordinate dims = U.Grid()->GlobalDimensions();
    Coordinate       x(Nd);
    SiteGauge        s, s0;
    RealD            dist = 0.;

    for (int64_t idx = 0; idx < U.Grid()->gSites(); ++idx)
    {
        Lexicographic::CoorFromIndex(x, idx, dims);
        if (x[tDir] != dims[tDir] - 1)
        {
            continue;
        }
        peekSite(s, U, x);
        peekSite(s0, U0, x);
        dist = std::max(dist, std::sqrt(norm2(s(tDir) - s0(tDir))));
    }

    return dist;
}

int main(int argc, char *argv[])
{
    Grid_init(&argc, &argv);

    const int    nt     = 6;
    const int    L      = 2*(nt - 1);
    const int    nStep  = 5;
    const double eps    = 0.02;
    const RealD  tol    = 1.0e-10;
    bool         pass   = true;

    GridCartesian    *gOpen = makeGrid(nt);
    GridCartesian    *gDbl  = makeGrid(L);
    GridParallelRNG  rng(gOpen);
    GaugeField       Uopen(gOpen), U0(gOpen), Uflow(gOpen), Udbl(gDbl);
    TimerArray       timer;

    rng.SeedFixedIntegers(std::vector<int>({11, 22, 33, 44}));
    SU<Nc>::HotConfiguration(rng, Uopen);
    U0 = Uopen;
    reflect(Udbl, Uopen);

    Evolution open(gOpen, 3.0, eps, -1.0, eps, GaugeBoundary::open, 1.0, &timer);
    Evolution periodic(gDbl, 3.0, eps, -1.0, eps, GaugeBoundary::periodic, 1.0, &timer);
    Evolution openCG(gOpen, 3.0, eps, -1.0, eps, GaugeBoundary::open, 0.5, &timer);

    Uflow = Uopen;
    for (int step = 0; step < nStep; ++step)
    {
        open.evolve_gauge(Uflow);
        periodic.evolve_gauge(Udbl);
    }

    const RealD dReflect = activeLinkDistance(Uflow, Udbl);
    const RealD dUnused  = unusedLinkDistance(Uflow, U0);

    std::cout << GridLogMessage << "open flow (cG = 1) vs reflected periodic flow: "
              << dReflect << std::endl;
    std::cout << GridLogMessage << "change of the unused links T -> 0: "
              << dUnused << std::endl;
    if (!(dReflect < tol))
    {
        std::cout << GridLogError << "FAIL: open flow differs from the reflected "
                  << "periodic flow" << std::endl;
        pass = false;
    }
    if (!(dUnused < tol))
    {
        std::cout << GridLogError << "FAIL: unused links changed" << std::endl;
        pass = false;
    }

    // sensitivity: flow cG = 1/2 must not reproduce the reflected flow
    Uflow = U0;
    for (int step = 0; step < nStep; ++step)
    {
        openCG.evolve_gauge(Uflow);
    }

    const RealD dCG = activeLinkDistance(Uflow, Udbl);

    std::cout << GridLogMessage << "open flow (cG = 1/2) vs reflected periodic flow: "
              << dCG << std::endl;
    if (!(dCG > 1.0e3*tol))
    {
        std::cout << GridLogError << "FAIL: the test does not resolve the flow cG"
                  << std::endl;
        pass = false;
    }

    std::cout << GridLogMessage << (pass ? "PASS" : "FAIL") << std::endl;
    Grid_finalize();

    return pass ? EXIT_SUCCESS : EXIT_FAILURE;
}
