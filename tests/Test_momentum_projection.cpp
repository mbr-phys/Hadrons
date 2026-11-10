/*
 * Test_momentum_projection.cpp, part of Hadrons
 *
 * Copyright (C) 2026
 *
 * Hadrons is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * Checks the batched momentum projection (MUtilities::MomentumProjectorList
 * plus the momProjector parameter) against the per-momentum path, for
 * MContraction::Meson (MSink::ScalarPoint sinks, both phase signs),
 * MContraction::StochasticMeson (same sinks with the source projected onto
 * the opposite momentum through sourceMom, non-zero source time; the
 * projected path builds its source phases separately) and
 * MContraction::ContractConservedSink
 * (moms), vector and axial, all directions. Meson projects all its gamma pairs
 * with one MomentumProjectMany call, so its comparison also checks that
 * packing; MomentumProjectMany is further checked directly, for scalar and
 * colour-matrix fields, against Project on each field and colour component.
 * The momenta include negative and
 * single-component ones, so a sign error in the phase is caught. Run it once
 * with the time direction split across ranks (e.g. --mpi 1.1.1.2) to check the
 * global sum in MomentumProject.
 */
#include <Hadrons/Application.hpp>
#include <Hadrons/Modules.hpp>

using namespace Grid;
using namespace Hadrons;

static const std::string dir   = "momentum_projection";
static const unsigned int traj = 0;

template <typename T>
static std::vector<T> readResult(const std::string stem, const std::string name)
{
    std::vector<T> v;
    ResultReader   reader(ModuleBase::resultFilename(dir + "/" + stem, traj));

    read(reader, name, v);

    return v;
}

static std::vector<int> momVec(const std::string &s)
{
    return strToVec<int>(s);
}

// |z| without std::abs, which does not accept Grid's Complex on GPU builds
static double complexAbs(const Complex &z)
{
    return std::sqrt(static_cast<double>(z.real()*z.real() + z.imag()*z.imag()));
}

// running max of |a - b| and of |a|
struct Diff
{
    double err{0.}, scale{0.};
    bool   missing{false};

    void add(const Complex &a, const Complex &b)
    {
        err   = std::max(err, complexAbs(a - b));
        scale = std::max(scale, complexAbs(a));
    }
    bool pass(const std::string what, const double tol)
    {
        const double rel = (scale > 0.) ? err/scale : err;

        LOG(Message) << what << ": max |projected - per momentum| / max |corr| = "
                     << rel << std::endl;
        if (missing)
        {
            LOG(Error) << what << ": a projected momentum or channel has no "
                       << "per-momentum counterpart" << std::endl;
        }

        return !missing and (rel < tol) and (scale > 0.);
    }
};

int main(int argc, char *argv[])
{
    // initialization //////////////////////////////////////////////////////////
    Grid_init(&argc, &argv);
    HadronsLogError.Active(GridLogError.isActive());
    HadronsLogWarning.Active(GridLogWarning.isActive());
    HadronsLogMessage.Active(GridLogMessage.isActive());
    HadronsLogIterative.Active(GridLogIterative.isActive());
    HadronsLogDebug.Active(GridLogDebug.isActive());
    LOG(Message) << "Grid initialized" << std::endl;

    const std::vector<std::string> momenta = {"0 0 0", "1 0 0", "0 -1 0",
                                              "1 1 0", "-1 0 1", "0 0 2"};
    const std::string  gammas     = "(Gamma5 Gamma5)(GammaX Gamma5)(GammaT GammaTGamma5)";
    const std::string  gammaSrcs  = "Gamma5 GammaX";
    const unsigned int sourceTime = 3;
    const double       tol        = 1.0e-10;

    // run setup ///////////////////////////////////////////////////////////////
    Application            application;
    Application::GlobalPar globalPar;

    globalPar.trajCounter.start = traj;
    globalPar.trajCounter.end   = traj + 1;
    globalPar.trajCounter.step  = 1;
    globalPar.runId             = "test_momentum_projection";
    application.setPar(globalPar);

    application.createModule<MGauge::Random>("gauge");

    MAction::Wilson::Par actionPar;
    actionPar.gauge    = "gauge";
    actionPar.mass     = 0.4;
    actionPar.boundary = "1 1 1 -1";
    actionPar.twist    = "0. 0. 0. 0.";
    application.createModule<MAction::Wilson>("wilson", actionPar);

    MSolver::RBPrecCG::Par solverPar;
    solverPar.action       = "wilson";
    solverPar.residual     = 1.0e-12;
    solverPar.maxIteration = 10000;
    application.createModule<MSolver::RBPrecCG>("cg", solverPar);

    // point-source propagator
    MSource::Point::Par ptPar;
    ptPar.position = "0 0 0 1";
    application.createModule<MSource::Point>("pt", ptPar);

    MFermion::GaugeProp::Par quarkPar;
    quarkPar.solver = "cg";
    quarkPar.source = "pt";
    application.createModule<MFermion::GaugeProp>("q", quarkPar);

    // independent Z2 wall sources and their solutions
    for (std::string n: {"A", "B"})
    {
        MSource::Z2::Par z2Par;
        z2Par.tA = sourceTime;
        z2Par.tB = sourceTime;
        application.createModule<MSource::Z2>("eta" + n, z2Par);

        MFermion::GaugeProp::Par phiPar;
        phiPar.solver = "cg";
        phiPar.source = "eta" + n;
        application.createModule<MFermion::GaugeProp>("phi" + n, phiPar);
    }

    // batched path
    // sign +1 is the MSink::Point phase; sign -1 with the negated momenta must
    // give the same correlators, labelled -p
    std::vector<std::string> negMomenta;
    for (auto &s: momenta)
    {
        std::vector<int> p = momVec(s);
        std::string      n;

        for (unsigned int mu = 0; mu < p.size(); ++mu)
        {
            n += ((mu > 0) ? " " : "") + std::to_string(-p[mu]);
        }
        negMomenta.push_back(n);
    }

    MUtilities::MomentumProjectorList::Par mpPar;
    mpPar.momenta = momenta;
    mpPar.sign    = 1;
    application.createModule<MUtilities::MomentumProjectorList>("mp", mpPar);
    mpPar.momenta = negMomenta;
    mpPar.sign    = -1;
    application.createModule<MUtilities::MomentumProjectorList>("mpNeg", mpPar);

    MContraction::Meson::Par mesPar;
    mesPar.q1           = "q";
    mesPar.q2           = "q";
    mesPar.gammas       = gammas;
    mesPar.momProjector = "mp";
    mesPar.output       = dir + "/meson_proj";
    application.createModule<MContraction::Meson>("meson_proj", mesPar);
    mesPar.momProjector = "mpNeg";
    mesPar.output       = dir + "/meson_projNeg";
    application.createModule<MContraction::Meson>("meson_projNeg", mesPar);

    MContraction::StochasticMeson::Par stPar;
    stPar.phiA         = "phiA";
    stPar.etaA         = "etaA";
    stPar.phiB         = "phiB";
    stPar.etaB         = "etaB";
    stPar.gammas       = gammas;
    stPar.sourceTime   = sourceTime;
    stPar.noiseA       = 0;
    stPar.noiseB       = 1;
    stPar.flowTime     = 0.;
    stPar.detSrc       = 0;
    stPar.momProjector = "mp";
    stPar.output       = dir + "/stoch_proj";
    application.createModule<MContraction::StochasticMeson>("stoch_proj", stPar);

    MContraction::ContractConservedSink::Par cPar;
    cPar.prop      = "q";
    cPar.action    = "wilson";
    cPar.source    = "pt";
    cPar.gammaSrcs = gammaSrcs;
    cPar.dirs      = {0, 1, 2, 3};
    for (auto current: {Current::Vector, Current::Axial})
    {
        const std::string c = (current == Current::Vector) ? "V" : "A";

        cPar.current      = current;
        cPar.moms         = {};
        cPar.momProjector = "mp";
        cPar.output       = dir + "/cons" + c + "_proj";
        application.createModule<MContraction::ContractConservedSink>("cons" + c + "_proj", cPar);

        // per-momentum path
        cPar.moms         = momenta;
        cPar.momProjector = "";
        cPar.output       = dir + "/cons" + c + "_mom";
        application.createModule<MContraction::ContractConservedSink>("cons" + c + "_mom", cPar);
    }

    // per-momentum path
    for (unsigned int m = 0; m < momenta.size(); ++m)
    {
        const std::string k = std::to_string(m);

        MSink::Point::Par sinkPar;
        sinkPar.mom = momenta[m];
        application.createModule<MSink::ScalarPoint>("snk" + k, sinkPar);

        mesPar.sink         = "snk" + k;
        mesPar.momProjector = "";
        mesPar.output       = dir + "/meson_" + k;
        application.createModule<MContraction::Meson>("meson_" + k, mesPar);

        stPar.sink         = "snk" + k;
        stPar.momProjector = "";
        stPar.sourceMom    = momenta[m];
        stPar.output       = dir + "/stoch_" + k;
        application.createModule<MContraction::StochasticMeson>("stoch_" + k, stPar);
    }

    // execution ///////////////////////////////////////////////////////////////
    application.saveParameterFile("MomentumProjection.xml");
    application.run();

    // comparison //////////////////////////////////////////////////////////////
    // every rank reads the files written by the boss rank
    Environment::getInstance().getGrid()->Barrier();

    bool pass = true;

    // MomentumProjectMany directly: several scalar fields, and colour-matrix
    // fields, in one projection each, against Project on every field and
    // colour component separately
    {
        GridCartesian   *grid = Environment::getInstance().getGrid();
        GridParallelRNG rng(grid);
        const Complex   i(0., 1.);
        LatticeComplex  coor(grid), ph(grid);
        std::vector<LatticeComplex> phases;

        rng.SeedFixedIntegers({5, 6, 7, 8});
        for (auto &s: momenta)
        {
            auto p = momVec(s);

            ph = Zero();
            for (unsigned int mu = 0; mu < p.size(); ++mu)
            {
                LatticeCoordinate(coor, mu);
                ph = ph + (static_cast<Real>(p[mu])/static_cast<Real>(grid->GlobalDimensions()[mu]))*coor;
            }
            phases.push_back(exp(static_cast<Real>(2*M_PI)*i*ph));
        }

        MomentumProjectMany<LatticeComplex, LatticeComplex> mp;
        mp.Allocate(momenta.size(), grid);
        mp.ImportMomenta(phases);

        const int nt = grid->GlobalDimensions()[Tp];
        std::vector<LatticeComplex>            scalars(3, grid);
        std::vector<LatticeColourMatrix>       matrices(2, grid);
        std::vector<std::vector<TComplex>>     projScalars;
        std::vector<std::vector<ColourMatrix>> projMatrices;
        std::vector<TComplex>                  ref;
        LatticeComplex                         comp(grid);
        Diff                                   dScalar, dMatrix;

        for (auto &f: scalars)  gaussian(rng, f);
        for (auto &f: matrices) gaussian(rng, f);
        mp.ProjectMany(scalars, projScalars);
        mp.ProjectMany(matrices, projMatrices);
        for (unsigned int f = 0; f < scalars.size(); ++f)
        {
            mp.Project(scalars[f], ref);
            for (unsigned int k = 0; k < momenta.size()*nt; ++k)
            {
                dScalar.add(TensorRemove(ref[k]), TensorRemove(projScalars[f][k]));
            }
        }
        for (unsigned int f = 0; f < matrices.size(); ++f)
        for (int a = 0; a < Nc; ++a)
        for (int b = 0; b < Nc; ++b)
        {
            comp = peekColour(matrices[f], a, b);
            mp.Project(comp, ref);
            for (unsigned int k = 0; k < momenta.size()*nt; ++k)
            {
                dMatrix.add(TensorRemove(ref[k]), projMatrices[f][k]()()(a, b));
            }
        }
        pass &= dScalar.pass("MomentumProjectMany, scalar fields", tol);
        pass &= dMatrix.pass("MomentumProjectMany, colour-matrix fields", tol);
    }

    {
        typedef MContraction::Meson::Result    Result;
        typedef MContraction::Meson::MomResult MomResult;
        auto proj    = readResult<MomResult>("meson_proj", "meson");
        auto projNeg = readResult<MomResult>("meson_projNeg", "meson");
        Diff d, dNeg;

        for (unsigned int m = 0; m < momenta.size(); ++m)
        {
            auto ref = readResult<Result>("meson_" + std::to_string(m), "meson");

            d.missing    |= (ref.size() != proj.size());
            dNeg.missing |= (ref.size() != projNeg.size());
            for (unsigned int g = 0; g < std::min(ref.size(), proj.size()); ++g)
            {
                d.missing |= (proj[g].momentum[m] != momVec(momenta[m]))
                             or (proj[g].gamma_snk != ref[g].gamma_snk)
                             or (proj[g].gamma_src != ref[g].gamma_src);
                for (unsigned int t = 0; t < ref[g].corr.size(); ++t)
                {
                    d.add(ref[g].corr[t], proj[g].corr[m][t]);
                }
            }
            for (unsigned int g = 0; g < std::min(ref.size(), projNeg.size()); ++g)
            {
                dNeg.missing |= (projNeg[g].momentum[m] != momVec(negMomenta[m]))
                                or (projNeg[g].gamma_snk != ref[g].gamma_snk)
                                or (projNeg[g].gamma_src != ref[g].gamma_src);
                for (unsigned int t = 0; t < ref[g].corr.size(); ++t)
                {
                    dNeg.add(ref[g].corr[t], projNeg[g].corr[m][t]);
                }
            }
        }
        pass &= d.pass("Meson, sign +1", tol);
        pass &= dNeg.pass("Meson, sign -1 at -p", tol);
    }
    {
        typedef MContraction::StochasticMeson::Result    Result;
        typedef MContraction::StochasticMeson::MomResult MomResult;
        auto proj = readResult<MomResult>("stoch_proj", "stochastic_meson");
        Diff d;

        for (unsigned int m = 0; m < momenta.size(); ++m)
        {
            auto ref = readResult<Result>("stoch_" + std::to_string(m), "stochastic_meson");

            d.missing |= (ref.size() != proj.size());
            for (unsigned int g = 0; g < std::min(ref.size(), proj.size()); ++g)
            {
                d.missing |= (proj[g].momentum[m] != momVec(momenta[m]))
                             or (proj[g].gamma_snk != ref[g].gamma_snk)
                             or (proj[g].gamma_src != ref[g].gamma_src)
                             or (proj[g].source_time != ref[g].source_time);
                for (unsigned int t = 0; t < ref[g].corr.size(); ++t)
                {
                    d.add(ref[g].corr[t], proj[g].corr[m][t]);
                }
            }
        }
        pass &= d.pass("StochasticMeson", tol);
    }
    for (std::string c: {"V", "A"})
    {
        typedef MContraction::ContractConservedSink::Result    Result;
        typedef MContraction::ContractConservedSink::MomResult MomResult;
        auto proj = readResult<MomResult>("cons" + c + "_proj", "conservedSink");
        auto ref  = readResult<Result>("cons" + c + "_mom", "conservedSink");
        Diff d;

        d.missing |= (ref.size() != proj.size()*momenta.size());
        for (auto &p: proj)
        for (unsigned int m = 0; m < momenta.size(); ++m)
        {
            auto r = std::find_if(ref.begin(), ref.end(), [&](const Result &x)
            {
                return (x.direction == p.direction) and (x.gamma_src == p.gamma_src)
                       and (momVec(x.momentum) == p.momentum[m]);
            });

            if (r == ref.end())
            {
                d.missing = true;
                continue;
            }
            for (unsigned int t = 0; t < r->corr.size(); ++t)
            {
                d.add(r->corr[t], p.corr[m][t]);
            }
        }
        pass &= d.pass("ContractConservedSink (" + c + ")", tol);
    }

    LOG(Message) << (pass ? "PASS" : "FAIL") << std::endl;
    LOG(Message) << "Grid is finalizing now" << std::endl;
    Grid_finalize();

    return pass ? EXIT_SUCCESS : EXIT_FAILURE;
}
