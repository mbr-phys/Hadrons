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
 * MContraction::Meson (MSink::ScalarPoint sinks), MContraction::StochasticMeson
 * (same sinks, non-zero source time) and MContraction::ContractConservedSink
 * (moms), vector and axial, all directions. The momenta include negative and
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
    MUtilities::MomentumProjectorList::Par mpPar;
    mpPar.momenta = momenta;
    application.createModule<MUtilities::MomentumProjectorList>("mp", mpPar);

    MContraction::Meson::Par mesPar;
    mesPar.q1           = "q";
    mesPar.q2           = "q";
    mesPar.gammas       = gammas;
    mesPar.momProjector = "mp";
    mesPar.output       = dir + "/meson_proj";
    application.createModule<MContraction::Meson>("meson_proj", mesPar);

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

    {
        typedef MContraction::Meson::Result    Result;
        typedef MContraction::Meson::MomResult MomResult;
        auto proj = readResult<MomResult>("meson_proj", "meson");
        Diff d;

        for (unsigned int m = 0; m < momenta.size(); ++m)
        {
            auto ref = readResult<Result>("meson_" + std::to_string(m), "meson");

            d.missing |= (ref.size() != proj.size());
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
        }
        pass &= d.pass("Meson", tol);
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
