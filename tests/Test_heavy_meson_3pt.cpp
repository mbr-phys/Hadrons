/*
 * Test_heavy_meson_3pt.cpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2026
 *
 * Author: Matthew Black <matthewkblack@protonmail.com>
 *
 * Hadrons is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * Hadrons is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Hadrons.  If not, see <http://www.gnu.org/licenses/>.
 *
 * See the full license in the file "LICENSE" in the top level distribution
 * directory.
 */

/*  END LEGAL */
#include <Hadrons/Application.hpp>
#include <Hadrons/Modules.hpp>

using namespace Grid;
using namespace Hadrons;

/*
 * Validation workflow for MContraction::HeavyMeson3pt (two-quark HQE
 * operators), MContraction::WeakNonEye3ptDeriv (four-quark operators with
 * gamma pairs and derivatives) and the HQE insertion modules
 * (MUtilities::CovariantDerivative, MGauge::FieldStrength,
 * MUtilities::GaugeInsertion). Requires Nt > 4 (the sink is at t = 4).
 *
 * Point sources are used at both ends, so every comparison below is an exact
 * identity configuration by configuration (no noise average). Gamma3pt
 * results are indexed corr[gamma][momentum][time], HeavyMeson3pt results
 * corr[triad][momentum][time], WeakNonEye3ptDeriv results
 * contractions[i].corr[momentum][time] and WeakNonEye3pt results
 * [i].corr[time] (entries ordered G, then trace 1, 2).
 *
 *  1. 'hm3pt' vertices Gamma*            == 'ref_gamma3pt' (same triads)
 *     'hm3pt_mom' vertices Gamma*        == 'ref_gamma3pt_mom', all momenta
 *  2. 'hm3pt' D4                         == 'ref_D4'
 *     'hm3pt' D4                         == 0.5*'ref_rhq_D4' except at the
 *                                           operator times t = 0, Nt-1 (RHQ
 *                                           insertions use periodic links)
 *  3. 'hm3pt' D2                         == 'ref_D2'
 *     'hm3pt_mom' at p = 0               == 'hm3pt' for every entry
 *  4. 'hm3pt' sigmaB, sigmaE, divE       == 'ref_sigmaB', 'ref_sigmaE', 'ref_divE'
 *  5. 'wne_DeltaF2', deriv 'none'        == 'ref_weakNonEye' (16 G x G,
 *                                           trace 1 and 2)
 *     'wne_DeltaF2', deriv 'qLeft@3'     == 'ref_wne_qLeft3'
 *     'wne_DeltaF2', deriv 'qBarRight@3' == 'ref_wne_qBarRight3'
 *     'wne_DeltaF2', deriv 'qBarLeft qLeft@0123'
 *                                        == sum_rho 'ref_wne_A_<rho>'
 *     'wne_mom_*' at p = 0               == 'wne_*' for every entry
 *  6. '*_g' (random gauge transformation) == the same instance on 'gauge',
 *     all entries (hm3pt_g, wne_g_DeltaF2, wne_g_DeltaF0_T1)
 *  7. 'hm3pt_unit' (unit gauge): sigmaB, sigmaE, sigmaBD4, sigmaED4, DiE, divE,
 *     DiSigmaE, gSigmaB, gSigmaE vanish
 *  8. 'hm3pt': D4 + D4L should be close to the time derivative of the
 *     Identity-vertex correlator; hermiticity/reality checks as in
 *     lifetimes_round2.md §7
 *  9. chirality-even vertices:
 *     'hm3pt' gD4, gD2, gD4D4            == 'ref_gD4', 'ref_gD2', 'ref_gD4D4'
 *     'hm3pt' gSigmaB, gSigmaE           == 'ref_gSigmaB', 'ref_gSigmaE'
 *     'hm3pt' gDi                        == sum_i 'ref_gDi_<i>'
 *     'hm3pt' gDiD4                      == sum_i 'ref_gDiD4_<i>'
 * 10. WeakEye3pt gamma pairs (loop = l_z2 adj(z2) from a Z2 wall at t = 2,
 *     loop_gT = loop gamma_4 from the noise gamma_4 z2):
 *     'we_all' diagonal entries          == 'we_diag' (both traces)
 *     'we_pair' (Identity GammaT), trace 1 and 2
 *                                        == 'we_diag_gT' Identity, trace 1 and 2
 *     'we_pair' (GammaT Identity), trace 1 and 2
 *                                        == 'we_diag_gT' GammaT, trace 1 and 2
 */

struct Chain
{
    std::string hSrc, hSnk, lSrc, lSrcF, lSnkF, flowU, ptSnk;
};

static const unsigned int tSnk = 4;

// propagators from two point sources on 'gauge', flowed to tau = 0.02
Chain makeChain(Application &application, const std::string &gauge,
                const std::string &sfx)
{
    std::string boundary = "1 1 1 -1";
    std::string twist    = "0. 0. 0. 0.";
    Chain       c;

    // point sources
    MSource::Point::Par ptPar;
    ptPar.position = "0 0 0 0";
    application.createModule<MSource::Point>("pt_src" + sfx, ptPar);
    ptPar.position = "1 1 1 " + std::to_string(tSnk);
    application.createModule<MSource::Point>("pt_snk" + sfx, ptPar);
    c.ptSnk = "pt_snk" + sfx;

    // heavy and light actions/solvers
    const std::vector<std::pair<std::string, double>> quarks = {{"h", 0.4}, {"l", 0.05}};
    for (auto &qm: quarks)
    {
        MAction::DWF::Par actionPar;
        actionPar.gauge    = gauge;
        actionPar.Ls       = 8;
        actionPar.M5       = 1.8;
        actionPar.mass     = qm.second;
        actionPar.boundary = boundary;
        actionPar.twist    = twist;
        application.createModule<MAction::DWF>("DWF_" + qm.first + sfx, actionPar);

        MSolver::RBPrecCG::Par solverPar;
        solverPar.action       = "DWF_" + qm.first + sfx;
        solverPar.residual     = 1.0e-10;
        solverPar.maxIteration = 10000;
        application.createModule<MSolver::RBPrecCG>("CG_" + qm.first + sfx, solverPar);

        for (std::string end: {"src", "snk"})
        {
            MFermion::GaugeProp::Par quarkPar;
            quarkPar.solver = "CG_" + qm.first + sfx;
            quarkPar.source = "pt_" + end + sfx;
            application.createModule<MFermion::GaugeProp>(qm.first + "_" + end + sfx, quarkPar);
        }
    }
    c.lSrc = "l_src" + sfx;

    // flow all four propagators and the gauge field
    MGradientFlow::FermionFlow::Par flowPar;
    flowPar.gauge         = gauge;
    flowPar.steps         = 2;
    flowPar.step_size     = 0.01;
    flowPar.meas_interval = 2;
    flowPar.props         = {"h_src" + sfx, "h_snk" + sfx, "l_src" + sfx, "l_snk" + sfx};
    flowPar.defaultType   = "PropagatorField";
    flowPar.outProps      = {"h_src_f" + sfx, "h_snk_f" + sfx, "l_src_f" + sfx, "l_snk_f" + sfx};
    flowPar.bc            = -1;
    application.createModule<MGradientFlow::FermionFlow>("flow" + sfx, flowPar);
    c.lSrcF = "l_src_f" + sfx;
    c.lSnkF = "l_snk_f" + sfx;
    c.hSrc  = "h_src_f" + sfx;
    c.hSnk  = "h_snk_f" + sfx;
    c.flowU = "flow" + sfx + "_U";

    return c;
}

// HeavyMeson3pt with all gamma vertices and all HQE operators
std::string makeHM3pt(Application &application, const Chain &c,
                      const std::string &name, const std::string &momProjector = "")
{
    MContraction::HeavyMeson3pt::Par par;
    par.gauge        = c.flowU;
    par.bc           = -1;
    par.flowTime     = 0.02;
    par.momProjector = momProjector;
    par.qSrc         = c.hSrc;
    par.qSnk         = c.hSnk;
    par.spectator    = c.lSrc;
    par.sinkNoise    = c.ptSnk;
    par.tSnk         = tSnk;
    par.gamma        = {"Gamma5 all Gamma5", "Gamma5 ops Gamma5"};
    application.createModule<MContraction::HeavyMeson3pt>(name, par);

    return name;
}

MContraction::WeakDerivSpec deriv(const std::string &legs, const std::string &dirs,
                                  const std::string &couple)
{
    MContraction::WeakDerivSpec d;

    d.legs   = legs;
    d.dirs   = dirs;
    d.couple = couple;

    return d;
}

// WeakNonEye3ptDeriv: DeltaF = 2 with derivative specifications, and DeltaF = 0
// (T1 legs) with off-diagonal gamma pairs; returns the two module names
std::vector<std::string> makeWNEDeriv(Application &application, const Chain &c,
                                      const std::string &prefix,
                                      const std::string &momProjector = "")
{
    MContraction::WeakNonEye3ptDeriv::Par par;
    par.gauge        = c.flowU;
    par.bc           = -1;
    par.momProjector = momProjector;
    par.gammaIn      = Gamma::Algebra::Gamma5;
    par.gammaOut     = Gamma::Algebra::Gamma5;

    par.qBarLeft    = c.hSrc;
    par.qLeft       = c.lSrcF;
    par.qBarRight   = c.lSnkF;
    par.qRight      = c.hSnk;
    par.gammaPairs  = "diag";
    par.derivatives = {deriv("none", "", ""),
                       // leading power: time derivative on the light legs
                       deriv("qLeft", "3", "none"),
                       deriv("qBarRight", "3", "none"),
                       // full R2/R3-type: derivative on qbar and q of one bilinear
                       deriv("qBarLeft qLeft", "0 1 2 3", "none"),
                       deriv("qRight qLeft", "0 1 2 3", "none"),
                       deriv("qBarRight qRight", "0 1 2 3", "none"),
                       deriv("qBarLeft qBarRight", "0 1 2 3", "none"),
                       // index coupled to the second gamma
                       deriv("qLeft", "0 1 2 3", "G2")};
    application.createModule<MContraction::WeakNonEye3ptDeriv>(prefix + "_DeltaF2", par);

    par.qBarLeft    = c.lSrcF;
    par.qLeft       = c.hSrc;
    par.qBarRight   = c.lSnkF;
    par.qRight      = c.hSnk;
    par.gammaPairs  = "(Identity GammaT)(GammaT Identity)";
    par.derivatives = {deriv("none", "", "")};
    application.createModule<MContraction::WeakNonEye3ptDeriv>(prefix + "_DeltaF0_T1", par);

    return {prefix + "_DeltaF2", prefix + "_DeltaF0_T1"};
}

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

    // run setup ///////////////////////////////////////////////////////////////
    Application              application;
    std::vector<std::string> results;

    // global parameters
    Application::GlobalPar globalPar;
    globalPar.trajCounter.start    = 1500;
    globalPar.trajCounter.end      = 1520;
    globalPar.trajCounter.step     = 20;
    globalPar.runId                = "test_heavy_meson_3pt";
    globalPar.genetic.maxGen       = 1000;
    globalPar.genetic.maxCstGen    = 200;
    globalPar.genetic.popSize      = 20;
    globalPar.genetic.mutationRate = 0.1;
    application.setPar(globalPar);

    // gauge fields: random, random gauge transform of it, unit
    application.createModule<MGauge::Random>("gauge");
    application.createModule<MGauge::RandomGaugeTransformation>("gt");
    MGauge::GaugeTransform::Par gtPar;
    gtPar.gauge     = "gauge";
    gtPar.transform = "gt";
    application.createModule<MGauge::GaugeTransform>("gauge_g", gtPar);
    application.createModule<MGauge::Unit>("gauge_unit");

    Chain c  = makeChain(application, "gauge", "");
    Chain cg = makeChain(application, "gauge_g", "_g");
    Chain cu = makeChain(application, "gauge_unit", "_unit");

    MUtilities::MPScalar::Par momProjPar;
    momProjPar.maxFourier = 1;
    application.createModule<MUtilities::MPScalar>("phases", momProjPar);

    results.push_back(makeHM3pt(application, c, "hm3pt"));
    results.push_back(makeHM3pt(application, c, "hm3pt_mom", "phases"));
    results.push_back(makeHM3pt(application, cg, "hm3pt_g"));
    results.push_back(makeHM3pt(application, cu, "hm3pt_unit"));
    for (auto &n: makeWNEDeriv(application, c, "wne"))
    {
        results.push_back(n);
    }
    for (auto &n: makeWNEDeriv(application, c, "wne_mom", "phases"))
    {
        results.push_back(n);
    }
    for (auto &n: makeWNEDeriv(application, cg, "wne_g"))
    {
        results.push_back(n);
    }

    // references on the random gauge field ////////////////////////////////////
    // spectator closed at the sink point
    MSink::SourceSink::Par ssPar;
    ssPar.source   = c.ptSnk;
    ssPar.sourceOp = "none";
    application.createModule<MSink::SourceSink>("snk_pt", ssPar);
    MSink::Smear::Par smPar;
    smPar.q    = c.lSrc;
    smPar.sink = "snk_pt";
    application.createModule<MSink::Smear>("l_src_snk", smPar);

    auto gamma3pt = [&](const std::string &name, const std::string &q3,
                        const std::string &gamma, const std::string &momProjector = "")
    {
        MContraction::Gamma3pt::Par g3Par;
        g3Par.q1           = "l_src_snk";
        g3Par.q2           = c.hSrc;
        g3Par.q3           = q3;
        g3Par.gamma        = {"Gamma5 " + gamma + " Gamma5"};
        g3Par.tSnk         = tSnk;
        g3Par.momProjector = momProjector;
        application.createModule<MContraction::Gamma3pt>(name, g3Par);
        results.push_back(name);
    };

    // 1. plain gamma insertions, at zero and nonzero momentum
    gamma3pt("ref_gamma3pt", c.hSnk, "all");
    gamma3pt("ref_gamma3pt_mom", c.hSnk, "all", "phases");

    // 2. time derivative (new module and RHQ insertion = 2 D_4 with periodic links)
    MUtilities::CovariantDerivativePropagator::Par cdPar;
    cdPar.q     = c.hSnk;
    cdPar.gauge = c.flowU;
    cdPar.type  = "symmetric";
    cdPar.dirs  = "3";
    cdPar.bc    = -1;
    application.createModule<MUtilities::CovariantDerivativePropagator>("h_snk_D4", cdPar);
    gamma3pt("ref_D4", "h_snk_D4", "Identity");

    MRHQ::RHQInsertionI::Par rhqPar;
    rhqPar.q     = c.hSnk;
    rhqPar.index = "3";
    rhqPar.gamma = Gamma::Algebra::Identity;
    rhqPar.gauge = c.flowU;
    application.createModule<MRHQ::RHQInsertionI>("h_snk_rhqD4", rhqPar);
    gamma3pt("ref_rhq_D4", "h_snk_rhqD4", "Identity");

    // 3. spatial Laplacian applied directly
    cdPar.type = "laplacian";
    cdPar.dirs = "0 1 2";
    application.createModule<MUtilities::CovariantDerivativePropagator>("h_snk_D2", cdPar);
    gamma3pt("ref_D2", "h_snk_D2", "Identity");

    // 4. gluonic insertions
    MUtilities::GaugeInsertionPropagator::Par giPar;
    giPar.q     = c.hSnk;
    giPar.gauge = c.flowU;
    giPar.mode  = "sigmaB";
    giPar.gamma = Gamma::Algebra::Identity;
    application.createModule<MUtilities::GaugeInsertionPropagator>("h_snk_sigmaB", giPar);
    gamma3pt("ref_sigmaB", "h_snk_sigmaB", "Identity");

    giPar.mode = "sigmaE";
    application.createModule<MUtilities::GaugeInsertionPropagator>("h_snk_sigmaE", giPar);
    gamma3pt("ref_sigmaE", "h_snk_sigmaE", "Identity");

    MGauge::FieldStrength::Par fsPar;
    fsPar.gauge = c.flowU;
    fsPar.mode  = "divE";
    application.createModule<MGauge::FieldStrength>("divE", fsPar);
    giPar.mode   = "matrix";
    giPar.matrix = "divE";
    application.createModule<MUtilities::GaugeInsertionPropagator>("h_snk_divE", giPar);
    gamma3pt("ref_divE", "h_snk_divE", "Identity");

    // 5. four-quark references: WeakNonEye3pt on (derived) legs
    auto weakNonEye = [&](const std::string &name, const std::string &qbl,
                          const std::string &ql, const std::string &qbr,
                          const std::string &qr)
    {
        MContraction::WeakNonEye3pt::Par wnePar;
        wnePar.qBarLeft  = qbl;
        wnePar.qLeft     = ql;
        wnePar.qBarRight = qbr;
        wnePar.qRight    = qr;
        wnePar.gammaIn   = Gamma::Algebra::Gamma5;
        wnePar.gammaOut  = Gamma::Algebra::Gamma5;
        application.createModule<MContraction::WeakNonEye3pt>(name, wnePar);
        results.push_back(name);
    };
    // symmetric derivative of a flowed leg in direction rho (created once)
    std::set<std::string> derived;
    auto derivLeg = [&](const std::string &q, const int rho)
    {
        std::string name = q + "_D" + std::to_string(rho);

        if (derived.insert(name).second)
        {
            cdPar.q    = q;
            cdPar.type = "symmetric";
            cdPar.dirs = std::to_string(rho);
            application.createModule<MUtilities::CovariantDerivativePropagator>(name, cdPar);
        }

        return name;
    };

    // 9. chirality-even vertices: gamma_4 on the derived legs above, and
    //    gamma_i on D_i (of q3 or of D_4 q3), one module per i
    gamma3pt("ref_gD4", "h_snk_D4", "GammaT");
    gamma3pt("ref_gD2", "h_snk_D2", "GammaT");
    gamma3pt("ref_gSigmaB", "h_snk_sigmaB", "GammaT");
    gamma3pt("ref_gSigmaE", "h_snk_sigmaE", "GammaT");
    cdPar.q    = c.hSnk;
    cdPar.type = "laplacian";
    cdPar.dirs = "3";
    application.createModule<MUtilities::CovariantDerivativePropagator>("h_snk_D4D4", cdPar);
    gamma3pt("ref_gD4D4", "h_snk_D4D4", "GammaT");
    const std::vector<std::string> gSpatial = {"GammaX", "GammaY", "GammaZ"};
    for (int i = 0; i < Tp; ++i)
    {
        gamma3pt("ref_gDi_" + std::to_string(i), derivLeg(c.hSnk, i), gSpatial[i]);
        gamma3pt("ref_gDiD4_" + std::to_string(i), derivLeg("h_snk_D4", i), gSpatial[i]);
    }

    weakNonEye("ref_weakNonEye", c.hSrc, c.lSrcF, c.lSnkF, c.hSnk);
    weakNonEye("ref_wne_qLeft3", c.hSrc, derivLeg(c.lSrcF, 3), c.lSnkF, c.hSnk);
    weakNonEye("ref_wne_qBarRight3", c.hSrc, c.lSrcF, derivLeg(c.lSnkF, 3), c.hSnk);
    // spec A: derivative on the adjoint leg qBarLeft and on qLeft, one module per rho
    for (int rho = 0; rho < Nd; ++rho)
    {
        weakNonEye("ref_wne_A_" + std::to_string(rho), derivLeg(c.hSrc, rho),
                   derivLeg(c.lSrcF, rho), c.lSnkF, c.hSnk);
    }

    // 10. eye contractions with gamma pairs
    MSource::Z2::Par z2Par;
    z2Par.tA = 2;
    z2Par.tB = 2;
    application.createModule<MSource::Z2>("z2_loop", z2Par);
    MFermion::GaugeProp::Par loopPropPar;
    loopPropPar.solver = "CG_l";
    loopPropPar.source = "z2_loop";
    application.createModule<MFermion::GaugeProp>("l_z2", loopPropPar);
    MContraction::Loop::Par loopPar;
    loopPar.q   = "l_z2";
    loopPar.eta = "z2_loop";
    application.createModule<MContraction::Loop>("loop", loopPar);
    // gamma_4 z2 as noise: loop_gT = l_z2 adj(gamma_4 z2) = loop gamma_4
    MSource::SeqGamma::Par sgPar;
    sgPar.q     = "z2_loop";
    sgPar.tA    = 0;
    sgPar.tB    = GridDefaultLatt()[Tp] - 1;
    sgPar.gamma = Gamma::Algebra::GammaT;
    sgPar.mom   = "0 0 0 0";
    application.createModule<MSource::SeqGamma>("z2_loop_gT", sgPar);
    loopPar.eta = "z2_loop_gT";
    application.createModule<MContraction::Loop>("loop_gT", loopPar);

    auto weakEye = [&](const std::string &name, const std::string &loop,
                       const std::string &gammaPairs)
    {
        MContraction::WeakEye3pt::Par wePar;
        wePar.qBarLeft   = c.hSrc;
        wePar.qBarRight  = c.hSnk;
        wePar.qSpectator = "l_src_snk";
        wePar.loop       = loop;
        wePar.tOut       = tSnk;
        wePar.gammaIn    = Gamma::Algebra::Gamma5;
        wePar.gammaOut   = Gamma::Algebra::Gamma5;
        wePar.gammaPairs = gammaPairs;
        application.createModule<MContraction::WeakEye3pt>(name, wePar);
        results.push_back(name);
    };
    weakEye("we_diag", "loop", "diag");
    weakEye("we_all", "loop", "all");
    weakEye("we_pair", "loop", "(Identity GammaT)(GammaT Identity)");
    weakEye("we_diag_gT", "loop_gT", "diag");

    // save data
    MIO::WriteResultGroup::Par wPar;
    wPar.results = results;
    wPar.output  = "heavy_meson_3pt";
    application.createModule<MIO::WriteResultGroup>("WriteToFile", wPar);

    // execution
    application.saveParameterFile("HeavyMeson3pt.xml");
    application.run();

    // epilogue
    LOG(Message) << "Grid is finalizing now" << std::endl;
    Grid_finalize();

    return EXIT_SUCCESS;
}
