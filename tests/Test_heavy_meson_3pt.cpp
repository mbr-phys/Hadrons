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
 * Validation workflow for MContraction::HeavyMeson3pt and the HQE insertion
 * modules (MUtilities::CovariantDerivative, MGauge::FieldStrength,
 * MUtilities::GaugeInsertion). Requires Nt > 4 (the sink is at t = 4).
 *
 * Point sources are used at both ends, so every comparison below is an exact
 * identity configuration by configuration (no noise average):
 *
 *  1. 'hm3pt' hqe 'gamma' entries        == 'ref_gamma3pt'
 *  2. 'hm3pt' hqe D4                     == 'ref_D4'
 *     'hm3pt' hqe D4                     == 0.5*'ref_rhq_D4' except at the
 *                                           operator times t = 0, Nt-1 (RHQ
 *                                           insertions use periodic links)
 *  3. 'hm3pt' hqe D2                     == 'ref_D2' (tests the summation by parts)
 *  4. 'hm3pt' hqe sigmaB, sigmaE, divE   == 'ref_sigmaB', 'ref_sigmaE', 'ref_divE'
 *  5. 'hm3pt' fourQuark 'DeltaF2', deriv 'none' == 'ref_weakNonEye' (16 G x G,
 *                                           trace 1 and 2)
 *  6. 'hm3pt_g' (random gauge transformation) == 'hm3pt', all entries
 *  7. 'hm3pt_unit' (unit gauge): sigmaB, sigmaE, sigmaBD4, sigmaED4, DiE, divE,
 *     DiSigmaE vanish
 *  8. 'hm3pt': D4 + D4L should be close to the time derivative of the
 *     'gamma' Identity correlator; hermiticity/reality checks as in
 *     lifetimes_round2.md §7
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

// HeavyMeson3pt with all HQE operators and a set of four-quark contractions
std::string makeHM3pt(Application &application, const Chain &c,
                      const std::string &name)
{
    MContraction::HeavyMeson3pt::Par par;
    par.gauge    = c.flowU;
    par.bc       = -1;
    par.flowTime = 0.02;

    par.hqe.qSrc      = c.hSrc;
    par.hqe.qSnk      = c.hSnk;
    par.hqe.spectator = c.lSrc;
    par.hqe.sinkNoise = c.ptSnk;
    par.hqe.tSnk      = tSnk;
    par.hqe.gammas    = "all";
    par.hqe.ops       = "all";

    MContraction::HeavyMeson3ptFourQuarkPar dB2;
    dB2.name      = "DeltaF2";
    dB2.qBarLeft  = c.hSrc;
    dB2.qLeft     = c.lSrcF;
    dB2.qBarRight = c.lSnkF;
    dB2.qRight    = c.hSnk;

    auto deriv = [](const std::string &legs, const std::string &dirs,
                    const std::string &couple)
    {
        MContraction::HeavyMeson3ptDerivPar d;
        d.legs   = legs;
        d.dirs   = dirs;
        d.couple = couple;

        return d;
    };
    dB2.derivatives = {deriv("", "", ""),
                       // leading power: time derivative on the light legs
                       deriv("qLeft", "3", ""),
                       deriv("qBarRight", "3", ""),
                       // full R2/R3-type: derivative on qbar and q of one bilinear
                       deriv("qBarLeft qLeft", "0 1 2 3", ""),
                       deriv("qRight qLeft", "0 1 2 3", ""),
                       deriv("qBarRight qRight", "0 1 2 3", ""),
                       deriv("qBarLeft qBarRight", "0 1 2 3", ""),
                       // index coupled to the second gamma
                       deriv("qLeft", "0 1 2 3", "G2")};

    MContraction::HeavyMeson3ptFourQuarkPar dB0;
    dB0.name       = "DeltaF0_T1";
    dB0.qBarLeft   = c.lSrcF;
    dB0.qLeft      = c.hSrc;
    dB0.qBarRight  = c.lSnkF;
    dB0.qRight     = c.hSnk;
    dB0.gammaPairs = "(Identity GammaT)(GammaT Identity)";

    par.fourQuark = {dB2, dB0};
    application.createModule<MContraction::HeavyMeson3pt>(name, par);

    return name;
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

    results.push_back(makeHM3pt(application, c, "hm3pt"));
    results.push_back(makeHM3pt(application, cg, "hm3pt_g"));
    results.push_back(makeHM3pt(application, cu, "hm3pt_unit"));

    // references on the random gauge field ////////////////////////////////////
    // spectator closed at the sink point
    MSink::SourceSink::Par ssPar;
    ssPar.source = c.ptSnk;
    application.createModule<MSink::SourceSink>("snk_pt", ssPar);
    MSink::Smear::Par smPar;
    smPar.q    = c.lSrc;
    smPar.sink = "snk_pt";
    application.createModule<MSink::Smear>("l_src_snk", smPar);

    auto gamma3pt = [&](const std::string &name, const std::string &q3,
                        const std::string &gamma)
    {
        MContraction::Gamma3pt::Par g3Par;
        g3Par.q1    = "l_src_snk";
        g3Par.q2    = c.hSrc;
        g3Par.q3    = q3;
        g3Par.gamma = gamma;
        g3Par.tSnk  = tSnk;
        application.createModule<MContraction::Gamma3pt>(name, g3Par);
        results.push_back(name);
    };

    // 1. plain gamma insertions
    gamma3pt("ref_gamma3pt", c.hSnk, "all");

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
    fsPar.mode  = "divergence";
    fsPar.nu    = 3;
    fsPar.dirs  = "0 1 2";
    application.createModule<MGauge::FieldStrength>("divE", fsPar);
    giPar.mode   = "matrix";
    giPar.matrix = "divE";
    application.createModule<MUtilities::GaugeInsertionPropagator>("h_snk_divE", giPar);
    gamma3pt("ref_divE", "h_snk_divE", "Identity");

    // 5. four-quark reference
    MContraction::WeakNonEye3pt::Par wnePar;
    wnePar.qBarLeft  = c.hSrc;
    wnePar.qLeft     = c.lSrcF;
    wnePar.qBarRight = c.lSnkF;
    wnePar.qRight    = c.hSnk;
    wnePar.gammaIn   = Gamma::Algebra::Gamma5;
    wnePar.gammaOut  = Gamma::Algebra::Gamma5;
    application.createModule<MContraction::WeakNonEye3pt>("ref_weakNonEye", wnePar);
    results.push_back("ref_weakNonEye");

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
