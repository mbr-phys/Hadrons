#include <Hadrons/Application.hpp>
#include <Hadrons/Modules.hpp>
#include "TestChecks.hpp"

using namespace Grid;
using namespace Hadrons;
using namespace HadronsTest;

/*
 * Validation of MContraction::HeavyMeson3pt (two-quark HQE operators),
 * MContraction::WeakNonEye3ptDeriv (four-quark operators with gamma pairs and
 * derivatives), MContraction::WeakEye3pt (gamma pairs) and the HQE insertion
 * modules (MUtilities::CovariantDerivative, MGauge::FieldStrength,
 * MUtilities::GaugeInsertion). Requires Nt > 4 (the sink is at t = 4).
 *
 * Point sources are used at both ends, so every comparison below is an
 * identity configuration by configuration (no noise average). Every checked
 * module writes its result to heavy_meson_3pt_check/<name>.1500.<ext>; after
 * application.run() the files are read back (TestChecks.hpp) and compared.
 * The test prints [PASS]/[FAIL]/[INFO] per check and returns EXIT_FAILURE if
 * any check fails. Tolerances (relative, norm over all times and momenta):
 *   exact identities (same propagators, only the contraction differs): 1e-10
 *   gauge covariance (propagators re-solved on the transformed field, so
 *   limited by the CG residual 1e-10):                                 1e-6
 *   vanishing on the unit gauge field (relative to C[Identity]):       1e-10
 *   approximate identities (O(a^2) or discretisation-dependent):  reported only
 * Gamma3pt results are corr[gamma][momentum][time], HeavyMeson3pt results
 * corr[triad][momentum][time], WeakNonEye3ptDeriv results
 * contractions[i].corr[momentum][time], WeakNonEye3pt and WeakEye3pt results
 * [i].corr[time]; entries are matched by their metadata, not by position.
 * Sigma<IJ> below is Grid's gamma_I gamma_J (no i), sigma_ij = i gamma_i gamma_j.
 *
 * Exact (failing) checks:
 *  1. 'hm3pt' vertices Gamma*            == 'ref_gamma3pt' (16 gammas)
 *     'hm3pt_mom' vertices Gamma*        == 'ref_gamma3pt_mom', all momenta
 *  2. 'hm3pt' D4                         == 'ref_D4' (CovariantDerivative)
 *     'hm3pt' D4                         == 0.5*'ref_rhq_D4' except at the
 *                                           operator times t = 0, Nt-1 (RHQ
 *                                           insertions use periodic links)
 *     'hm3pt' D4L                        == 'ref_D4L' (D_4 on q2 in Gamma3pt)
 *  3. 'hm3pt' D2, D4D4, D4D4D4, D2D4     == 'ref_D2', 'ref_D4D4', 'ref_D4D4D4',
 *                                           'ref_D2D4'
 *     'hm3pt_mom' at p = 0               == 'hm3pt' for every entry
 *  4. 'hm3pt' sigmaB, sigmaE, divE       == 'ref_sigmaB', 'ref_sigmaE', 'ref_divE'
 *     'hm3pt' sigmaBD4, sigmaED4         == 'ref_sigmaBD4', 'ref_sigmaED4'
 *     'hm3pt' DiE                        == sum_i 'ref_DiE_<i>' (D_i of E_i q3)
 *     'hm3pt' DiSigmaE                   == i sum_{i != j} 'ref_DiSigmaE_<i><j>'
 *  5. 'wne_DeltaF2', deriv 'none'        == 'ref_weakNonEye' (16 G x G,
 *                                           trace 1 and 2)
 *     'wne_DeltaF2', deriv 'qLeft@3'     == 'ref_wne_qLeft3'
 *     'wne_DeltaF2', deriv 'qBarRight@3' == 'ref_wne_qBarRight3'
 *     'wne_DeltaF2', specs A-D (same-bilinear pairs over 0123)
 *                                        == sum_rho 'ref_wne_{A,B,C,D}_<rho>'
 *     'wne_DeltaF2', deriv 'qLeft@0123:G2'
 *                                        == sum_rho 'ref_wne_cG2_<rho>'
 *                                           (legs none, qLeft -> D_rho qLeft,
 *                                           explicit pairs (G, gamma_rho G))
 *     'wne_mom_*' at p = 0               == 'wne_*' for every entry
 *  6. '*_g' (random gauge transformation) == the same instance on 'gauge',
 *     all entries (hm3pt_g, wne_g_DeltaF2, wne_g_DeltaF0_T1); tolerance 1e-6
 *  7. 'hm3pt_unit' (unit gauge): sigmaB, sigmaE, sigmaBD4, sigmaED4, DiE, divE,
 *     DiSigmaE, sigmaEDi, sigmaCurlE, gSigmaB, gSigmaE vanish
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
 * 13. per-leg gauge fields ('wne_legs', heavy legs on 'flowV', a stout
 *     smearing of the flowed links, light legs on the flowed links):
 *     deriv 'none', 'qLeft@3'            == the same entries of 'wne_DeltaF2'
 *     deriv 'qBarLeft@3', 'qRight@3'     == 'ref_wneV_qBarLeft3', 'ref_wneV_qRight3'
 *     specs A-D                          == sum_rho 'ref_wneV_{A,B,C,D}_<rho>'
 *                                           (heavy legs derived on 'flowV')
 * 11. spin-orbit vertices:
 *     'hm3pt' sigmaEDi                   == i sum_{i != j} 'ref_sigmaEDi_<i><j>'
 *     'hm3pt' sigmaCurlE                 == i sum_{i != j} 'ref_sigmaCurlE_<i><j>'
 *
 * Approximate identities (reported only):
 *  8. 'hm3pt' D4 + D4L ~ symmetric time derivative of C[Identity]
 * 11. 'hm3pt' DiSigmaE - sigmaEDi ~ sigmaCurlE (point-split, up to O(a^2))
 * 12. 'hm3pt' DiE ~ divE + sum_i 'ref_EiDi_<i>' (E_i D_i q3). In the continuum
 *     D_i(E_i q) = [D_i, E_i] q + E_i D_i q exactly; on the lattice DiE is
 *     point-split and divE is a colour-matrix insertion, so they agree only
 *     up to discretisation effects. Note: the physical (Hermitian) Darwin
 *     part satisfies Re<DiE> = 1/2 <divE> only after the ensemble average,
 *     never configuration by configuration.
 *
 * Not checkable on one configuration: hermiticity/reality (right- vs
 * left-acting insertions agree up to conjugation only in the ensemble
 * average), the Hermitian part 1/2 (DiSigmaE + sigmaEDi), and the physical
 * sign checks of lifetimes_round2.md section 7.
 */

static const std::string  checkDir = "heavy_meson_3pt_check";
static const unsigned int trajNum  = 1500;
static Checker            chk(checkDir, trajNum);

typedef MContraction::Gamma3pt::Result                    G3Result;
typedef MContraction::HeavyMeson3pt::Result               HMResult;
typedef MContraction::WeakNonEye3ptDeriv::Result          WNEDResult;
typedef std::vector<MContraction::WeakNonEye3pt::Result>  WNEResult;
typedef std::vector<MContraction::WeakEye3pt::Result>     WEResult;

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
    par.output       = chk.stem(name);
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
    par.output       = chk.stem(prefix + "_DeltaF2");
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
    par.output      = chk.stem(prefix + "_DeltaF0_T1");
    application.createModule<MContraction::WeakNonEye3ptDeriv>(prefix + "_DeltaF0_T1", par);

    return {prefix + "_DeltaF2", prefix + "_DeltaF0_T1"};
}

// result accessors (return an empty vector if the entry is missing, which the
// checks report as a failure) ////////////////////////////////////////////////
std::string gammaName(const Gamma::Algebra g)
{
    std::ostringstream s;

    s << g;

    return s.str();
}

// the 16 positive-sign gamma matrices, as expanded by "all"
std::vector<Gamma::Algebra> allGammas(void)
{
    std::vector<Gamma::Algebra> v;

    for (unsigned int i = 1; i < Gamma::nGamma; i += 2)
    {
        v.push_back(static_cast<Gamma::Algebra>(i));
    }

    return v;
}

// HeavyMeson3pt, sink = source = Gamma5, vertex v; mom < 0: all momenta
CVec hmCorr(const HMResult &r, const std::string &v, const int mom = -1)
{
    for (unsigned int i = 0; i < r.gamma.size(); ++i)
    {
        if ((r.gamma[i].vertex == v) and (r.gamma[i].sink == Gamma::Algebra::Gamma5)
            and (r.gamma[i].source == Gamma::Algebra::Gamma5))
        {
            return (mom < 0) ? flat(r.corr[i]) : r.corr[i][mom];
        }
    }

    return CVec();
}

// Gamma3pt, sink = source = Gamma5, vertex v
CVec g3Corr(const G3Result &r, const Gamma::Algebra v)
{
    for (unsigned int i = 0; i < r.gamma.size(); ++i)
    {
        if ((r.gamma[i].vertex == v) and (r.gamma[i].sink == Gamma::Algebra::Gamma5)
            and (r.gamma[i].source == Gamma::Algebra::Gamma5))
        {
            return flat(r.corr[i]);
        }
    }

    return CVec();
}

// WeakNonEye3ptDeriv entry; mom < 0: all momenta
CVec wnedCorr(const WNEDResult &r, const std::string &deriv, const Gamma::Algebra g1,
              const Gamma::Algebra g2, const unsigned int trace, const int mom = -1)
{
    for (auto &c: r.contractions)
    {
        if ((c.deriv == deriv) and (c.gamma1 == g1) and (c.gamma2 == g2) and (c.trace == trace))
        {
            return (mom < 0) ? flat(c.corr) : c.corr[mom];
        }
    }

    return CVec();
}

// WeakNonEye3pt entry
CVec wneCorr(const WNEResult &r, const Gamma::Algebra g, const unsigned int trace)
{
    for (auto &c: r)
    {
        if ((c.info.op == g) and (c.info.trace == trace))
        {
            return c.corr;
        }
    }

    return CVec();
}

// WeakEye3pt entry
CVec weCorr(const WEResult &r, const Gamma::Algebra g1, const Gamma::Algebra g2,
            const unsigned int trace)
{
    for (auto &c: r)
    {
        if ((c.info.gamma1 == g1) and (c.info.gamma2 == g2) and (c.info.trace == trace))
        {
            return c.corr;
        }
    }

    return CVec();
}

// index of the zero momentum, -1 if absent
int zeroMomentum(const std::vector<std::vector<int>> &mom)
{
    for (unsigned int m = 0; m < mom.size(); ++m)
    {
        bool zero = true;

        for (auto p: mom[m])
        {
            zero = zero and (p == 0);
        }
        if (zero)
        {
            return m;
        }
    }

    return -1;
}

// all entries of a HeavyMeson3pt / WeakNonEye3ptDeriv result, in result order
CVec allEntries(const HMResult &r)
{
    CVec f;

    for (auto &c: r.corr)
    {
        CVec x = flat(c);

        f.insert(f.end(), x.begin(), x.end());
    }

    return f;
}

CVec allEntries(const WNEDResult &r)
{
    CVec f;

    for (auto &c: r.contractions)
    {
        CVec x = flat(c.corr);

        f.insert(f.end(), x.begin(), x.end());
    }

    return f;
}

// times 1..nt-2 of a single-momentum correlator
CVec innerTimes(const CVec &c)
{
    return (c.size() > 2) ? CVec(c.begin() + 1, c.end() - 1) : CVec();
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
    const Complex            iUnit(0., 1.);
    const double             tolExact = 1.0e-10, tolGauge = 1.0e-6;

    // global parameters
    Application::GlobalPar globalPar;
    globalPar.trajCounter.start    = trajNum;
    globalPar.trajCounter.end      = trajNum + 20;
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

    // Gamma3pt reference: q3 (and optionally q2) replaced by a derived leg
    auto gamma3pt = [&](const std::string &name, const std::string &q3,
                        const std::string &gamma, const std::string &momProjector = "",
                        const std::string &q2 = "")
    {
        MContraction::Gamma3pt::Par g3Par;
        g3Par.q1           = "l_src_snk";
        g3Par.q2           = q2.empty() ? c.hSrc : q2;
        g3Par.q3           = q3;
        g3Par.gamma        = {"Gamma5 " + gamma + " Gamma5"};
        g3Par.tSnk         = tSnk;
        g3Par.momProjector = momProjector;
        g3Par.output       = chk.stem(name);
        application.createModule<MContraction::Gamma3pt>(name, g3Par);
        results.push_back(name);
    };

    MUtilities::CovariantDerivativePropagator::Par cdPar;
    cdPar.gauge = c.flowU;
    cdPar.bc    = -1;
    // symmetric derivative of a flowed field in direction rho (created once), on the
    // flowed links or on the given gauge field
    std::set<std::string> derived;
    auto derivLeg = [&](const std::string &q, const int rho, const std::string &gauge = "")
    {
        std::string name = q + "_D" + std::to_string(rho) + (gauge.empty() ? "" : "_" + gauge);

        if (derived.insert(name).second)
        {
            auto par = cdPar;

            par.gauge = gauge.empty() ? c.flowU : gauge;
            par.q     = q;
            par.type  = "symmetric";
            par.dirs  = std::to_string(rho);
            application.createModule<MUtilities::CovariantDerivativePropagator>(name, par);
        }

        return name;
    };
    // Laplacian of a flowed field over dirs
    auto laplaceLeg = [&](const std::string &name, const std::string &q, const std::string &dirs)
    {
        cdPar.q    = q;
        cdPar.type = "laplacian";
        cdPar.dirs = dirs;
        application.createModule<MUtilities::CovariantDerivativePropagator>(name, cdPar);

        return name;
    };

    MUtilities::GaugeInsertionPropagator::Par giPar;
    giPar.gauge = c.flowU;
    // gluonic insertion: mode sigmaB/sigmaE, or matrix with the colour matrix M
    auto insertion = [&](const std::string &name, const std::string &q,
                         const std::string &mode, const std::string &M,
                         const Gamma::Algebra gamma)
    {
        giPar.q      = q;
        giPar.mode   = mode;
        giPar.matrix = M;
        giPar.gamma  = gamma;
        application.createModule<MUtilities::GaugeInsertionPropagator>(name, giPar);

        return name;
    };

    MGauge::FieldStrength::Par fsPar;
    fsPar.gauge = c.flowU;

    // 1. plain gamma insertions, at zero and nonzero momentum
    gamma3pt("ref_gamma3pt", c.hSnk, "all");
    gamma3pt("ref_gamma3pt_mom", c.hSnk, "all", "phases");

    // 2. time derivative (new module and RHQ insertion = 2 D_4 with periodic links),
    //    and D_4 acting on qbar (derivative of q2 in Gamma3pt)
    const std::string hSnkD4 = derivLeg(c.hSnk, 3);
    gamma3pt("ref_D4", hSnkD4, "Identity");

    MRHQ::RHQInsertionI::Par rhqPar;
    rhqPar.q     = c.hSnk;
    rhqPar.index = "3";
    rhqPar.gamma = Gamma::Algebra::Identity;
    rhqPar.gauge = c.flowU;
    application.createModule<MRHQ::RHQInsertionI>("h_snk_rhqD4", rhqPar);
    gamma3pt("ref_rhq_D4", "h_snk_rhqD4", "Identity");
    gamma3pt("ref_D4L", c.hSnk, "Identity", "", derivLeg(c.hSrc, 3));

    // 3. Laplacians and nested time derivatives applied directly
    const std::string hSnkD2   = laplaceLeg("h_snk_D2", c.hSnk, "0 1 2");
    const std::string hSnkD4D4 = laplaceLeg("h_snk_D4D4", c.hSnk, "3");
    gamma3pt("ref_D2", hSnkD2, "Identity");
    gamma3pt("ref_D4D4", hSnkD4D4, "Identity");
    gamma3pt("ref_D4D4D4", derivLeg(hSnkD4D4, 3), "Identity");
    gamma3pt("ref_D2D4", laplaceLeg("h_snk_D2D4", hSnkD4, "0 1 2"), "Identity");

    // 4. gluonic insertions
    insertion("h_snk_sigmaB", c.hSnk, "sigmaB", "", Gamma::Algebra::Identity);
    insertion("h_snk_sigmaE", c.hSnk, "sigmaE", "", Gamma::Algebra::Identity);
    insertion("h_snk_sigmaBD4", hSnkD4, "sigmaB", "", Gamma::Algebra::Identity);
    insertion("h_snk_sigmaED4", hSnkD4, "sigmaE", "", Gamma::Algebra::Identity);
    gamma3pt("ref_sigmaB", "h_snk_sigmaB", "Identity");
    gamma3pt("ref_sigmaE", "h_snk_sigmaE", "Identity");
    gamma3pt("ref_sigmaBD4", "h_snk_sigmaBD4", "Identity");
    gamma3pt("ref_sigmaED4", "h_snk_sigmaED4", "Identity");

    fsPar.mode = "divE";
    fsPar.mu   = 0;
    fsPar.nu   = 0;
    application.createModule<MGauge::FieldStrength>("divE", fsPar);
    insertion("h_snk_divE", c.hSnk, "matrix", "divE", Gamma::Algebra::Identity);
    gamma3pt("ref_divE", "h_snk_divE", "Identity");

    // E_j = F_{j t}
    for (int j = 0; j < Tp; ++j)
    {
        fsPar.mode = "clover";
        fsPar.mu   = j;
        fsPar.nu   = Tp;
        application.createModule<MGauge::FieldStrength>("E_" + std::to_string(j), fsPar);
    }

    // DiE = sum_i D_i (E_i q3), and the reference for sum_i E_i D_i q3
    for (int i = 0; i < Tp; ++i)
    {
        const std::string si = std::to_string(i);

        insertion("h_snk_E_" + si, c.hSnk, "matrix", "E_" + si, Gamma::Algebra::Identity);
        gamma3pt("ref_DiE_" + si, derivLeg("h_snk_E_" + si, i), "Identity");
        insertion("h_snk_EiDi_" + si, derivLeg(c.hSnk, i), "matrix", "E_" + si,
                  Gamma::Algebra::Identity);
        gamma3pt("ref_EiDi_" + si, "h_snk_EiDi_" + si, "Identity");
    }

    // 9. chirality-even vertices: gamma_4 on the derived legs above, and
    //    gamma_i on D_i (of q3 or of D_4 q3), one module per i
    gamma3pt("ref_gD4", hSnkD4, "GammaT");
    gamma3pt("ref_gD2", hSnkD2, "GammaT");
    gamma3pt("ref_gSigmaB", "h_snk_sigmaB", "GammaT");
    gamma3pt("ref_gSigmaE", "h_snk_sigmaE", "GammaT");
    gamma3pt("ref_gD4D4", hSnkD4D4, "GammaT");
    const std::vector<std::string> gSpatial = {"GammaX", "GammaY", "GammaZ"};
    for (int i = 0; i < Tp; ++i)
    {
        gamma3pt("ref_gDi_" + std::to_string(i), derivLeg(c.hSnk, i), gSpatial[i]);
        gamma3pt("ref_gDiD4_" + std::to_string(i), derivLeg(hSnkD4, i), gSpatial[i]);
    }

    // 11. spin-orbit vertices, one module per (i, j), i != j:
    //     Sigma(i, j) E_j D_i q3, D_i (Sigma(i, j) E_j q3), Sigma(i, j) [D_i, E_j] q3
    auto sigmaName = [](const int i, const int j)
    {
        const std::string xyz = "XYZ";
        const int         lo = std::min(i, j), hi = std::max(i, j);

        return std::string(i < j ? "" : "Minus") + "Sigma" + xyz[lo] + xyz[hi];
    };
    for (int i = 0; i < Tp; ++i)
    for (int j = 0; j < Tp; ++j)
    {
        if (i == j)
        {
            continue;
        }

        const std::string    ij    = std::to_string(i) + std::to_string(j);
        const std::string    ej    = "E_" + std::to_string(j);
        const Gamma::Algebra sigma = strToVec<Gamma::Algebra>(sigmaName(i, j))[0];

        insertion("h_snk_sigmaEDi_" + ij, derivLeg(c.hSnk, i), "matrix", ej, sigma);
        gamma3pt("ref_sigmaEDi_" + ij, "h_snk_sigmaEDi_" + ij, "Identity");

        insertion("h_snk_sigmaE_" + ij, c.hSnk, "matrix", ej, sigma);
        gamma3pt("ref_DiSigmaE_" + ij, derivLeg("h_snk_sigmaE_" + ij, i), "Identity");

        fsPar.mode = "DE";
        fsPar.mu   = i;
        fsPar.nu   = j;
        application.createModule<MGauge::FieldStrength>("DE_" + ij, fsPar);
        insertion("h_snk_sigmaCurlE_" + ij, c.hSnk, "matrix", "DE_" + ij, sigma);
        gamma3pt("ref_sigmaCurlE_" + ij, "h_snk_sigmaCurlE_" + ij, "Identity");
    }

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
        wnePar.output    = chk.stem(name);
        application.createModule<MContraction::WeakNonEye3pt>(name, wnePar);
        results.push_back(name);
    };

    weakNonEye("ref_weakNonEye", c.hSrc, c.lSrcF, c.lSnkF, c.hSnk);
    weakNonEye("ref_wne_qLeft3", c.hSrc, derivLeg(c.lSrcF, 3), c.lSnkF, c.hSnk);
    weakNonEye("ref_wne_qBarRight3", c.hSrc, c.lSrcF, derivLeg(c.lSnkF, 3), c.hSnk);
    // specs A-D: derivatives on the two legs of one bilinear, one module per rho
    for (int rho = 0; rho < Nd; ++rho)
    {
        const std::string r = std::to_string(rho);

        weakNonEye("ref_wne_A_" + r, derivLeg(c.hSrc, rho), derivLeg(c.lSrcF, rho), c.lSnkF, c.hSnk);
        weakNonEye("ref_wne_B_" + r, c.hSrc, derivLeg(c.lSrcF, rho), c.lSnkF, derivLeg(c.hSnk, rho));
        weakNonEye("ref_wne_C_" + r, c.hSrc, c.lSrcF, derivLeg(c.lSnkF, rho), derivLeg(c.hSnk, rho));
        weakNonEye("ref_wne_D_" + r, derivLeg(c.hSrc, rho), c.lSrcF, derivLeg(c.lSnkF, rho), c.hSnk);
    }
    // 13. per-leg gauge fields: heavy legs on a stout smearing of the flowed links
    MGauge::StoutSmearing::Par stPar;
    stPar.gauge     = c.flowU;
    stPar.steps     = 1;
    stPar.orthogDim = "";
    stPar.rho       = 0.1;
    application.createModule<MGauge::StoutSmearing>("flowV", stPar);
    {
        MContraction::WeakNonEye3ptDeriv::Par par;

        par.qBarLeft     = c.hSrc;
        par.qLeft        = c.lSrcF;
        par.qBarRight    = c.lSnkF;
        par.qRight       = c.hSnk;
        par.gammaIn      = Gamma::Algebra::Gamma5;
        par.gammaOut     = Gamma::Algebra::Gamma5;
        par.gammaPairs   = "diag";
        par.derivatives  = {deriv("none", "", ""),
                            deriv("qLeft", "3", "none"),
                            deriv("qBarLeft", "3", "none"),
                            deriv("qRight", "3", "none"),
                            deriv("qBarLeft qLeft", "0 1 2 3", "none"),
                            deriv("qRight qLeft", "0 1 2 3", "none"),
                            deriv("qBarRight qRight", "0 1 2 3", "none"),
                            deriv("qBarLeft qBarRight", "0 1 2 3", "none")};
        par.gauge        = c.flowU;
        par.legGauges    = {"flowV", "", "", "flowV"};
        par.bc           = -1;
        par.momProjector = "";
        par.output       = chk.stem("wne_legs");
        application.createModule<MContraction::WeakNonEye3ptDeriv>("wne_legs", par);
        results.push_back("wne_legs");
    }
    weakNonEye("ref_wneV_qBarLeft3", derivLeg(c.hSrc, 3, "flowV"), c.lSrcF, c.lSnkF, c.hSnk);
    weakNonEye("ref_wneV_qRight3", c.hSrc, c.lSrcF, c.lSnkF, derivLeg(c.hSnk, 3, "flowV"));
    for (int rho = 0; rho < Nd; ++rho)
    {
        const std::string r = std::to_string(rho);

        weakNonEye("ref_wneV_A_" + r, derivLeg(c.hSrc, rho, "flowV"), derivLeg(c.lSrcF, rho), c.lSnkF, c.hSnk);
        weakNonEye("ref_wneV_B_" + r, c.hSrc, derivLeg(c.lSrcF, rho), c.lSnkF, derivLeg(c.hSnk, rho, "flowV"));
        weakNonEye("ref_wneV_C_" + r, c.hSrc, c.lSrcF, derivLeg(c.lSnkF, rho), derivLeg(c.hSnk, rho, "flowV"));
        weakNonEye("ref_wneV_D_" + r, derivLeg(c.hSrc, rho, "flowV"), c.lSrcF, derivLeg(c.lSnkF, rho), c.hSnk);
    }

    // coupled spec qLeft@0123:G2: underived WeakNonEye3ptDeriv with qLeft -> D_rho qLeft
    // and the explicit pairs (G, gamma_rho G), one module per rho
    for (int rho = 0; rho < Nd; ++rho)
    {
        const std::string                 name = "ref_wne_cG2_" + std::to_string(rho);
        std::string                       pairs;
        MContraction::WeakNonEye3ptDeriv::Par par;

        for (auto &G: Gamma::gall)
        {
            pairs += "(" + gammaName(G.g) + " " + gammaName((Gamma::gmu[rho]*G).g) + ")";
        }
        par.qBarLeft     = c.hSrc;
        par.qLeft        = derivLeg(c.lSrcF, rho);
        par.qBarRight    = c.lSnkF;
        par.qRight       = c.hSnk;
        par.gammaIn      = Gamma::Algebra::Gamma5;
        par.gammaOut     = Gamma::Algebra::Gamma5;
        par.gammaPairs   = pairs;
        par.derivatives  = {deriv("none", "", "")};
        par.gauge        = c.flowU;
        par.bc           = -1;
        par.momProjector = "";
        par.output       = chk.stem(name);
        application.createModule<MContraction::WeakNonEye3ptDeriv>(name, par);
        results.push_back(name);
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
        wePar.output     = chk.stem(name);
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

    // checks //////////////////////////////////////////////////////////////////
    const unsigned int nt = GridDefaultLatt()[Tp];
    auto readG3 = [](const std::string &name)
    {
        return chk.read<G3Result>(name, "gamma3pt");
    };
    // single-vertex Gamma3pt reference (Identity or the given vertex)
    auto ref = [&](const std::string &name, const Gamma::Algebra v = Gamma::Algebra::Identity)
    {
        return g3Corr(readG3(name), v);
    };

    HMResult hm     = chk.read<HMResult>("hm3pt", "heavyMeson3pt");
    HMResult hmMom  = chk.read<HMResult>("hm3pt_mom", "heavyMeson3pt");
    HMResult hmG    = chk.read<HMResult>("hm3pt_g", "heavyMeson3pt");
    HMResult hmUnit = chk.read<HMResult>("hm3pt_unit", "heavyMeson3pt");

    // 1. gamma vertices
    {
        G3Result r    = readG3("ref_gamma3pt");
        G3Result rMom = readG3("ref_gamma3pt_mom");

        if (hmMom.momentum != rMom.momentum)
        {
            chk.fail("1. hm3pt_mom vs ref_gamma3pt_mom", "momentum lists differ");
        }
        for (auto g: allGammas())
        {
            chk.exact("1. hm3pt " + gammaName(g) + " vs ref_gamma3pt",
                      hmCorr(hm, gammaName(g)), g3Corr(r, g), tolExact);
            chk.exact("1. hm3pt_mom " + gammaName(g) + " vs ref_gamma3pt_mom",
                      hmCorr(hmMom, gammaName(g)), g3Corr(rMom, g), tolExact);
        }
    }

    // 2.-4. derivative and gluonic vertices
    {
        CVec half;

        axpy(half, 0.5, ref("ref_rhq_D4"));
        chk.exact("2. D4 vs ref_D4", hmCorr(hm, "D4"), ref("ref_D4"), tolExact);
        chk.exact("2. D4 vs 0.5 ref_rhq_D4 (t = 1..Nt-2)", innerTimes(hmCorr(hm, "D4")),
                  innerTimes(half), tolExact);
        chk.exact("2. D4L vs ref_D4L", hmCorr(hm, "D4L"), ref("ref_D4L"), tolExact);
    }
    const std::vector<std::pair<std::string, std::string>> direct = {
        {"D2", "ref_D2"}, {"D4D4", "ref_D4D4"}, {"D4D4D4", "ref_D4D4D4"},
        {"D2D4", "ref_D2D4"}, {"sigmaB", "ref_sigmaB"}, {"sigmaE", "ref_sigmaE"},
        {"sigmaBD4", "ref_sigmaBD4"}, {"sigmaED4", "ref_sigmaED4"}, {"divE", "ref_divE"}};
    for (auto &d: direct)
    {
        chk.exact("3./4. " + d.first + " vs " + d.second, hmCorr(hm, d.first),
                  ref(d.second), tolExact);
    }
    {
        CVec diE, eiDi, diSigmaE, sigmaEDi, sigmaCurlE, gDi, gDiD4;

        for (int i = 0; i < Tp; ++i)
        {
            const std::string si = std::to_string(i);

            axpy(diE, 1., ref("ref_DiE_" + si));
            axpy(eiDi, 1., ref("ref_EiDi_" + si));
            axpy(gDi, 1., ref("ref_gDi_" + si, strToVec<Gamma::Algebra>(gSpatial[i])[0]));
            axpy(gDiD4, 1., ref("ref_gDiD4_" + si, strToVec<Gamma::Algebra>(gSpatial[i])[0]));
            for (int j = 0; j < Tp; ++j)
            {
                if (i == j)
                {
                    continue;
                }

                const std::string ij = si + std::to_string(j);

                axpy(diSigmaE, iUnit, ref("ref_DiSigmaE_" + ij));
                axpy(sigmaEDi, iUnit, ref("ref_sigmaEDi_" + ij));
                axpy(sigmaCurlE, iUnit, ref("ref_sigmaCurlE_" + ij));
            }
        }
        chk.exact("4. DiE vs sum_i ref_DiE_i", hmCorr(hm, "DiE"), diE, tolExact);
        chk.exact("4. DiSigmaE vs i sum_ij ref_DiSigmaE_ij", hmCorr(hm, "DiSigmaE"), diSigmaE, tolExact);
        chk.exact("11. sigmaEDi vs i sum_ij ref_sigmaEDi_ij", hmCorr(hm, "sigmaEDi"), sigmaEDi, tolExact);
        chk.exact("11. sigmaCurlE vs i sum_ij ref_sigmaCurlE_ij", hmCorr(hm, "sigmaCurlE"), sigmaCurlE, tolExact);
        chk.exact("9. gDi vs sum_i ref_gDi_i", hmCorr(hm, "gDi"), gDi, tolExact);
        chk.exact("9. gDiD4 vs sum_i ref_gDiD4_i", hmCorr(hm, "gDiD4"), gDiD4, tolExact);

        // 12. Darwin: DiE ~ divE + sum_i E_i D_i (approximate on the lattice;
        //     Re<DiE> = 1/2 <divE> only holds after the ensemble average)
        CVec darwin = hmCorr(hm, "divE");

        axpy(darwin, 1., eiDi);
        chk.report("12. DiE vs divE + sum_i E_i D_i", hmCorr(hm, "DiE"), darwin);

        // 11. DiSigmaE - sigmaEDi ~ sigmaCurlE (point split)
        CVec split = hmCorr(hm, "DiSigmaE");

        axpy(split, -1., hmCorr(hm, "sigmaEDi"));
        chk.report("11. DiSigmaE - sigmaEDi vs sigmaCurlE", split, hmCorr(hm, "sigmaCurlE"));
    }

    // 9. chirality-even vertices with a single reference
    chk.exact("9. gD4 vs ref_gD4", hmCorr(hm, "gD4"), ref("ref_gD4", Gamma::Algebra::GammaT), tolExact);
    chk.exact("9. gD2 vs ref_gD2", hmCorr(hm, "gD2"), ref("ref_gD2", Gamma::Algebra::GammaT), tolExact);
    chk.exact("9. gD4D4 vs ref_gD4D4", hmCorr(hm, "gD4D4"), ref("ref_gD4D4", Gamma::Algebra::GammaT), tolExact);
    chk.exact("9. gSigmaB vs ref_gSigmaB", hmCorr(hm, "gSigmaB"), ref("ref_gSigmaB", Gamma::Algebra::GammaT), tolExact);
    chk.exact("9. gSigmaE vs ref_gSigmaE", hmCorr(hm, "gSigmaE"), ref("ref_gSigmaE", Gamma::Algebra::GammaT), tolExact);

    // 3. momentum projection at p = 0 vs sliceSum, every entry
    {
        const int p0 = zeroMomentum(hmMom.momentum);

        if (p0 < 0)
        {
            chk.fail("3. hm3pt_mom p = 0", "no zero momentum in the projector list");
        }
        else
        {
            for (auto &t: hm.gamma)
            {
                if ((t.sink == Gamma::Algebra::Gamma5) and (t.source == Gamma::Algebra::Gamma5))
                {
                    chk.exact("3. hm3pt_mom(p = 0) vs hm3pt " + t.vertex,
                              hmCorr(hmMom, t.vertex, p0), hmCorr(hm, t.vertex), tolExact);
                }
            }
        }
    }

    // 8. D4 + D4L ~ symmetric time derivative of C[Identity] (approximate)
    {
        CVec c1 = hmCorr(hm, "Identity"), dt, sum;

        for (unsigned int t = 1; (t + 1 < nt) and (c1.size() == nt); ++t)
        {
            dt.push_back(0.5*(c1[t + 1] - c1[t - 1]));
        }
        axpy(sum, 1., innerTimes(hmCorr(hm, "D4")));
        axpy(sum, 1., innerTimes(hmCorr(hm, "D4L")));
        chk.report("8. D4 + D4L vs (C[1](t+1) - C[1](t-1))/2", sum, dt);
    }

    // 6. gauge covariance, all entries
    chk.exact("6. hm3pt_g vs hm3pt (all entries)", allEntries(hmG), allEntries(hm), tolGauge);

    // 7. unit gauge: field-strength vertices vanish
    {
        const double scale = norm(hmCorr(hmUnit, "Identity"));

        for (std::string v: {"sigmaB", "sigmaE", "sigmaBD4", "sigmaED4", "DiE", "divE",
                             "DiSigmaE", "sigmaEDi", "sigmaCurlE", "gSigmaB", "gSigmaE"})
        {
            chk.vanish("7. hm3pt_unit " + v + " vanishes", hmCorr(hmUnit, v), scale, tolExact);
        }
    }

    // 5. four-quark operators
    {
        WNEDResult wne     = chk.read<WNEDResult>("wne_DeltaF2", "weakNonEye3ptDeriv");
        WNEDResult wneT1   = chk.read<WNEDResult>("wne_DeltaF0_T1", "weakNonEye3ptDeriv");
        WNEDResult wneMom  = chk.read<WNEDResult>("wne_mom_DeltaF2", "weakNonEye3ptDeriv");
        WNEDResult wneMomT1 = chk.read<WNEDResult>("wne_mom_DeltaF0_T1", "weakNonEye3ptDeriv");
        WNEDResult wneG    = chk.read<WNEDResult>("wne_g_DeltaF2", "weakNonEye3ptDeriv");
        WNEDResult wneGT1  = chk.read<WNEDResult>("wne_g_DeltaF0_T1", "weakNonEye3ptDeriv");
        auto readWNE = [](const std::string &name)
        {
            return chk.read<WNEResult>(name, "weakNonEye3pt");
        };
        WNEResult  refNone = readWNE("ref_weakNonEye");
        WNEResult  refQL3  = readWNE("ref_wne_qLeft3");
        WNEResult  refQBR3 = readWNE("ref_wne_qBarRight3");
        const std::vector<std::pair<std::string, std::string>> specs = {
            {"A", "qBarLeft qLeft@0123"}, {"B", "qRight qLeft@0123"},
            {"C", "qBarRight qRight@0123"}, {"D", "qBarLeft qBarRight@0123"}};
        std::map<std::string, std::vector<WNEResult>> refSpec;
        std::vector<WNEDResult>                       refCG2;

        for (int rho = 0; rho < Nd; ++rho)
        {
            for (auto &sp: specs)
            {
                refSpec[sp.first].push_back(readWNE("ref_wne_" + sp.first + "_" + std::to_string(rho)));
            }
            refCG2.push_back(chk.read<WNEDResult>("ref_wne_cG2_" + std::to_string(rho),
                                                  "weakNonEye3ptDeriv"));
        }
        for (auto &G: Gamma::gall)
        for (unsigned int tr = 1; tr <= 2; ++tr)
        {
            const Gamma::Algebra g   = G.g;
            const std::string    tag = " " + gammaName(g) + " trace " + std::to_string(tr);

            chk.exact("5. wne none vs ref_weakNonEye" + tag,
                      wnedCorr(wne, "none", g, g, tr), wneCorr(refNone, g, tr), tolExact);
            chk.exact("5. wne qLeft@3 vs ref_wne_qLeft3" + tag,
                      wnedCorr(wne, "qLeft@3", g, g, tr), wneCorr(refQL3, g, tr), tolExact);
            chk.exact("5. wne qBarRight@3 vs ref_wne_qBarRight3" + tag,
                      wnedCorr(wne, "qBarRight@3", g, g, tr), wneCorr(refQBR3, g, tr), tolExact);
            for (auto &sp: specs)
            {
                CVec sum;

                for (int rho = 0; rho < Nd; ++rho)
                {
                    axpy(sum, 1., wneCorr(refSpec[sp.first][rho], g, tr));
                }
                chk.exact("5. wne " + sp.second + " vs sum_rho ref_wne_" + sp.first + tag,
                          wnedCorr(wne, sp.second, g, g, tr), sum, tolExact);
            }

            CVec sumCG2;

            for (int rho = 0; rho < Nd; ++rho)
            {
                const Gamma::Algebra g2 = (Gamma::gmu[rho]*G).g;

                axpy(sumCG2, 1., wnedCorr(refCG2[rho], "none", g, g2, tr));
            }
            chk.exact("5. wne qLeft@0123:G2 vs sum_rho ref_wne_cG2" + tag,
                      wnedCorr(wne, "qLeft@0123:G2", g, g, tr), sumCG2, tolExact);
        }

        // momentum projection at p = 0, every entry
        const int p0 = zeroMomentum(wneMom.momentum);

        if (p0 < 0)
        {
            chk.fail("5. wne_mom p = 0", "no zero momentum in the projector list");
        }
        else
        {
            std::vector<std::pair<const WNEDResult *, const WNEDResult *>> momPairs =
                {{&wneMom, &wne}, {&wneMomT1, &wneT1}};

            for (auto &mp: momPairs)
            {
                CVec a, b;

                for (auto &ct: mp.first->contractions)
                {
                    a.insert(a.end(), ct.corr[p0].begin(), ct.corr[p0].end());
                }
                b = allEntries(*mp.second);
                chk.exact("5. wne_mom(p = 0) vs wne (all entries)", a, b, tolExact);
            }
        }

        // 13. per-leg gauge fields
        {
            WNEDResult wneLegs = chk.read<WNEDResult>("wne_legs", "weakNonEye3ptDeriv");
            WNEResult  refQBL3 = readWNE("ref_wneV_qBarLeft3");
            WNEResult  refQR3  = readWNE("ref_wneV_qRight3");
            std::map<std::string, std::vector<WNEResult>> refV;

            for (int rho = 0; rho < Nd; ++rho)
            {
                for (auto &sp: specs)
                {
                    refV[sp.first].push_back(readWNE("ref_wneV_" + sp.first + "_" + std::to_string(rho)));
                }
            }
            for (auto &G: Gamma::gall)
            for (unsigned int tr = 1; tr <= 2; ++tr)
            {
                const Gamma::Algebra g   = G.g;
                const std::string    tag = " " + gammaName(g) + " trace " + std::to_string(tr);

                chk.exact("13. wne_legs none vs wne" + tag,
                          wnedCorr(wneLegs, "none", g, g, tr), wnedCorr(wne, "none", g, g, tr), tolExact);
                chk.exact("13. wne_legs qLeft@3 vs wne" + tag,
                          wnedCorr(wneLegs, "qLeft@3", g, g, tr), wnedCorr(wne, "qLeft@3", g, g, tr), tolExact);
                chk.exact("13. wne_legs qBarLeft@3 vs ref_wneV_qBarLeft3" + tag,
                          wnedCorr(wneLegs, "qBarLeft@3", g, g, tr), wneCorr(refQBL3, g, tr), tolExact);
                chk.exact("13. wne_legs qRight@3 vs ref_wneV_qRight3" + tag,
                          wnedCorr(wneLegs, "qRight@3", g, g, tr), wneCorr(refQR3, g, tr), tolExact);
                for (auto &sp: specs)
                {
                    CVec sum;

                    for (int rho = 0; rho < Nd; ++rho)
                    {
                        axpy(sum, 1., wneCorr(refV[sp.first][rho], g, tr));
                    }
                    chk.exact("13. wne_legs " + sp.second + " vs sum_rho ref_wneV_" + sp.first + tag,
                              wnedCorr(wneLegs, sp.second, g, g, tr), sum, tolExact);
                }
            }
        }

        // 6. gauge covariance
        chk.exact("6. wne_g_DeltaF2 vs wne_DeltaF2 (all entries)", allEntries(wneG),
                  allEntries(wne), tolGauge);
        chk.exact("6. wne_g_DeltaF0_T1 vs wne_DeltaF0_T1 (all entries)", allEntries(wneGT1),
                  allEntries(wneT1), tolGauge);
    }

    // 10. eye gamma pairs
    {
        auto readWE = [](const std::string &name)
        {
            return chk.read<WEResult>(name, "weakEye3pt");
        };
        WEResult             diag   = readWE("we_diag");
        WEResult             all    = readWE("we_all");
        WEResult             pair   = readWE("we_pair");
        WEResult             diagGT = readWE("we_diag_gT");
        const Gamma::Algebra one = Gamma::Algebra::Identity, gT = Gamma::Algebra::GammaT;

        for (auto &G: Gamma::gall)
        for (unsigned int tr = 1; tr <= 2; ++tr)
        {
            chk.exact("10. we_all vs we_diag " + gammaName(G.g) + " trace " + std::to_string(tr),
                      weCorr(all, G.g, G.g, tr), weCorr(diag, G.g, G.g, tr), tolExact);
        }
        for (unsigned int tr = 1; tr <= 2; ++tr)
        {
            const std::string t = " trace " + std::to_string(tr);

            chk.exact("10. we_pair (Identity GammaT) vs we_diag_gT Identity" + t,
                      weCorr(pair, one, gT, tr), weCorr(diagGT, one, one, tr), tolExact);
            chk.exact("10. we_pair (GammaT Identity) vs we_diag_gT GammaT" + t,
                      weCorr(pair, gT, one, tr), weCorr(diagGT, gT, gT, tr), tolExact);
        }
    }

    const int status = chk.summary();

    // epilogue
    LOG(Message) << "Grid is finalizing now" << std::endl;
    Grid_finalize();

    return status;
}
