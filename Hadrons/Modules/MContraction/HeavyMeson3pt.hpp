/*
 * HeavyMeson3pt.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
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
#ifndef Hadrons_MContraction_HeavyMeson3pt_hpp_
#define Hadrons_MContraction_HeavyMeson3pt_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Serialization.hpp>
#include <Hadrons/Modules/MUtilities/DerivHelpers.hpp>
#include <Hadrons/Modules/MUtilities/HQEOperators.hpp>
#include <Grid/algorithms/blas/MomentumProject.h>

BEGIN_HADRONS_NAMESPACE

/*
 * Heavy-meson three-point functions with HQE two-quark operator insertions,
 * for one (source wall, sink wall) pair at one flow time. Same contraction
 * and gamma conventions as MContraction::Gamma3pt. Four-quark operators are
 * in MContraction::WeakNonEye3ptDeriv.
 *
 * All operator fields (heavy propagators at the operator, gauge links, field
 * strength) are taken at the same flow time; the spectator line is unflowed.
 * Conventions for derivatives are those of MUtilities/DerivHelpers.hpp, for
 * field strength and sigma_{mu nu} those of MUtilities/HQEOperators.hpp.
 *
 *                   q2           q3
 *              /----<------*------<----¬
 *             /            X            \
 *            /                           \
 *   i, gSrc *                            * f, gSnk
 *            \                          /
 *             \                        /
 *              \----------->----------/
 *                       spectator
 *
 *   C[gSnk, X, gSrc](t) = sum_x tr[ gSnk q1Snk (adj(gSrc) g5) adj(q2)(x) (g5 X) q3(x) ]
 *   q1Snk = sum_{z, z_t = tSnk} adj(sinkNoise)(z) spectator(z)
 *   (for a plain wall sum, pass a wall source as sinkNoise)
 *
 *   q2 = qSrc (heavy from the source wall), q3 = qSnk (heavy from the sink
 *   wall), spectator = light from the source wall (unflowed). For a gamma
 *   vertex X = gV this is exactly MContraction::Gamma3pt with q1 the
 *   spectator sinked with adj(sinkNoise).
 *
 *   gamma: vector of "sink vertex source" triads, as in Gamma3pt. Sink and
 *   source are gamma matrices or "all"; the vertex is
 *     - a gamma matrix gV or "all" (X = gV, as in Gamma3pt), or
 *     - an HQE operator O (X = O), or "ops" for all of them:
 *       D2       sum_{i<3} D^-_i D^+_i              nearest-neighbour Laplacian
 *       D4       D_4                                symmetric
 *       D4L      left-acting D_4 on qbar            (D4 + D4L ~ d/dt C[Identity])
 *       D4D4     D^-_4 D^+_4
 *       D4D4D4   D_4 D^-_4 D^+_4
 *       D2D4     sum_i D^-_i D^+_i D_4
 *       sigmaB   sum_{i<j<3} sigma_ij F_ij
 *       sigmaE   sum_{j<3} sigma_j3 F_j3
 *       sigmaBD4 sigmaB D_4
 *       sigmaED4 sigmaE D_4
 *       DiE      sum_i D_i F_i3                     (derivative acts on F q)
 *       divE     sum_i [D_i, F_i3]                  (colour-matrix insertion)
 *       DiSigmaE sum_{i,j} D_i sigma_ij F_j3        (derivative acts on F q)
 *     chirality-even partners of the operators above, needed for the
 *     positive-flow-time mixing (gf_lifetime_notes, operator basis):
 *       gD4      gamma_4 D_4                        (twist-2 T_44)
 *       gDi      sum_{i<3} gamma_i D_i              (with gD4: Dslash)
 *       gD4D4    gamma_4 D^-_4 D^+_4
 *       gD2      gamma_4 sum_{i<3} D^-_i D^+_i
 *       gDiD4    sum_{i<3} gamma_i D_i D_4
 *       gSigmaB  gamma_4 sigmaB
 *       gSigmaE  gamma_4 sigmaE
 *   Repeated triads are removed. Every operator except D4L is applied
 *   directly to q3.
 *
 * Output:
 *   momentum: list of momenta ({0, 0, 0} without momProjector)
 *   gamma:    list of {sink, vertex, source} triads
 *   corr:     [triad][momentum][time]
 *
 * Parameters:
 * - gauge:        flowed gauge field (links and field strength)
 * - bc:           fermion temporal boundary phase for time derivatives (+1, -1)
 * - flowTime:     label stored with the result
 * - momProjector: MomentumProject module (e.g. MUtilities::MPScalar); empty for
 *                 a zero-momentum sum only (as in Gamma3pt)
 * - qSrc:         flowed heavy propagator from the source wall
 * - qSnk:         flowed heavy propagator from the sink wall
 * - spectator:    unflowed light propagator from the source wall
 * - sinkNoise:    source field of the sink wall
 * - tSnk:         sink time
 * - gamma:        "sink vertex source" triads
 * - output:       result file stem (empty = no file)
 *
 * There are no default values: every parameter except momProjector and
 * output must be set.
 */

/******************************************************************************
 *                             HeavyMeson3pt                                  *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MContraction)

class HeavyMeson3ptPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(HeavyMeson3ptPar,
                                    std::string,              gauge,
                                    int,                      bc,
                                    double,                   flowTime,
                                    std::string,              momProjector,
                                    std::string,              qSrc,
                                    std::string,              qSnk,
                                    std::string,              spectator,
                                    std::string,              sinkNoise,
                                    unsigned int,             tSnk,
                                    std::vector<std::string>, gamma,
                                    std::string,              output);
};

template <typename FImpl, typename GImpl>
class THeavyMeson3pt: public Module<HeavyMeson3ptPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    typedef typename GImpl::GaugeLinkField GaugeLinkField;
    using MPType = MomentumProject<ComplexField, ComplexField>;
    typedef std::vector<std::vector<Complex>> MomCorr;
    class HQETriad: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(HQETriad,
                                        Gamma::Algebra, sink,
                                        std::string,    vertex,
                                        Gamma::Algebra, source);
    };
    class Result: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(Result,
                                        double,                        flowTime,
                                        unsigned int,                  tSnk,
                                        std::vector<std::vector<int>>, momentum,
                                        std::vector<HQETriad>,         gamma,
                                        std::vector<MomCorr>,          corr);
    };
public:
    // constructor
    THeavyMeson3pt(const std::string name);
    // destructor
    virtual ~THeavyMeson3pt(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    virtual std::vector<std::string> getOutputFiles(void);
protected:
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    void computeHQE(Result &res, const GaugeField &U);
    void parseHQEGammas(void);
    bool usedPair(const std::string &v, const Gamma::Algebra src) const;
    void project(MomCorr &acc, ComplexField &c);
private:
    std::vector<HQETriad>       triads_;
    std::vector<Gamma::Algebra> sinks_, sources_;
    std::vector<std::string>    vertices_;
    bool                        zeroMomentumOnly_, needF_;
    unsigned int                nMom_;
};

MODULE_REGISTER_TMP(HeavyMeson3pt, ARG(THeavyMeson3pt<FIMPL, GIMPL>), MContraction);

/******************************************************************************
 *                     THeavyMeson3pt implementation                          *
 ******************************************************************************/
static const std::vector<std::string> heavyMeson3ptOps = {"D2", "D4", "D4L", "D4D4", "D4D4D4", "D2D4", "sigmaB", "sigmaE", "sigmaBD4", "sigmaED4", "DiE", "divE", "DiSigmaE",
                                                          "gD4", "gDi", "gD4D4", "gD2", "gDiD4", "gSigmaB", "gSigmaE"};

// operators involving the field strength
static const std::vector<std::string> heavyMeson3ptFOps = {"sigmaB", "sigmaE", "sigmaBD4", "sigmaED4", "DiE", "divE", "DiSigmaE", "gSigmaB", "gSigmaE"};

template <typename T>
inline bool heavyMeson3ptContains(const std::vector<T> &v, const T &x)
{
    return std::find(v.begin(), v.end(), x) != v.end();
}

// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl, typename GImpl>
THeavyMeson3pt<FImpl, GImpl>::THeavyMeson3pt(const std::string name)
: Module<HeavyMeson3ptPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl, typename GImpl>
std::vector<std::string> THeavyMeson3pt<FImpl, GImpl>::getInput(void)
{
    std::vector<std::string> in = {par().gauge};
    auto add = [&in](const std::string &s)
    {
        if (!s.empty() and !heavyMeson3ptContains(in, s))
        {
            in.push_back(s);
        }
    };

    add(par().qSrc);
    add(par().qSnk);
    add(par().spectator);
    add(par().sinkNoise);
    if (!par().momProjector.empty())
    {
        add(par().momProjector);
        add(par().momProjector + "_momList");
    }

    return in;
}

template <typename FImpl, typename GImpl>
std::vector<std::string> THeavyMeson3pt<FImpl, GImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};

    return out;
}

template <typename FImpl, typename GImpl>
std::vector<std::string> THeavyMeson3pt<FImpl, GImpl>::getOutputFiles(void)
{
    std::vector<std::string> output;

    if (!par().output.empty())
    {
        output.push_back(resultFilename(par().output));
    }

    return output;
}

// parsers /////////////////////////////////////////////////////////////////////
// Fill triads_ (and the sink/vertex/source sets) from par().gamma, with
// the same "sink vertex source" format and "all" expansion as Gamma3pt.
template <typename FImpl, typename GImpl>
void THeavyMeson3pt<FImpl, GImpl>::parseHQEGammas(void)
{
    auto gammaName = [](const Gamma::Algebra g)
    {
        std::ostringstream s;

        s << g;

        return s.str();
    };
    auto expandGamma = [](const std::string &tok)
    {
        std::vector<Gamma::Algebra> v;

        if (tok == "all")
        {
            for (unsigned int i = 1; i < Gamma::nGamma; i += 2)
            {
                v.push_back(static_cast<Gamma::Algebra>(i));
            }
        }
        else
        {
            v = strToVec<Gamma::Algebra>(tok);
            if ((v.size() != 1) or (v[0] == Gamma::Algebra::undef))
            {
                HADRONS_ERROR(Argument, "invalid gamma matrix '" + tok + "'");
            }
        }

        return v;
    };
    auto expandVertex = [&](const std::string &tok)
    {
        std::vector<std::string> v;

        if (tok == "ops")
        {
            v = heavyMeson3ptOps;
        }
        else if (heavyMeson3ptContains(heavyMeson3ptOps, tok))
        {
            v.push_back(tok);
        }
        else
        {
            for (auto g: expandGamma(tok))
            {
                v.push_back(gammaName(g));
            }
        }

        return v;
    };

    triads_.clear();
    sinks_.clear();
    sources_.clear();
    vertices_.clear();
    for (auto &str: par().gamma)
    {
        std::vector<std::string> tok = strToVec<std::string>(str);

        if (tok.size() != 3)
        {
            HADRONS_ERROR(Argument, "gamma entry '" + str
                          + "' must be a 'sink vertex source' triad");
        }

        auto snks = expandGamma(tok[0]);
        auto vtxs = expandVertex(tok[1]);
        auto srcs = expandGamma(tok[2]);

        for (auto src: srcs)
        for (auto snk: snks)
        for (auto &vtx: vtxs)
        {
            HQETriad t;

            t.sink   = snk;
            t.vertex = vtx;
            t.source = src;
            if (!heavyMeson3ptContains(triads_, t))
            {
                triads_.push_back(t);
            }
            if (!heavyMeson3ptContains(sinks_, snk))
            {
                sinks_.push_back(snk);
            }
            if (!heavyMeson3ptContains(sources_, src))
            {
                sources_.push_back(src);
            }
            if (!heavyMeson3ptContains(vertices_, vtx))
            {
                vertices_.push_back(vtx);
            }
        }
    }
}

template <typename FImpl, typename GImpl>
bool THeavyMeson3pt<FImpl, GImpl>::usedPair(const std::string &v,
                                             const Gamma::Algebra src) const
{
    for (auto &t: triads_)
    {
        if ((t.vertex == v) and (t.source == src))
        {
            return true;
        }
    }

    return false;
}

// add the momentum projection of c to acc[momentum][time]
template <typename FImpl, typename GImpl>
void THeavyMeson3pt<FImpl, GImpl>::project(MomCorr &acc, ComplexField &c)
{
    const unsigned int                                    nt = env().getDim(Tp);
    std::vector<typename ComplexField::scalar_object>     modes;

    startTimer("FT");
    if (zeroMomentumOnly_)
    {
        sliceSum(c, modes, Tp);
    }
    else
    {
        auto &mp = envGet(MPType, par().momProjector);

        mp.Project(c, modes);
    }
    stopTimer("FT");
    if (acc.empty())
    {
        acc.assign(nMom_, std::vector<Complex>(nt, Complex(0., 0.)));
    }
    for (unsigned int m = 0; m < nMom_; ++m)
    for (unsigned int t = 0; t < nt; ++t)
    {
        acc[m][t] += TensorRemove(modes[nt*m + t]);
    }
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl, typename GImpl>
void THeavyMeson3pt<FImpl, GImpl>::setup(void)
{
    GridBase *grid = envGetGrid(PropagatorField);

    if (par().gauge.empty())
    {
        HADRONS_ERROR(Argument, "gauge must be set");
    }
    if ((par().bc != 1) and (par().bc != -1))
    {
        HADRONS_ERROR(Argument, "bc must be +1 or -1");
    }
    if (par().flowTime < 0.)
    {
        HADRONS_ERROR(Argument, "flowTime must be non-negative");
    }
    if (par().qSrc.empty() or par().qSnk.empty() or par().spectator.empty()
        or par().sinkNoise.empty() or par().gamma.empty())
    {
        HADRONS_ERROR(Argument, "qSrc, qSnk, spectator, sinkNoise and gamma must be set");
    }
    if (par().tSnk >= env().getDim(Tp))
    {
        HADRONS_ERROR(Argument, "tSnk out of range");
    }
    zeroMomentumOnly_ = par().momProjector.empty();
    needF_            = false;
    parseHQEGammas();
    for (auto &v: vertices_)
    {
        needF_ = needF_ or heavyMeson3ptContains(heavyMeson3ptFOps, v);
    }
    envTmpLat(PropagatorField, "tmp");
    envTmpLat(PropagatorField, "tmp2");
    envTmpLat(PropagatorField, "Y");
    envTmpLat(PropagatorField, "dq");
    envTmpLat(PropagatorField, "d4q3");
    envTmp(std::vector<PropagatorField>, "prodSrc", 1, sources_.size(), grid);
    envTmp(std::vector<SpinMatrixField>, "Z", 1, sources_.size(), grid);
    envTmpLat(GaugeLinkField, "linkT");
    envTmpLat(GaugeLinkField, "link");
    if (needF_)
    {
        envTmp(std::vector<GaugeLinkField>, "Es", 1, 3, grid);
        envTmp(std::vector<GaugeLinkField>, "Bs", 1, 3, grid);
        envTmpLat(GaugeLinkField, "divE");
    }
    envTmpLat(ComplexField, "contraction");
    envCreate(HadronsSerializable, getName(), 1, 0);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl, typename GImpl>
void THeavyMeson3pt<FImpl, GImpl>::execute(void)
{
    auto   &U = envGet(GaugeField, par().gauge);
    Result result;

    LOG(Message) << "Computing heavy-meson 3pt contractions at flow time "
                 << par().flowTime << " with gauge field '" << par().gauge
                 << "'" << std::endl;
    result.flowTime = par().flowTime;
    result.tSnk     = par().tSnk;
    if (zeroMomentumOnly_)
    {
        result.momentum = {{0, 0, 0}};
    }
    else
    {
        result.momentum = envGet(std::vector<std::vector<int>>,
                                 par().momProjector + "_momList");
    }
    nMom_ = result.momentum.size();
    computeHQE(result, U);
    startTimer("I/O");
    saveResult(par().output, "heavyMeson3pt", result);
    auto &out = envGet(HadronsSerializable, getName());
    out = result;
    stopTimer("I/O");
}

// two-quark (HQE) part ////////////////////////////////////////////////////////
template <typename FImpl, typename GImpl>
void THeavyMeson3pt<FImpl, GImpl>::computeHQE(Result &res, const GaugeField &U)
{
    const HeavyMeson3ptPar    &hp = par();
    auto                      &q2 = envGet(PropagatorField, hp.qSrc);
    auto                      &q3 = envGet(PropagatorField, hp.qSnk);
    auto                      &ql = envGet(PropagatorField, hp.spectator);
    Gamma                     g5(Gamma::Algebra::Gamma5);
    const unsigned int        nSrc = sources_.size();

    envGetTmp(PropagatorField, tmp);
    envGetTmp(PropagatorField, tmp2);
    envGetTmp(PropagatorField, Y);
    envGetTmp(PropagatorField, dq);
    envGetTmp(PropagatorField, d4q3);
    envGetTmp(std::vector<PropagatorField>, prodSrc);
    envGetTmp(std::vector<SpinMatrixField>, Z);
    envGetTmp(GaugeLinkField, linkT);
    envGetTmp(GaugeLinkField, link);
    envGetTmp(ComplexField, contraction);

    LOG(Message) << "qSrc '" << hp.qSrc << "', qSnk '" << hp.qSnk
                 << "', spectator '" << hp.spectator << "', sink noise '"
                 << hp.sinkNoise << "', tSnk " << hp.tSnk << ", "
                 << triads_.size() << " correlator(s)" << std::endl;

    res.gamma = triads_;
    res.corr.assign(triads_.size(), MomCorr());

    // project tr[gSnk Z] for every sink paired with (vertex, source)
    auto emit = [&](const std::string &v, const unsigned int s, const SpinMatrixField &Zs)
    {
        for (auto snk: sinks_)
        {
            HQETriad cur;

            cur.sink   = snk;
            cur.vertex = v;
            cur.source = sources_[s];

            auto loc = std::find(triads_.begin(), triads_.end(), cur);

            if (loc == triads_.end())
            {
                continue;
            }
            startTimer("Trace");
            contraction = traceSpin(Gamma(snk)*Zs);
            stopTimer("Trace");
            project(res.corr[loc - triads_.begin()], contraction);
        }
    };

    // spectator closed at the sink wall
    startTimer("Sink");
    SlicedPropagator snkSlice;
    auto             &eta = envGet(PropagatorField, hp.sinkNoise);

    tmp = adj(eta)*ql;
    sliceSum(tmp, snkSlice, Tp);
    SitePropagator q1Snk = snkSlice[hp.tSnk];
    stopTimer("Sink");

    // left blocks q1Snk (adj(gSrc) g5) adj(q2), one per source gamma
    startTimer("Trace");
    tmp = adj(q2);
    for (unsigned int s = 0; s < nSrc; ++s)
    {
        Gamma gSrc(sources_[s]);

        prodSrc[s] = q1Snk*(adj(gSrc)*g5)*tmp;
    }
    stopTimer("Trace");

    // shared building blocks
    startTimer("Insertion");
    boundaryLink<GImpl>(linkT, U, Tp, par().bc);

    std::vector<GaugeLinkField> *E = nullptr, *B = nullptr;
    if (needF_)
    {
        envGetTmp(std::vector<GaugeLinkField>, Es);
        envGetTmp(std::vector<GaugeLinkField>, Bs);
        HQE::chromoE<GImpl>(Es, U);
        HQE::chromoB<GImpl>(Bs, U);
        E = &Es;
        B = &Bs;
    }

    bool haveD4 = false;
    auto ensureD4 = [&](void)
    {
        if (!haveD4)
        {
            covDeriv<GImpl>(d4q3, q3, linkT, Tp, DerivType::symmetric);
            haveD4 = true;
        }
    };
    stopTimer("Insertion");

    const std::vector<int> spatial = {0, 1, 2};

    // Y = O q3 for a vertex applied directly to q3
    auto applyDirect = [&](const std::string &v)
    {
        startTimer("Insertion");
        if (v == "D4")
        {
            ensureD4();
            Y = d4q3;
        }
        else if (v == "D4D4")
        {
            laplacian<GImpl>(Y, q3, U, {Tp}, par().bc);
        }
        else if (v == "D4D4D4")
        {
            laplacian<GImpl>(tmp, q3, U, {Tp}, par().bc);
            covDeriv<GImpl>(Y, tmp, linkT, Tp, DerivType::symmetric);
        }
        else if (v == "sigmaB")
        {
            HQE::applySigmaB<GImpl>(Y, q3, *B);
        }
        else if (v == "sigmaE")
        {
            HQE::applySigmaE<GImpl>(Y, q3, *E);
        }
        else if (v == "sigmaBD4")
        {
            ensureD4();
            HQE::applySigmaB<GImpl>(Y, d4q3, *B);
        }
        else if (v == "sigmaED4")
        {
            ensureD4();
            HQE::applySigmaE<GImpl>(Y, d4q3, *E);
        }
        else if (v == "divE")
        {
            envGetTmp(GaugeLinkField, divE);
            HQE::divE<GImpl>(divE, *E, U);
            Y = divE*q3;
        }
        else if (v == "D2")
        {
            laplacian<GImpl>(Y, q3, U, spatial, par().bc);
        }
        else if (v == "D2D4")
        {
            ensureD4();
            laplacian<GImpl>(Y, d4q3, U, spatial, par().bc);
        }
        else if (v == "DiE")
        {
            Y = Zero();
            for (int i: spatial)
            {
                link = PeekIndex<LorentzIndex>(U, i);
                tmp  = (*E)[i]*q3;
                covDeriv<GImpl>(tmp2, tmp, link, i, DerivType::symmetric);
                Y   += tmp2;
            }
        }
        else if (v == "DiSigmaE")
        {
            Y = Zero();
            for (int i: spatial)
            {
                link = PeekIndex<LorentzIndex>(U, i);
                tmp  = Zero();
                for (int j: spatial)
                {
                    if (j != i)
                    {
                        HQE::applySigma(tmp2, q3, i, j);
                        tmp += (*E)[j]*tmp2;
                    }
                }
                covDeriv<GImpl>(tmp2, tmp, link, i, DerivType::symmetric);
                Y += tmp2;
            }
        }
        else if (v == "gD4")
        {
            ensureD4();
            Y = Gamma::gmu[Tp]*d4q3;
        }
        else if (v == "gD4D4")
        {
            laplacian<GImpl>(tmp, q3, U, {Tp}, par().bc);
            Y = Gamma::gmu[Tp]*tmp;
        }
        else if (v == "gD2")
        {
            laplacian<GImpl>(tmp, q3, U, spatial, par().bc);
            Y = Gamma::gmu[Tp]*tmp;
        }
        else if ((v == "gDi") or (v == "gDiD4"))
        {
            // D_i acts on q3 (gDi) or on D_4 q3 (gDiD4)
            if (v == "gDiD4")
            {
                ensureD4();
            }
            const PropagatorField &f = (v == "gDi") ? q3 : d4q3;

            Y = Zero();
            for (int i: spatial)
            {
                link = PeekIndex<LorentzIndex>(U, i);
                covDeriv<GImpl>(tmp, f, link, i, DerivType::symmetric);
                Y += Gamma::gmu[i]*tmp;
            }
        }
        else if (v == "gSigmaB")
        {
            HQE::applySigmaB<GImpl>(tmp, q3, *B);
            Y = Gamma::gmu[Tp]*tmp;
        }
        else if (v == "gSigmaE")
        {
            HQE::applySigmaE<GImpl>(tmp, q3, *E);
            Y = Gamma::gmu[Tp]*tmp;
        }
        else
        {
            // gamma vertex: tr[... g5 gV q3], as in Gamma3pt
            Gamma gV(strToVec<Gamma::Algebra>(v)[0]);

            Y = gV*q3;
        }
        stopTimer("Insertion");
    };

    // vertices applied directly to q3 (or, for D4L, to the left block)
    for (auto &v: vertices_)
    {
        if (v == "D4L")
        {
            startTimer("Insertion");
            covDeriv<GImpl>(dq, q2, linkT, Tp, DerivType::symmetric);
            tmp = adj(dq);
            stopTimer("Insertion");
            for (unsigned int s = 0; s < nSrc; ++s)
            {
                if (!usedPair(v, sources_[s]))
                {
                    continue;
                }

                Gamma gSrc(sources_[s]);

                startTimer("Trace");
                tmp2 = q1Snk*(adj(gSrc)*g5)*tmp;
                Z[s] = traceColour(tmp2*(g5*q3));
                stopTimer("Trace");
                emit(v, s, Z[s]);
            }
            continue;
        }
        applyDirect(v);
        for (unsigned int s = 0; s < nSrc; ++s)
        {
            if (!usedPair(v, sources_[s]))
            {
                continue;
            }
            startTimer("Trace");
            Z[s] = traceColour(prodSrc[s]*(g5*Y));
            stopTimer("Trace");
            emit(v, s, Z[s]);
        }
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_HeavyMeson3pt_hpp_
