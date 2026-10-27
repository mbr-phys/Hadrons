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
#include <Hadrons/Modules/MContraction/HQEOperators.hpp>

BEGIN_HADRONS_NAMESPACE

/*
 * Heavy-meson three-point functions with HQE two-quark and four-quark
 * operator insertions, for one (source wall, sink wall) pair at one flow time.
 *
 * All operator fields (propagators at the operator, gauge links, field
 * strength) are taken at the same flow time; the spectator line of the
 * two-quark part is unflowed. Derivatives of each propagator leg are built
 * once per direction and shared by all operators; conventions are those of
 * MContraction/HQEOperators.hpp.
 *
 * ---------------------------------------------------------------------------
 * Two-quark (HQE) part
 * ---------------------------------------------------------------------------
 *
 *   C[X](t) = sum_{x, x_t = t} tr[ W(x) X q3(x) ]
 *   W(x)    = gA K gB g5 q2(x)^dag g5
 *   K       = sum_{z, z_t = tSnk} eta(z)^dag S_l(z)       (sinkNoise = eta)
 *           = sum_{z, z_t = tSnk} S_l(z)                  (sinkNoise empty)
 *
 *   q2 = qSrc (heavy from the source wall), q3 = qSnk (heavy from the sink
 *   wall), S_l = spectator (light from the source wall, unflowed). With
 *   gA = gB = Gamma5 and X = Gamma this is MContraction::Gamma3pt. The noise
 *   sink uses eta^dag, as required for complex Z2 wall noise.
 *
 *   gammas: plain insertions X = Gamma ("all" = the 16 Gamma matrices)
 *   ops (space separated, or "all"):
 *     D2       sum_{i<3} (D^-_i D^+_i)            nearest-neighbour Laplacian
 *     D4       D_4                                 symmetric
 *     D4L      left-acting D_4 on qbar             (D4 + D4L ~ d/dt C[1], up to O(a^2))
 *     D4D4     D^-_4 D^+_4
 *     D4D4D4   D_4 D^-_4 D^+_4
 *     D2D4     sum_i D^-_i D^+_i D_4
 *     sigmaB   sum_{i<j<3} sigma_ij F_ij
 *     sigmaE   sum_{j<3} sigma_j3 F_j3
 *     sigmaBD4 sigmaB D_4
 *     sigmaED4 sigmaE D_4
 *     DiE      sum_i D_i F_i3                      (derivative acts on F q)
 *     divE     sum_i [D_i, F_i3]                   (colour-matrix insertion)
 *     DiSigmaE sum_{i,j} D_i sigma_ij F_j3         (derivative acts on F q)
 *   Spatial derivatives are applied by exact summation by parts on the time
 *   slice (periodic spatial boundaries), e.g.
 *     sum_x tr[W D^-_i D^+_i q3] = -sum_x tr[(W <-D^+_i)(D^+_i q3)],
 *   so each leg is differentiated at most once per direction.
 *
 * ---------------------------------------------------------------------------
 * Four-quark part (generalises MContraction::WeakNonEye3pt)
 * ---------------------------------------------------------------------------
 *
 *   L = qLeft adj(gIn) g5 adj(qBarLeft) g5
 *   R = qBarRight gOut g5 adj(qRight) g5
 *   one trace : tr[L G1 R G2]
 *   two traces: tr[L G1] tr[R G2]
 *
 *   With G1 = G2 = G this is exactly WeakNonEye3pt. gammaPairs: "diag" (the 16
 *   G x G, default), "all" (256 pairs) or an explicit list "(G1 G2)(G1 G2)...".
 *
 *   derivatives: list of specifications {legs, dirs, couple}
 *     legs:   subset of "qBarLeft qLeft qBarRight qRight" (empty = no derivative)
 *     dirs:   summed directions rho (default "0 1 2 3")
 *     couple: ""   -> sum_rho contraction with D_rho on every listed leg
 *             "G1" -> as above with G1 -> gamma_rho G1
 *             "G2" -> as above with G2 -> gamma_rho G2
 *   Derivatives are symmetric; on the adjoint legs (qBarLeft, qRight) they act
 *   on the antiquark field to the left. The summation over rho is Euclidean.
 *   If no specification is given, only the underived contraction is computed.
 *
 * Parameters:
 * - gauge:     flowed gauge field (links and field strength)
 * - bc:        fermion temporal boundary phase for time derivatives (+1, -1)
 * - flowTime:  label stored with the result
 * - hqe:       two-quark part (skipped if qSnk is empty)
 * - fourQuark: list of four-quark contractions
 * - output:    result file stem
 */

/******************************************************************************
 *                             HeavyMeson3pt                                  *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MContraction)

class HeavyMeson3ptHQEPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(HeavyMeson3ptHQEPar,
                                    std::string,    qSrc,
                                    std::string,    qSnk,
                                    std::string,    spectator,
                                    std::string,    sinkNoise,
                                    unsigned int,   tSnk,
                                    Gamma::Algebra, gammaSnk,
                                    Gamma::Algebra, gammaSrc,
                                    std::string,    gammas,
                                    std::string,    ops);
    HeavyMeson3ptHQEPar(void)
    : tSnk(0), gammaSnk(Gamma::Algebra::Gamma5)
    , gammaSrc(Gamma::Algebra::Gamma5) {}
};

class HeavyMeson3ptDerivPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(HeavyMeson3ptDerivPar,
                                    std::string, legs,
                                    std::string, dirs,
                                    std::string, couple);
};

class HeavyMeson3ptFourQuarkPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(HeavyMeson3ptFourQuarkPar,
                                    std::string,    name,
                                    std::string,    qBarLeft,
                                    std::string,    qLeft,
                                    std::string,    qBarRight,
                                    std::string,    qRight,
                                    Gamma::Algebra, gammaIn,
                                    Gamma::Algebra, gammaOut,
                                    std::string,    gammaPairs,
                                    std::vector<HeavyMeson3ptDerivPar>, derivatives);
    HeavyMeson3ptFourQuarkPar(void)
    : gammaIn(Gamma::Algebra::Gamma5), gammaOut(Gamma::Algebra::Gamma5)
    , gammaPairs("diag") {}
};

class HeavyMeson3ptPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(HeavyMeson3ptPar,
                                    std::string,         gauge,
                                    int,                 bc,
                                    double,              flowTime,
                                    HeavyMeson3ptHQEPar, hqe,
                                    std::vector<HeavyMeson3ptFourQuarkPar>, fourQuark,
                                    std::string,         output);
    HeavyMeson3ptPar(void): bc(-1), flowTime(0.) {}
};

template <typename FImpl, typename GImpl>
class THeavyMeson3pt: public Module<HeavyMeson3ptPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    typedef typename GImpl::GaugeLinkField GaugeLinkField;
    typedef std::pair<Gamma::Algebra, Gamma::Algebra> GammaPair;
    class HQEResult: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(HQEResult,
                                        std::string,          op,
                                        Gamma::Algebra,       gamma,
                                        std::vector<Complex>, corr);
    };
    class FourQuarkResult: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(FourQuarkResult,
                                        std::string,          name,
                                        std::string,          deriv,
                                        Gamma::Algebra,       gamma1,
                                        Gamma::Algebra,       gamma2,
                                        unsigned int,         trace,
                                        std::vector<Complex>, corr);
    };
    class Result: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(Result,
                                        double,                       flowTime,
                                        unsigned int,                 tSnk,
                                        std::vector<HQEResult>,       hqe,
                                        std::vector<FourQuarkResult>, fourQuark);
    };
private:
    // parsed derivative specification
    struct DerivSpec
    {
        std::string      label;
        bool             leg[4];
        std::vector<int> dirs;
        int              couple; // 0: none, 1: G1, 2: G2
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
    void computeHQE(std::vector<HQEResult> &res, const GaugeField &U);
    void computeFourQuark(std::vector<FourQuarkResult> &res, const GaugeField &U);
    std::vector<std::string> parseOps(const std::string &ops) const;
    std::vector<GammaPair>   parseGammaPairs(const std::string &str) const;
    std::vector<DerivSpec>   parseDerivs(const HeavyMeson3ptFourQuarkPar &fq) const;
    void addSlice(std::vector<Complex> &acc, const ComplexField &c) const;
};

MODULE_REGISTER_TMP(HeavyMeson3pt, ARG(THeavyMeson3pt<FIMPL, GIMPL>), MContraction);

/******************************************************************************
 *                     THeavyMeson3pt implementation                          *
 ******************************************************************************/
static const std::vector<std::string> heavyMeson3ptOps =
    {"D2", "D4", "D4L", "D4D4", "D4D4D4", "D2D4", "sigmaB", "sigmaE",
     "sigmaBD4", "sigmaED4", "DiE", "divE", "DiSigmaE"};

static const std::vector<std::string> heavyMeson3ptLegs =
    {"qBarLeft", "qLeft", "qBarRight", "qRight"};

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
        if (!s.empty() and (std::find(in.begin(), in.end(), s) == in.end()))
        {
            in.push_back(s);
        }
    };

    if (!par().hqe.qSnk.empty())
    {
        add(par().hqe.qSrc);
        add(par().hqe.qSnk);
        add(par().hqe.spectator);
        add(par().hqe.sinkNoise);
    }
    for (auto &fq: par().fourQuark)
    {
        add(fq.qBarLeft);
        add(fq.qLeft);
        add(fq.qBarRight);
        add(fq.qRight);
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
template <typename FImpl, typename GImpl>
std::vector<std::string> THeavyMeson3pt<FImpl, GImpl>::parseOps(const std::string &ops) const
{
    if (ops == "all")
    {
        return heavyMeson3ptOps;
    }

    std::vector<std::string> list = strToVec<std::string>(ops);

    for (auto &o: list)
    {
        if (std::find(heavyMeson3ptOps.begin(), heavyMeson3ptOps.end(), o)
            == heavyMeson3ptOps.end())
        {
            HADRONS_ERROR(Argument, "unknown HQE operator '" + o + "'");
        }
    }

    return list;
}

template <typename FImpl, typename GImpl>
std::vector<typename THeavyMeson3pt<FImpl, GImpl>::GammaPair>
THeavyMeson3pt<FImpl, GImpl>::parseGammaPairs(const std::string &str) const
{
    std::vector<GammaPair> pairs;

    if (str.empty() or (str == "diag"))
    {
        for (auto &g: Gamma::gall)
        {
            pairs.push_back({g.g, g.g});
        }
    }
    else if (str == "all")
    {
        for (auto &g1: Gamma::gall)
        for (auto &g2: Gamma::gall)
        {
            pairs.push_back({g1.g, g2.g});
        }
    }
    else
    {
        pairs = strToVec<GammaPair>(str);
        for (auto &p: pairs)
        {
            if ((p.first == Gamma::Algebra::undef) or (p.second == Gamma::Algebra::undef))
            {
                HADRONS_ERROR(Argument, "wrong gamma pair list '" + str + "'");
            }
        }
    }

    return pairs;
}

template <typename FImpl, typename GImpl>
std::vector<typename THeavyMeson3pt<FImpl, GImpl>::DerivSpec>
THeavyMeson3pt<FImpl, GImpl>::parseDerivs(const HeavyMeson3ptFourQuarkPar &fq) const
{
    std::vector<DerivSpec> specs;

    for (auto &d: fq.derivatives)
    {
        DerivSpec                s;
        std::vector<std::string> legs = strToVec<std::string>(d.legs);

        for (unsigned int k = 0; k < 4; ++k)
        {
            s.leg[k] = false;
        }
        for (auto &l: legs)
        {
            auto it = std::find(heavyMeson3ptLegs.begin(), heavyMeson3ptLegs.end(), l);

            if (it == heavyMeson3ptLegs.end())
            {
                HADRONS_ERROR(Argument, "unknown leg '" + l + "' (expected qBarLeft, "
                              "qLeft, qBarRight or qRight)");
            }
            s.leg[it - heavyMeson3ptLegs.begin()] = true;
        }
        s.dirs = d.dirs.empty() ? std::vector<int>({0, 1, 2, 3}) : strToVec<int>(d.dirs);
        for (int rho: s.dirs)
        {
            if ((rho < 0) or (rho >= Nd))
            {
                HADRONS_ERROR(Argument, "direction out of range: " + d.dirs);
            }
        }
        if (d.couple.empty())
        {
            s.couple = 0;
        }
        else if (d.couple == "G1")
        {
            s.couple = 1;
        }
        else if (d.couple == "G2")
        {
            s.couple = 2;
        }
        else
        {
            HADRONS_ERROR(Argument, "couple must be empty, 'G1' or 'G2'");
        }
        if (legs.empty())
        {
            s.label = "none";
            s.dirs.clear();
            s.couple = 0;
        }
        else
        {
            std::string dl;

            for (int rho: s.dirs)
            {
                dl += std::to_string(rho);
            }
            s.label = d.legs + "@" + dl;
            if (s.couple > 0)
            {
                s.label += ":" + d.couple;
            }
        }
        specs.push_back(s);
    }
    if (specs.empty())
    {
        DerivSpec s;

        s.label  = "none";
        s.couple = 0;
        for (unsigned int k = 0; k < 4; ++k)
        {
            s.leg[k] = false;
        }
        specs.push_back(s);
    }

    return specs;
}

template <typename FImpl, typename GImpl>
void THeavyMeson3pt<FImpl, GImpl>::addSlice(std::vector<Complex> &acc,
                                             const ComplexField &c) const
{
    SlicedComplex buf;

    sliceSum(c, buf, Tp);
    if (acc.empty())
    {
        acc.assign(buf.size(), Complex(0., 0.));
    }
    for (unsigned int t = 0; t < buf.size(); ++t)
    {
        acc[t] += TensorRemove(buf[t]);
    }
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl, typename GImpl>
void THeavyMeson3pt<FImpl, GImpl>::setup(void)
{
    if ((par().bc != 1) and (par().bc != -1))
    {
        HADRONS_ERROR(Argument, "bc must be +1 or -1");
    }
    if (!par().hqe.qSnk.empty())
    {
        if (par().hqe.qSrc.empty() or par().hqe.spectator.empty())
        {
            HADRONS_ERROR(Argument, "hqe needs qSrc, qSnk and spectator");
        }
        if (par().hqe.tSnk >= env().getDim(Tp))
        {
            HADRONS_ERROR(Argument, "hqe.tSnk out of range");
        }
        parseOps(par().hqe.ops);
    }
    for (auto &fq: par().fourQuark)
    {
        parseGammaPairs(fq.gammaPairs);
        parseDerivs(fq);
    }
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
    result.tSnk     = par().hqe.tSnk;
    if (!par().hqe.qSnk.empty())
    {
        computeHQE(result.hqe, U);
    }
    if (!par().fourQuark.empty())
    {
        computeFourQuark(result.fourQuark, U);
    }
    saveResult(par().output, "heavyMeson3pt", result);
    auto &out = envGet(HadronsSerializable, getName());
    out = result;
}

// two-quark (HQE) part ////////////////////////////////////////////////////////
template <typename FImpl, typename GImpl>
void THeavyMeson3pt<FImpl, GImpl>::computeHQE(std::vector<HQEResult> &res,
                                               const GaugeField &U)
{
    const HeavyMeson3ptHQEPar &hp  = par().hqe;
    auto                      &q2  = envGet(PropagatorField, hp.qSrc);
    auto                      &q3  = envGet(PropagatorField, hp.qSnk);
    auto                      &ql  = envGet(PropagatorField, hp.spectator);
    GridBase                  *grid = q3.Grid();
    Gamma                     g5(Gamma::Algebra::Gamma5);
    Gamma                     gA(hp.gammaSnk), gB(hp.gammaSrc);
    std::vector<std::string>  ops = parseOps(hp.ops);
    auto has = [&ops](const std::string &o)
    {
        return std::find(ops.begin(), ops.end(), o) != ops.end();
    };
    auto add = [this, &res](const std::string &op, const Gamma::Algebra g,
                            const ComplexField &c)
    {
        HQEResult r;

        r.op    = op;
        r.gamma = g;
        addSlice(r.corr, c);
        res.push_back(r);
    };

    LOG(Message) << "HQE part: qSrc '" << hp.qSrc << "', qSnk '" << hp.qSnk
                 << "', spectator '" << hp.spectator << "', sink noise '"
                 << hp.sinkNoise << "', tSnk " << hp.tSnk << std::endl;

    // spectator closed at the sink wall
    SlicedPropagator snk;
    {
        PropagatorField tmp(grid);

        if (hp.sinkNoise.empty())
        {
            tmp = ql;
        }
        else
        {
            auto &eta = envGet(PropagatorField, hp.sinkNoise);

            tmp = adj(eta)*ql;
        }
        sliceSum(tmp, snk, Tp);
    }
    SitePropagator K = gA*snk[hp.tSnk]*gB;

    PropagatorField W(grid), Wd(grid), qd(grid), tmp(grid), tmp2(grid);
    ComplexField    c(grid);
    GaugeLinkField  linkT(grid), link(grid);

    // W(x) and its left derivatives W <-D = K g5 (D q2)^dag g5
    W = K*(g5*adj(q2)*g5);
    HQE::boundaryLink<GImpl>(linkT, U, Tp, par().bc);

    // plain gamma insertions
    if (!hp.gammas.empty())
    {
        std::vector<Gamma::Algebra> gl;

        if (hp.gammas == "all")
        {
            for (auto &g: Gamma::gall)
            {
                gl.push_back(g.g);
            }
        }
        else
        {
            gl = strToVec<Gamma::Algebra>(hp.gammas);
        }
        for (auto g: gl)
        {
            c = trace(W*(Gamma(g)*q3));
            add("gamma", g, c);
        }
    }

    // field strength
    std::vector<GaugeLinkField> F;
    if (has("sigmaB") or has("sigmaE") or has("sigmaBD4") or has("sigmaED4")
        or has("DiE") or has("divE") or has("DiSigmaE"))
    {
        HQE::allFieldStrengths<GImpl>(F, U);
    }

    // time derivatives
    PropagatorField d4q3(grid);
    if (has("D4") or has("D2D4") or has("D4D4D4") or has("sigmaBD4")
        or has("sigmaED4"))
    {
        HQE::covDeriv<GImpl>(d4q3, q3, linkT, Tp, HQE::DerivType::symmetric);
    }
    if (has("D4"))
    {
        c = trace(W*d4q3);
        add("D4", Gamma::Algebra::Identity, c);
    }
    if (has("D4L"))
    {
        HQE::covDeriv<GImpl>(qd, q2, linkT, Tp, HQE::DerivType::symmetric);
        Wd = K*(g5*adj(qd)*g5);
        c  = trace(Wd*q3);
        add("D4L", Gamma::Algebra::Identity, c);
    }
    if (has("D4D4") or has("D4D4D4"))
    {
        HQE::laplacian<GImpl>(tmp2, q3, U, {Tp}, par().bc);
        if (has("D4D4"))
        {
            c = trace(W*tmp2);
            add("D4D4", Gamma::Algebra::Identity, c);
        }
        if (has("D4D4D4"))
        {
            HQE::covDeriv<GImpl>(tmp, tmp2, linkT, Tp, HQE::DerivType::symmetric);
            c = trace(W*tmp);
            add("D4D4D4", Gamma::Algebra::Identity, c);
        }
    }

    // local gluonic insertions
    if (has("sigmaB"))
    {
        HQE::applySigmaB<GImpl>(tmp, q3, F);
        c = trace(W*tmp);
        add("sigmaB", Gamma::Algebra::Identity, c);
    }
    if (has("sigmaE"))
    {
        HQE::applySigmaE<GImpl>(tmp, q3, F);
        c = trace(W*tmp);
        add("sigmaE", Gamma::Algebra::Identity, c);
    }
    if (has("sigmaBD4"))
    {
        HQE::applySigmaB<GImpl>(tmp, d4q3, F);
        c = trace(W*tmp);
        add("sigmaBD4", Gamma::Algebra::Identity, c);
    }
    if (has("sigmaED4"))
    {
        HQE::applySigmaE<GImpl>(tmp, d4q3, F);
        c = trace(W*tmp);
        add("sigmaED4", Gamma::Algebra::Identity, c);
    }
    if (has("divE"))
    {
        GaugeLinkField divE(grid);

        HQE::fieldStrengthDivergence<GImpl>(divE, F, U, Tp, {0, 1, 2});
        tmp = divE*q3;
        c   = trace(W*tmp);
        add("divE", Gamma::Algebra::Identity, c);
    }

    // spatial derivatives, by parts, one direction at a time
    const bool fwd = has("D2") or has("D2D4");
    const bool sym = has("DiE") or has("DiSigmaE");
    if (fwd or sym)
    {
        ComplexField cD2(grid), cD2D4(grid), cDiE(grid), cDiSE(grid);

        cD2 = Zero(); cD2D4 = Zero(); cDiE = Zero(); cDiSE = Zero();
        for (int i = 0; i < Tp; ++i)
        {
            link = PeekIndex<LorentzIndex>(U, i);
            if (fwd)
            {
                HQE::covDeriv<GImpl>(qd, q2, link, i, HQE::DerivType::forward);
                Wd = K*(g5*adj(qd)*g5);
                if (has("D2"))
                {
                    HQE::covDeriv<GImpl>(tmp, q3, link, i, HQE::DerivType::forward);
                    cD2 -= trace(Wd*tmp);
                }
                if (has("D2D4"))
                {
                    HQE::covDeriv<GImpl>(tmp, d4q3, link, i, HQE::DerivType::forward);
                    cD2D4 -= trace(Wd*tmp);
                }
            }
            if (sym)
            {
                HQE::covDeriv<GImpl>(qd, q2, link, i, HQE::DerivType::symmetric);
                Wd = K*(g5*adj(qd)*g5);
                if (has("DiE"))
                {
                    tmp   = F[HQE::fsIndex(i, Tp)]*q3;
                    cDiE -= trace(Wd*tmp);
                }
                if (has("DiSigmaE"))
                {
                    tmp2 = Zero();
                    for (int j = 0; j < Tp; ++j)
                    {
                        if (j != i)
                        {
                            HQE::applySigma(tmp, q3, i, j);
                            tmp2 += F[HQE::fsIndex(j, Tp)]*tmp;
                        }
                    }
                    cDiSE -= trace(Wd*tmp2);
                }
            }
        }
        if (has("D2"))
        {
            add("D2", Gamma::Algebra::Identity, cD2);
        }
        if (has("D2D4"))
        {
            add("D2D4", Gamma::Algebra::Identity, cD2D4);
        }
        if (has("DiE"))
        {
            add("DiE", Gamma::Algebra::Identity, cDiE);
        }
        if (has("DiSigmaE"))
        {
            add("DiSigmaE", Gamma::Algebra::Identity, cDiSE);
        }
    }
}

// four-quark part /////////////////////////////////////////////////////////////
template <typename FImpl, typename GImpl>
void THeavyMeson3pt<FImpl, GImpl>::computeFourQuark(std::vector<FourQuarkResult> &res,
                                                     const GaugeField &U)
{
    GridBase *grid = U.Grid();
    Gamma    g5(Gamma::Algebra::Gamma5);
    unsigned int nEntry = par().fourQuark.size();

    // parsed input and accumulators acc[entry][spec][pair][trace]
    std::vector<std::vector<GammaPair>>  pairs(nEntry);
    std::vector<std::vector<DerivSpec>>  specs(nEntry);
    std::vector<std::array<const PropagatorField *, 4>> legs(nEntry);
    std::vector<std::vector<std::vector<std::array<std::vector<Complex>, 2>>>> acc(nEntry);
    std::vector<std::string> derivNames;

    for (unsigned int e = 0; e < nEntry; ++e)
    {
        auto &fq = par().fourQuark[e];
        const std::string names[4] = {fq.qBarLeft, fq.qLeft, fq.qBarRight, fq.qRight};

        pairs[e] = parseGammaPairs(fq.gammaPairs);
        specs[e] = parseDerivs(fq);
        for (unsigned int k = 0; k < 4; ++k)
        {
            legs[e][k] = &envGet(PropagatorField, names[k]);
        }
        acc[e].resize(specs[e].size());
        for (unsigned int s = 0; s < specs[e].size(); ++s)
        {
            acc[e][s].resize(pairs[e].size());
            for (unsigned int k = 0; k < 4; ++k)
            {
                if (specs[e][s].leg[k] and (std::find(derivNames.begin(),
                    derivNames.end(), names[k]) == derivNames.end()))
                {
                    derivNames.push_back(names[k]);
                }
            }
        }
        LOG(Message) << "Four-quark '" << fq.name << "': " << pairs[e].size()
                     << " gamma pair(s), " << specs[e].size()
                     << " derivative specification(s)" << std::endl;
    }

    PropagatorField L(grid), R(grid);
    ComplexField    c(grid);

    // contract one (entry, spec) with the given legs, rho < 0: no coupling
    auto contract = [&](const unsigned int e, const unsigned int s,
                        const std::array<const PropagatorField *, 4> &q,
                        const int rho)
    {
        auto  &fq = par().fourQuark[e];
        Gamma gIn(fq.gammaIn), gOut(fq.gammaOut);
        const DerivSpec &sp = specs[e][s];

        L = (*q[1])*adj(gIn)*g5*adj(*q[0])*g5;
        R = (*q[2])*gOut*g5*adj(*q[3])*g5;
        for (unsigned int p = 0; p < pairs[e].size(); ++p)
        {
            Gamma G1(pairs[e][p].first), G2(pairs[e][p].second);

            if ((sp.couple == 1) and (rho >= 0))
            {
                G1 = Gamma::gmu[rho]*G1;
            }
            else if ((sp.couple == 2) and (rho >= 0))
            {
                G2 = Gamma::gmu[rho]*G2;
            }
            c = trace(L*G1*R*G2);
            addSlice(acc[e][s][p][0], c);
            c = trace(L*G1)*trace(R*G2);
            addSlice(acc[e][s][p][1], c);
        }
    };

    // underived contractions
    for (unsigned int e = 0; e < nEntry; ++e)
    for (unsigned int s = 0; s < specs[e].size(); ++s)
    {
        if (specs[e][s].dirs.empty())
        {
            contract(e, s, legs[e], -1);
        }
    }

    // derived contractions, one direction at a time, each leg derived once
    if (!derivNames.empty())
    {
        std::vector<PropagatorField> dq(derivNames.size(), PropagatorField(grid));
        GaugeLinkField               link(grid);

        for (int rho = 0; rho < Nd; ++rho)
        {
            bool needed = false;

            for (unsigned int e = 0; e < nEntry; ++e)
            for (auto &sp: specs[e])
            {
                needed = needed or (std::find(sp.dirs.begin(), sp.dirs.end(), rho)
                                    != sp.dirs.end());
            }
            if (!needed)
            {
                continue;
            }
            HQE::boundaryLink<GImpl>(link, U, rho, par().bc);
            for (unsigned int n = 0; n < derivNames.size(); ++n)
            {
                auto &q = envGet(PropagatorField, derivNames[n]);

                HQE::covDeriv<GImpl>(dq[n], q, link, rho, HQE::DerivType::symmetric);
            }
            for (unsigned int e = 0; e < nEntry; ++e)
            {
                auto &fq = par().fourQuark[e];
                const std::string names[4] = {fq.qBarLeft, fq.qLeft, fq.qBarRight, fq.qRight};

                for (unsigned int s = 0; s < specs[e].size(); ++s)
                {
                    auto &sp = specs[e][s];

                    if (std::find(sp.dirs.begin(), sp.dirs.end(), rho) == sp.dirs.end())
                    {
                        continue;
                    }

                    std::array<const PropagatorField *, 4> q = legs[e];

                    for (unsigned int k = 0; k < 4; ++k)
                    {
                        if (sp.leg[k])
                        {
                            auto n = std::find(derivNames.begin(), derivNames.end(),
                                               names[k]) - derivNames.begin();

                            q[k] = &dq[n];
                        }
                    }
                    contract(e, s, q, rho);
                }
            }
        }
    }

    // store
    for (unsigned int e = 0; e < nEntry; ++e)
    for (unsigned int s = 0; s < specs[e].size(); ++s)
    for (unsigned int p = 0; p < pairs[e].size(); ++p)
    for (unsigned int tr = 0; tr < 2; ++tr)
    {
        FourQuarkResult r;

        r.name   = par().fourQuark[e].name;
        r.deriv  = specs[e][s].label;
        r.gamma1 = pairs[e][p].first;
        r.gamma2 = pairs[e][p].second;
        r.trace  = tr + 1;
        r.corr   = acc[e][s][p][tr];
        res.push_back(r);
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_HeavyMeson3pt_hpp_
