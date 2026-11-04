/*
 * WeakNonEye3ptDeriv.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
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
#ifndef Hadrons_MContraction_WeakNonEye3ptDeriv_hpp_
#define Hadrons_MContraction_WeakNonEye3ptDeriv_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Serialization.hpp>
#include <Hadrons/Modules/MUtilities/DerivHelpers.hpp>
#include <Grid/algorithms/blas/MomentumProject.h>

BEGIN_HADRONS_NAMESPACE

/*
 * Weak Hamiltonian meson 3-pt diagrams, non-eye topologies, generalised:
 * arbitrary gamma pairs G1 x G2, covariant derivatives on the legs and
 * optional momentum projection. Same legs and diagrams as
 * MContraction::WeakNonEye3pt:
 *
 *   L = qLeft adj(gIn) g5 adj(qBarLeft) g5
 *   R = qBarRight gOut g5 adj(qRight) g5
 *   one trace : tr[L G1 R G2]
 *   two traces: tr[L G1] tr[R G2]
 *
 * With gammaPairs = "diag" and the derivative specification legs = "none"
 * this is exactly WeakNonEye3pt.
 *
 * gammaPairs: "diag" (the 16 G x G), "all" (256 pairs) or an explicit list
 *             "(G1 G2)(G1 G2)..."
 *
 * derivatives: non-empty list of specifications {legs, dirs, couple}
 *   legs:   "none" (underived contraction; dirs and couple must be empty) or
 *           a subset of "qBarLeft qLeft qBarRight qRight"
 *   dirs:   summed directions rho, e.g. "0 1 2 3"
 *   couple: "none" -> sum_rho contraction with D_rho on every listed leg
 *           "G1"   -> as above with G1 -> gamma_rho G1
 *           "G2"   -> as above with G2 -> gamma_rho G2
 *   Derivatives are symmetric (MUtilities/DerivHelpers.hpp); on the adjoint
 *   legs (qBarLeft, qRight) they act on the antiquark field to the left.
 *   The summation over rho is Euclidean; rho = 3 uses the temporal
 *   boundary phase bc. Each distinct leg is differentiated once per rho.
 *   Labels: "none" or "<legs>@<dirs>[:G1|:G2]", e.g. "qBarLeft qLeft@0123".
 *
 * Parameters:
 * - qLeft, qBarLeft, qRight, qBarRight: legs, as in WeakNonEye3pt
 * - gammaIn, gammaOut: source and sink gammas, as in WeakNonEye3pt
 * - gammaPairs:   operator gamma pairs (see above)
 * - derivatives:  derivative specifications (see above)
 * - gauge:        gauge field for the derivatives (flowed, at the flow time
 *                 of the legs)
 * - bc:           fermion temporal boundary phase (+1 or -1)
 * - momProjector: MomentumProject module (e.g. MUtilities::MPScalar); empty
 *                 for a zero-momentum sum only
 * - output:       result file stem (empty = no file)
 *
 * There are no default values: every parameter except momProjector and
 * output must be set.
 */

/******************************************************************************
 *                          WeakNonEye3ptDeriv                                *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MContraction)

class WeakDerivSpec: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(WeakDerivSpec,
                                    std::string, legs,
                                    std::string, dirs,
                                    std::string, couple);
};

class WeakNonEye3ptDerivPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(WeakNonEye3ptDerivPar,
                                    std::string,                qLeft,
                                    std::string,                qBarLeft,
                                    std::string,                qRight,
                                    std::string,                qBarRight,
                                    Gamma::Algebra,             gammaIn,
                                    Gamma::Algebra,             gammaOut,
                                    std::string,                gammaPairs,
                                    std::vector<WeakDerivSpec>, derivatives,
                                    std::string,                gauge,
                                    int,                        bc,
                                    std::string,                momProjector,
                                    std::string,                output);
};

template <typename FImpl, typename GImpl>
class TWeakNonEye3ptDeriv: public Module<WeakNonEye3ptDerivPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    typedef typename GImpl::GaugeLinkField GaugeLinkField;
    typedef std::pair<Gamma::Algebra, Gamma::Algebra> GammaPair;
    using MPType = MomentumProject<ComplexField, ComplexField>;
    typedef std::vector<std::vector<Complex>> MomCorr;
    class Contraction: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(Contraction,
                                        std::string,    deriv,
                                        Gamma::Algebra, gamma1,
                                        Gamma::Algebra, gamma2,
                                        unsigned int,   trace,
                                        MomCorr,        corr);
    };
    class Result: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(Result,
                                        Gamma::Algebra,                in,
                                        Gamma::Algebra,                out,
                                        std::vector<std::vector<int>>, momentum,
                                        std::vector<Contraction>,      contractions);
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
    TWeakNonEye3ptDeriv(const std::string name);
    // destructor
    virtual ~TWeakNonEye3ptDeriv(void) {};
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
    std::vector<std::string> legNames(void) const;
    std::vector<GammaPair>   parseGammaPairs(void) const;
    std::vector<DerivSpec>   parseDerivs(void) const;
    void project(MomCorr &acc, ComplexField &c);
private:
    std::vector<std::string> derivNames_;
    bool                     zeroMomentumOnly_;
    unsigned int             nMom_;
};

MODULE_REGISTER_TMP(WeakNonEye3ptDeriv, ARG(TWeakNonEye3ptDeriv<FIMPL, GIMPL>), MContraction);

/******************************************************************************
 *                   TWeakNonEye3ptDeriv implementation                       *
 ******************************************************************************/
// legs in the order of the derivative specification
static const std::vector<std::string> weakDerivLegs =
    {"qBarLeft", "qLeft", "qBarRight", "qRight"};

// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl, typename GImpl>
TWeakNonEye3ptDeriv<FImpl, GImpl>::TWeakNonEye3ptDeriv(const std::string name)
: Module<WeakNonEye3ptDerivPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl, typename GImpl>
std::vector<std::string> TWeakNonEye3ptDeriv<FImpl, GImpl>::getInput(void)
{
    std::vector<std::string> in = {par().qLeft, par().qBarLeft,
                                   par().qRight, par().qBarRight, par().gauge};

    if (!par().momProjector.empty())
    {
        in.push_back(par().momProjector);
        in.push_back(par().momProjector + "_momList");
    }

    return in;
}

template <typename FImpl, typename GImpl>
std::vector<std::string> TWeakNonEye3ptDeriv<FImpl, GImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};

    return out;
}

template <typename FImpl, typename GImpl>
std::vector<std::string> TWeakNonEye3ptDeriv<FImpl, GImpl>::getOutputFiles(void)
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
std::vector<std::string> TWeakNonEye3ptDeriv<FImpl, GImpl>::legNames(void) const
{
    return {par().qBarLeft, par().qLeft, par().qBarRight, par().qRight};
}

template <typename FImpl, typename GImpl>
std::vector<typename TWeakNonEye3ptDeriv<FImpl, GImpl>::GammaPair>
TWeakNonEye3ptDeriv<FImpl, GImpl>::parseGammaPairs(void) const
{
    const std::string      &str = par().gammaPairs;
    std::vector<GammaPair> pairs;

    if (str.empty())
    {
        HADRONS_ERROR(Argument, "gammaPairs must be set ('diag', 'all' or a list of pairs)");
    }
    else if (str == "diag")
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
        if (pairs.empty())
        {
            HADRONS_ERROR(Argument, "wrong gamma pair list '" + str + "'");
        }
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
std::vector<typename TWeakNonEye3ptDeriv<FImpl, GImpl>::DerivSpec>
TWeakNonEye3ptDeriv<FImpl, GImpl>::parseDerivs(void) const
{
    std::vector<DerivSpec> specs;

    if (par().derivatives.empty())
    {
        HADRONS_ERROR(Argument, "derivatives must contain at least one specification "
                      "(legs = 'none' for the underived contraction)");
    }
    for (auto &d: par().derivatives)
    {
        std::vector<std::string> legs = strToVec<std::string>(d.legs);
        DerivSpec                s;

        s.label  = "none";
        s.couple = 0;
        for (unsigned int k = 0; k < 4; ++k)
        {
            s.leg[k] = false;
        }
        if (legs.empty())
        {
            HADRONS_ERROR(Argument, "derivative legs must be set ('none' or a list of legs)");
        }
        if (d.legs == "none")
        {
            if (!d.dirs.empty() or !d.couple.empty())
            {
                HADRONS_ERROR(Argument, "dirs and couple must be empty when legs = 'none'");
            }
            specs.push_back(s);
            continue;
        }
        for (auto &l: legs)
        {
            auto it = std::find(weakDerivLegs.begin(), weakDerivLegs.end(), l);

            if (it == weakDerivLegs.end())
            {
                HADRONS_ERROR(Argument, "unknown leg '" + l + "' (expected qBarLeft, "
                              "qLeft, qBarRight or qRight)");
            }
            s.leg[it - weakDerivLegs.begin()] = true;
        }
        s.dirs = strToVec<int>(d.dirs);
        if (s.dirs.empty())
        {
            HADRONS_ERROR(Argument, "derivative dirs must be set for legs '" + d.legs + "'");
        }
        for (int rho: s.dirs)
        {
            if ((rho < 0) or (rho >= Nd))
            {
                HADRONS_ERROR(Argument, "direction out of range: " + d.dirs);
            }
        }
        if (d.couple == "none")
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
            HADRONS_ERROR(Argument, "couple must be 'none', 'G1' or 'G2'");
        }

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
        specs.push_back(s);
    }

    return specs;
}

// add the momentum projection of c to acc[momentum][time]
template <typename FImpl, typename GImpl>
void TWeakNonEye3ptDeriv<FImpl, GImpl>::project(MomCorr &acc, ComplexField &c)
{
    const unsigned int                                nt = env().getDim(Tp);
    std::vector<typename ComplexField::scalar_object> modes;

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
void TWeakNonEye3ptDeriv<FImpl, GImpl>::setup(void)
{
    GridBase *grid = envGetGrid(PropagatorField);
    auto     names = legNames();

    for (unsigned int k = 0; k < 4; ++k)
    {
        if (names[k].empty())
        {
            HADRONS_ERROR(Argument, weakDerivLegs[k] + " must be set");
        }
    }
    if ((par().gammaIn == Gamma::Algebra::undef) or (par().gammaOut == Gamma::Algebra::undef))
    {
        HADRONS_ERROR(Argument, "gammaIn and gammaOut must be set");
    }
    if (par().gauge.empty())
    {
        HADRONS_ERROR(Argument, "gauge must be set");
    }
    if ((par().bc != 1) and (par().bc != -1))
    {
        HADRONS_ERROR(Argument, "bc must be +1 or -1");
    }
    parseGammaPairs();
    zeroMomentumOnly_ = par().momProjector.empty();

    // distinct differentiated legs
    derivNames_.clear();
    for (auto &sp: parseDerivs())
    {
        for (unsigned int k = 0; k < 4; ++k)
        {
            if (sp.leg[k] and (std::find(derivNames_.begin(), derivNames_.end(), names[k])
                               == derivNames_.end()))
            {
                derivNames_.push_back(names[k]);
            }
        }
    }

    envTmpLat(PropagatorField, "L");
    envTmpLat(PropagatorField, "R");
    envTmpLat(GaugeLinkField, "link");
    if (!derivNames_.empty())
    {
        envTmp(std::vector<PropagatorField>, "dq", 1, derivNames_.size(), grid);
    }
    envTmpLat(ComplexField, "contraction");
    envCreate(HadronsSerializable, getName(), 1, 0);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl, typename GImpl>
void TWeakNonEye3ptDeriv<FImpl, GImpl>::execute(void)
{
    LOG(Message) << "Computing mesonic weak 3pt contractions, non-eye topologies, "
                 << "with gamma pairs and derivatives" << std::endl;
    LOG(Message) << "gIn : " << par().gammaIn << std::endl;
    LOG(Message) << "gOut: " << par().gammaOut << std::endl;
    LOG(Message) << "ql  : " << par().qLeft << std::endl;
    LOG(Message) << "qbl : " << par().qBarLeft << std::endl;
    LOG(Message) << "qr  : " << par().qRight << std::endl;
    LOG(Message) << "qbr : " << par().qBarRight << std::endl;

    auto                         &U    = envGet(GaugeField, par().gauge);
    auto                         names = legNames();
    std::vector<GammaPair>       pairs = parseGammaPairs();
    std::vector<DerivSpec>       specs = parseDerivs();
    Gamma                        gIn(par().gammaIn), gOut(par().gammaOut);
    Gamma                        g5(Gamma::Algebra::Gamma5);
    std::array<const PropagatorField *, 4> legs;
    Result                       result;

    envGetTmp(PropagatorField, L);
    envGetTmp(PropagatorField, R);
    envGetTmp(GaugeLinkField, link);
    envGetTmp(ComplexField, contraction);

    result.in  = par().gammaIn;
    result.out = par().gammaOut;
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
    for (unsigned int k = 0; k < 4; ++k)
    {
        legs[k] = &envGet(PropagatorField, names[k]);
    }
    LOG(Message) << pairs.size() << " gamma pair(s), " << specs.size()
                 << " derivative specification(s)" << std::endl;

    // accumulators acc[spec][pair][trace][momentum][time]
    std::vector<std::vector<std::array<MomCorr, 2>>> acc(specs.size(),
        std::vector<std::array<MomCorr, 2>>(pairs.size()));

    // contract one spec with the given legs, rho < 0: no coupling
    auto contract = [&](const unsigned int s,
                        const std::array<const PropagatorField *, 4> &q,
                        const int rho)
    {
        const DerivSpec &sp = specs[s];

        startTimer("Trace");
        L = (*q[1])*adj(gIn)*g5*adj(*q[0])*g5;
        R = (*q[2])*gOut*g5*adj(*q[3])*g5;
        stopTimer("Trace");
        for (unsigned int p = 0; p < pairs.size(); ++p)
        {
            Gamma G1(pairs[p].first), G2(pairs[p].second);

            if ((sp.couple == 1) and (rho >= 0))
            {
                G1 = Gamma::gmu[rho]*G1;
            }
            else if ((sp.couple == 2) and (rho >= 0))
            {
                G2 = Gamma::gmu[rho]*G2;
            }
            startTimer("Trace");
            contraction = trace(L*G1*R*G2);
            stopTimer("Trace");
            project(acc[s][p][0], contraction);
            startTimer("Trace");
            contraction = trace(L*G1)*trace(R*G2);
            stopTimer("Trace");
            project(acc[s][p][1], contraction);
        }
    };

    // underived contractions
    for (unsigned int s = 0; s < specs.size(); ++s)
    {
        if (specs[s].dirs.empty())
        {
            contract(s, legs, -1);
        }
    }

    // derived contractions, one direction at a time, each leg derived once
    if (!derivNames_.empty())
    {
        envGetTmp(std::vector<PropagatorField>, dq);
        for (int rho = 0; rho < Nd; ++rho)
        {
            bool needed = false;

            for (auto &sp: specs)
            {
                needed = needed or (std::find(sp.dirs.begin(), sp.dirs.end(), rho)
                                    != sp.dirs.end());
            }
            if (!needed)
            {
                continue;
            }
            startTimer("Insertion");
            boundaryLink<GImpl>(link, U, rho, par().bc);
            for (unsigned int n = 0; n < derivNames_.size(); ++n)
            {
                auto &q = envGet(PropagatorField, derivNames_[n]);

                covDeriv<GImpl>(dq[n], q, link, rho, DerivType::symmetric);
            }
            stopTimer("Insertion");
            for (unsigned int s = 0; s < specs.size(); ++s)
            {
                auto &sp = specs[s];

                if (std::find(sp.dirs.begin(), sp.dirs.end(), rho) == sp.dirs.end())
                {
                    continue;
                }

                std::array<const PropagatorField *, 4> q = legs;

                for (unsigned int k = 0; k < 4; ++k)
                {
                    if (sp.leg[k])
                    {
                        auto n = std::find(derivNames_.begin(), derivNames_.end(),
                                           names[k]) - derivNames_.begin();

                        q[k] = &dq[n];
                    }
                }
                contract(s, q, rho);
            }
        }
    }

    // store
    for (unsigned int s = 0; s < specs.size(); ++s)
    for (unsigned int p = 0; p < pairs.size(); ++p)
    for (unsigned int tr = 0; tr < 2; ++tr)
    {
        Contraction c;

        c.deriv  = specs[s].label;
        c.gamma1 = pairs[p].first;
        c.gamma2 = pairs[p].second;
        c.trace  = tr + 1;
        c.corr   = acc[s][p][tr];
        result.contractions.push_back(c);
    }
    startTimer("I/O");
    saveResult(par().output, "weakNonEye3ptDeriv", result);
    auto &out = envGet(HadronsSerializable, getName());
    out = result;
    stopTimer("I/O");
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_WeakNonEye3ptDeriv_hpp_
