/*
 * WeakEye3pt.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2023
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Fionn O hOgain <fionn.o.hogain@ed.ac.uk>
 * Author: Lanny91 <andrew.lawson@gmail.com>
 * Author: Ryan Hill <rchrys.hill@gmail.com>
 * Author: ferben <ferben@debian.felix.com>
 * Author: fionnoh <fionnoh@gmail.com>
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
#ifndef Hadrons_MContraction_WeakEye3pt_hpp_
#define Hadrons_MContraction_WeakEye3pt_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Serialization.hpp>

BEGIN_HADRONS_NAMESPACE

/*
 * Weak Hamiltonian meson 3-pt diagrams, eye topologies.
 * 
 * Schematics:       loop                 |                  
 *                  /-<-¬                 |                             
 *                 /     \                |            qbl     G     qbr
 *                 \     /                |        /----<------*------<----¬         
 *            qbl   \   /    qbr          |       /          /-*-¬          \
 *       /-----<-----* *-----<----¬       |      /          / G2  \          \
 *  gIn *           G1 G2          * gOut | gIn *           \     /  loop    * gOut
 *       \                        /       |      \           \->-/          /
 *        \                      /        |       \                        /
 *         \---------->---------/         |        \----------->----------/
 *                   qs                   |                   qs
 *                                        |
 *                one trace               |                two traces
 *
 * one trace : tr(qbr*gOut*qs*adj(gIn)*g5*adj(qbl)*g5*G1*loop*G2)
 * two traces: tr(qbr*gOut*qs*adj(gIn)*g5*adj(qbl)*g5*G1)*tr(loop*G2)
 *
 * gammaPairs: operator gamma pairs G1 x G2, as in WeakNonEye3ptDeriv:
 *             "diag" (the 16 G x G, the original WeakEye3pt), "all" (256
 *             pairs) or an explicit list "(G1 G2)(G1 G2)...".
 *             Off-diagonal pairs are needed, e.g., for penguin-type
 *             operators after colour and Fierz identities.
 * Results: one entry per pair and trace (pair order, trace 1 then 2), with
 *          metadata {in, out, gamma1, gamma2, trace}.
 *
 * There are no default values: every parameter except output must be set.
 */

BEGIN_MODULE_NAMESPACE(MContraction)

class WeakEye3ptPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(WeakEye3ptPar,
                                    std::string,    qBarLeft,
                                    std::string,    qBarRight,
                                    std::string,    qSpectator,
                                    std::string,    loop,
                                    unsigned int,   tOut,
                                    Gamma::Algebra, gammaIn,
                                    Gamma::Algebra, gammaOut,
                                    std::string,    gammaPairs,
                                    std::string,    output);
};

template <typename FImpl>
class TWeakEye3pt: public Module<WeakEye3ptPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    typedef std::pair<Gamma::Algebra, Gamma::Algebra> GammaPair;
    class Metadata: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(Metadata,
                                        Gamma::Algebra, in,
                                        Gamma::Algebra, out,
                                        Gamma::Algebra, gamma1,
                                        Gamma::Algebra, gamma2,
                                        unsigned int,   trace);
    };
    typedef Correlator<Metadata> Result;
public:
    // constructor
    TWeakEye3pt(const std::string name);
    // destructor
    virtual ~TWeakEye3pt(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    virtual std::vector<std::string> getOutputFiles(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    void parseGammaPairs(void);
private:
    std::vector<GammaPair>      pairs_;
    std::vector<Gamma::Algebra> g1_, g2_;
};

MODULE_REGISTER_TMP(WeakEye3pt, TWeakEye3pt<FIMPL>, MContraction);

/******************************************************************************
 *                        TWeakEye3pt implementation                          *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TWeakEye3pt<FImpl>::TWeakEye3pt(const std::string name)
: Module<WeakEye3ptPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TWeakEye3pt<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().qBarLeft, par().qBarRight, 
                                   par().qSpectator, par().loop};
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TWeakEye3pt<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

template <typename FImpl>
std::vector<std::string> TWeakEye3pt<FImpl>::getOutputFiles(void)
{
    std::vector<std::string> output;
    
    if (!par().output.empty())
        output.push_back(resultFilename(par().output));
    
    return output;
}

// parser //////////////////////////////////////////////////////////////////////
// fill pairs_ and the distinct first (g1_) and second (g2_) gammas
template <typename FImpl>
void TWeakEye3pt<FImpl>::parseGammaPairs(void)
{
    const std::string &str = par().gammaPairs;

    pairs_.clear();
    if (str.empty())
    {
        HADRONS_ERROR(Argument, "gammaPairs must be set ('diag', 'all' or a list of pairs)");
    }
    else if (str == "diag")
    {
        for (auto &g: Gamma::gall)
        {
            pairs_.push_back({g.g, g.g});
        }
    }
    else if (str == "all")
    {
        for (auto &g1: Gamma::gall)
        for (auto &g2: Gamma::gall)
        {
            pairs_.push_back({g1.g, g2.g});
        }
    }
    else
    {
        pairs_ = strToVec<GammaPair>(str);
        if (pairs_.empty())
        {
            HADRONS_ERROR(Argument, "wrong gamma pair list '" + str + "'");
        }
        for (auto &p: pairs_)
        {
            if ((p.first == Gamma::Algebra::undef) or (p.second == Gamma::Algebra::undef))
            {
                HADRONS_ERROR(Argument, "wrong gamma pair list '" + str + "'");
            }
        }
    }
    g1_.clear();
    g2_.clear();
    for (auto &p: pairs_)
    {
        if (std::find(g1_.begin(), g1_.end(), p.first) == g1_.end())
        {
            g1_.push_back(p.first);
        }
        if (std::find(g2_.begin(), g2_.end(), p.second) == g2_.end())
        {
            g2_.push_back(p.second);
        }
    }
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TWeakEye3pt<FImpl>::setup(void)
{
    GridBase *grid = envGetGrid(PropagatorField);

    if (par().qBarLeft.empty() or par().qBarRight.empty()
        or par().qSpectator.empty() or par().loop.empty())
    {
        HADRONS_ERROR(Argument, "qBarLeft, qBarRight, qSpectator and loop must be set");
    }
    if ((par().gammaIn == Gamma::Algebra::undef) or (par().gammaOut == Gamma::Algebra::undef))
    {
        HADRONS_ERROR(Argument, "gammaIn and gammaOut must be set");
    }
    if (par().tOut >= env().getDim(Tp))
    {
        HADRONS_ERROR(Argument, "tOut out of range");
    }
    parseGammaPairs();
    envTmpLat(PropagatorField, "M");
    envTmpLat(ComplexField, "corr");
    envTmp(std::vector<ComplexField>, "trM", 1, g1_.size(), grid);
    envTmp(std::vector<ComplexField>, "trL", 1, g2_.size(), grid);
    envCreate(HadronsSerializable, getName(), 1, 0);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TWeakEye3pt<FImpl>::execute(void)
{
    LOG(Message) << "Computing mesonic weak 3pt contractions, eye topologies" << std::endl;
    LOG(Message) << "gIn : " << par().gammaIn << std::endl;
    LOG(Message) << "gOut: " << par().gammaOut << std::endl;
    LOG(Message) << "tOut: " << par().tOut << std::endl;
    LOG(Message) << "qbl : " << par().qBarLeft << std::endl;
    LOG(Message) << "qbr : " << par().qBarRight << std::endl;
    LOG(Message) << "qs  : " << par().qSpectator << std::endl;
    LOG(Message) << "loop: " << par().loop << std::endl;

    std::vector<Result> result;
    Result              r;
    auto                &qbl  = envGet(PropagatorField, par().qBarLeft);
    auto                &qbr  = envGet(PropagatorField, par().qBarRight);
    auto                &loop = envGet(PropagatorField, par().loop);
    auto                &qs   = envGet(SlicedPropagator, par().qSpectator);
    auto                qst   = qs[par().tOut];
    Gamma               gIn(par().gammaIn), gOut(par().gammaOut);
    Gamma               g5(Gamma::Algebra::Gamma5);

    envGetTmp(PropagatorField, M);
    envGetTmp(ComplexField, corr);
    envGetTmp(std::vector<ComplexField>, trM);
    envGetTmp(std::vector<ComplexField>, trL);

    // heavy line, and the single traces of the two-trace topology
    M = qbr*gOut*qst*adj(gIn)*g5*adj(qbl)*g5;
    for (unsigned int i = 0; i < g1_.size(); ++i)
    {
        Gamma G1(g1_[i]);

        trM[i] = trace(M*G1);
    }
    for (unsigned int i = 0; i < g2_.size(); ++i)
    {
        Gamma G2(g2_[i]);

        trL[i] = trace(loop*G2);
    }

    auto push = [&](void)
    {
        SlicedComplex buf;

        sliceSum(corr, buf, Tp);
        r.corr.clear();
        for (unsigned int t = 0; t < buf.size(); ++t)
        {
            r.corr.push_back(TensorRemove(buf[t]));
        }
        result.push_back(r);
    };

    r.info.in  = par().gammaIn;
    r.info.out = par().gammaOut;
    for (auto &p: pairs_)
    {
        Gamma        G1(p.first), G2(p.second);
        unsigned int i1 = std::find(g1_.begin(), g1_.end(), p.first) - g1_.begin();
        unsigned int i2 = std::find(g2_.begin(), g2_.end(), p.second) - g2_.begin();

        r.info.gamma1 = p.first;
        r.info.gamma2 = p.second;
        // one trace
        corr = trace(M*G1*loop*G2);
        r.info.trace = 1;
        push();
        // two traces
        corr = trM[i1]*trL[i2];
        r.info.trace = 2;
        push();
    }
    saveResult(par().output, "weakEye3pt", result);
    auto &out = envGet(HadronsSerializable, getName());
    out = result;
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_WeakEye3pt_hpp_
