/*
 * StochasticMeson.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
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

#ifndef Hadrons_MContraction_StochasticMeson_hpp_
#define Hadrons_MContraction_StochasticMeson_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Serialization.hpp>

BEGIN_HADRONS_NAMESPACE

/*
 * Stochastic meson correlator 
 * -----------------------------------------------------
 * Computes connected two-point function with bilinears:
 * 
 *   C_{Gamma1 Gamma2}(t) = sum_{x,y} tr[
 *       (gamma5 Gamma1) * (phi_1A * eta_A^dag)(t,x; t0,y) *
 *       (Gamma2^dag gamma5) * (eta_B * phi_2B^dag)(t0,y; t,x)
 *   ]
 * 
 * The reverse-oriented second line is represented with gamma5 hermiticity.
 * For independent wall-noise fields eta_A and eta_B, the local factors are
 *   P(y) = eta_A(y)^dag Gamma_source^dag gamma5 eta_B(y),
 *   Q(x) = phi_2B(x)^dag gamma5 Gamma_sink phi_1A(x).
 * This keeps both noise factors at the source and both solution fields at
 * the sink, so the estimator has a nonzero zero-flow wall-source limit.
 * P is summed over space on the explicitly supplied source time slice; Q is
 * projected by the supplied MSink object and stored at absolute sink time,
 * following Meson.hpp.
 * 
 * Parameters:
 * - phi1A: solution field for flavour 1, noise A
 * - etaA:  noise field for noise A
 * - phi2B: solution field for flavour 2, noise B
 * - etaB:  noise field for noise B
 * - gammas: gamma matrix pairs "(Gamma_sink Gamma_source)..."
 * - sink:   spatial momentum projector (MSink module name)
 * - sourceTime: centre time of the source bilinear
 * - noiseA/noiseB: labels from one common per-source-time noise pool
 * - flowTime: fermion flow time recorded with the result
 * - output: output file stem
 * 
 * The module computes one noise pair (A,B) per invocation.
 * Independent ensembles A and B are required (etaA != etaB).
 */

BEGIN_MODULE_NAMESPACE(MContraction)

typedef std::pair<Gamma::Algebra, std::vector<Gamma::Algebra>> GammaMapEntry;

class StochasticMesonPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(StochasticMesonPar,
                                    std::string, phi1A,
                                    std::string, etaA,
                                    std::string, phi2B,
                                    std::string, etaB,
                                    std::string, gammas,
                                    std::string, sink,
                                    unsigned int, sourceTime,
                                    unsigned int, noiseA,
                                    unsigned int, noiseB,
                                    double, flowTime,
                                    std::string, output);
};

template <typename FImpl1, typename FImpl2>
class TStochasticMeson: public Module<StochasticMesonPar>
{
public:
    FERM_TYPE_ALIASES(FImpl1, 1);
    FERM_TYPE_ALIASES(FImpl2, 2);
    BASIC_TYPE_ALIASES(ScalarImplCR, Scalar);
    SINK_TYPE_ALIASES(Scalar);
    
    class Result: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(Result,
                                        Gamma::Algebra, gamma_snk,
                                        Gamma::Algebra, gamma_src,
                                        double, flow_time,
                                        unsigned int, source_time,
                                        unsigned int, noise_a,
                                        unsigned int, noise_b,
                                        std::vector<Complex>, corr);
    };
    
public:
    // constructor
    TStochasticMeson(const std::string name);
    // destructor
    virtual ~TStochasticMeson(void) {};
    // parse gamma string
    virtual int parseGammaString(std::map<Gamma::Algebra, std::vector<Gamma::Algebra>> &gammaMap);
    // dependencies/products
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    virtual std::vector<std::string> getOutputFiles(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
    
};

MODULE_REGISTER_TMP(StochasticMeson, ARG(TStochasticMeson<FIMPL, FIMPL>), MContraction);

/******************************************************************************
 *                        TStochasticMeson implementation                     *
 ******************************************************************************/

// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2>
TStochasticMeson<FImpl1, FImpl2>::TStochasticMeson(const std::string name)
: Module<StochasticMesonPar>(name)
{}

// parse arguments /////////////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2>
int TStochasticMeson<FImpl1, FImpl2>::parseGammaString(std::map<Gamma::Algebra, std::vector<Gamma::Algebra>> &gammaMap)
{
    int nGammaPairs = 0;
    
    if (par().gammas.compare("all") == 0)
    {
        for (Gamma gammaA: Gamma::gall) {
            std::vector<Gamma::Algebra> gMapTmp;
            for (Gamma gammaB: Gamma::gall) {      
                gMapTmp.push_back(gammaB.g);
                nGammaPairs++;
            }
            gammaMap[gammaA.g] = gMapTmp;
        }
    }
    else
    {
        std::vector<std::pair<Gamma::Algebra, Gamma::Algebra>> tmp;
        tmp = strToVec<std::pair<Gamma::Algebra, Gamma::Algebra>>(par().gammas);
        
        for (unsigned int j = 0; j < tmp.size(); j++)
        {
            if ((tmp[j].first == Gamma::Algebra::undef) || 
                (tmp[j].second == Gamma::Algebra::undef))
            {
                HADRONS_ERROR(Argument, "Wrong Argument for Gamma matrices. " + par().gammas); 
            }
        }
        for (Gamma gammaA: Gamma::gall) {
            std::vector<Gamma::Algebra> gMapTmp;
            for (unsigned int j = 0; j < tmp.size(); j++)
            {
                if (tmp[j].first == gammaA.g)
                {
                    gMapTmp.push_back(tmp[j].second);
                    nGammaPairs++;
                }
            }
            gammaMap[gammaA.g] = gMapTmp;
        }
    }
    
    return nGammaPairs; 
}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2>
std::vector<std::string> TStochasticMeson<FImpl1, FImpl2>::getInput(void)
{
    std::vector<std::string> input = {par().phi1A, par().etaA, par().phi2B, par().etaB};
    
    if (!par().sink.empty())
    {
        input.push_back(par().sink);
    }
    
    return input;
}

template <typename FImpl1, typename FImpl2>
std::vector<std::string> TStochasticMeson<FImpl1, FImpl2>::getOutput(void)
{
    std::vector<std::string> output = {getName()};
    
    return output;
}

template <typename FImpl1, typename FImpl2>
std::vector<std::string> TStochasticMeson<FImpl1, FImpl2>::getOutputFiles(void)
{
    std::vector<std::string> output;
    
    if (!par().output.empty())
        output.push_back(resultFilename(par().output));
    
    return output;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2>
void TStochasticMeson<FImpl1, FImpl2>::setup(void)
{
    envTmpLat(LatticeComplex, "c");
    envCreate(HadronsSerializable, getName(), 1, 0);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2>
void TStochasticMeson<FImpl1, FImpl2>::execute(void)
{
    if (par().etaA == par().etaB)
    {
        HADRONS_ERROR(Argument, "etaA and etaB must be distinct noise fields");
    }

    LOG(Message) << "Computing stochastic meson contraction '" << getName() << "' using"
                 << " phi1A='" << par().phi1A << "', etaA='" << par().etaA << "'"
                 << " phi2B='" << par().phi2B << "', etaB='" << par().etaB << "'"
                 << std::endl;

    std::vector<Result> result;
    Gamma g5(Gamma::Algebra::Gamma5);
    int nt = env().getDim(Tp);
    if (par().sourceTime >= static_cast<unsigned int>(nt))
    {
        HADRONS_ERROR(Argument, "sourceTime is outside the temporal lattice extent");
    }
    LOG(Message) << ", source_time=" << par().sourceTime
                 << ", noise_a=" << par().noiseA
                 << ", noise_b=" << par().noiseB
                 << std::endl;

    std::map<Gamma::Algebra, std::vector<Gamma::Algebra>> gammaMap;
    int nGammas = parseGammaString(gammaMap);
    result.resize(nGammas);
    for (unsigned int i = 0; i < result.size(); ++i)
    {
        result[i].corr.resize(nt);
        result[i].source_time = par().sourceTime;
        result[i].noise_a = par().noiseA;
        result[i].noise_b = par().noiseB;
        result[i].flow_time = par().flowTime;
    }
    
    // Get input fields
    auto &phi1A = envGet(PropagatorField1, par().phi1A);
    auto &etaA  = envGet(PropagatorField1, par().etaA);
    auto &phi2B = envGet(PropagatorField2, par().phi2B);
    auto &etaB  = envGet(PropagatorField2, par().etaB);
    
    if (par().sink.empty())
    {
        HADRONS_ERROR(Definition, "no sink provided");
    }

    LOG(Message) << "(using sink '" << par().sink << "')" << std::endl;

    std::string sinkNs = vm().getModuleNamespace(env().getObjectModule(par().sink));
    if (sinkNs != "MSink")
    {
        HADRONS_ERROR(Definition, "sink must be an MSink module");
    }
    SinkFnScalar &sink = envGet(SinkFnScalar, par().sink);

    // Use global lattice coordinates, rather than local grid indices, so the
    // source-time restriction and sink slicing work with a distributed time direction.
    Lattice<iScalar<vInteger>> time(phi1A.Grid());
    LatticeCoordinate(time, Tp);
    PropagatorField1 source(phi1A.Grid());
    envGetTmp(LatticeComplex, c);

    unsigned int i = 0;
    for (auto &ss : gammaMap)
    {
        Gamma::Algebra gammaSink = ss.first;
        Gamma gSnk(gammaSink);

        for (Gamma::Algebra &gammaSource : ss.second)
        {
            Gamma gSrc(gammaSource);

            // This follows Meson.hpp: the source bilinear is Hermitian
            // conjugated, while the sink gamma is not.
            // P = sum_{y, y0=sourceTime} etaA^dag adj(Gamma_source) gamma5 etaB.
            source = adj(etaA) * adj(gSrc) * g5 * etaB;
            // This locates the flowed source operator at sourceTime; it does
            // not restrict the temporal smearing already contained in its fields.
            source = where(time == static_cast<int>(par().sourceTime), source, 0.*source);
            auto P = sum(source);

            // Q = phi2B^dag g5 Gamma_sink phi1A.  Gamma5 hermiticity puts
            // both wall-noise factors at the source rather than at the sink.
            c = trace(adj(phi2B) * g5 * gSnk * phi1A * P);
            std::vector<TComplex> buf = sink(c);
            for (int dt = 0; dt < nt; ++dt)
            {
                result[i].corr[dt] = TensorRemove(buf[dt]);
            }
            
            result[i].gamma_snk = gSnk.g;
            result[i].gamma_src = gSrc.g;
            i++;
        }
    }
    
    startTimer("I/O");
    saveResult(par().output, "stochastic_meson", result);
    stopTimer("I/O");
    
    auto &out = envGet(HadronsSerializable, getName());
    out = result;
    
    LOG(Message) << "Stochastic meson contraction complete." << std::endl;
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_StochasticMeson_hpp_
