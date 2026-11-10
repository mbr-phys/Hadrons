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
#include <Grid/algorithms/blas/MomentumProject.h>
#include <Hadrons/MomentumProjectMany.hpp>

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
 * Momentum. The source is projected onto the momentum opposite to the sink
 * one: with sink phase exp(+i p.x) the source carries exp(-i p.y), so that
 * the result is sum_{x,y} exp(i p.(x - y)) C(x, y), the translation average
 * of the point-to-all correlator at momentum p, times the spatial volume.
 * - With momProjector, every momentum of the projector is computed at once,
 *   the source phase being the complex conjugate of the projector phase.
 *   This path needs spin-colour diagonal noise, eta = 1_spin x M with M a
 *   colour matrix field (MNoise::*SpinColorDiagonal, also after fermion flow,
 *   which acts trivially in spin); this is checked. Then
 *     P_p = (Gamma_source^dag gamma5) x N_p,
 *     N_p = sum_{y, y0 = sourceTime} exp(-i p.y) M_A(y)^dag M_B(y),
 *   is independent of the gamma structures, and with the colour matrix field
 *     Y(x) = tr_spin[phi_2B(x)^dag gamma5 Gamma_sink phi_1A(x) Gamma_source^dag gamma5]
 *   the correlator is C(p, t) = tr_colour[Y(p, t) N_p]. The Y fields of all
 *   gamma pairs, and N, are each projected with one GEMM
 *   (MomentumProjectMany).
 * - With sink, the source is projected onto -p for sourceMom = p (empty:
 *   zero momentum), matching an MSink::Point sink at momentum p.
 * 
 * Parameters:
 * - phiA: solution field for flavour 1, noise A
 * - etaA:  noise field for noise A
 * - phiB: solution field for flavour 2, noise B
 * - etaB:  noise field for noise B
 * - gammas: gamma matrix pairs "(Gamma_sink Gamma_source)..."
 * - sink:   spatial momentum projector (MSink module name)
 * - momProjector: momentum projector module (MUtilities::MomentumProjectorList
 *                 or MPScalar), used instead of sink to project onto all its
 *                 momenta at once. Exclusive with sink. The result then holds
 *                 the momenta (integer vectors) and corr[momentum][time].
 * - sourceMom: source momentum (integer vector) for the sink path, empty for
 *              zero; not used with momProjector
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
                                    std::string, phiA,
                                    std::string, etaA,
                                    std::string, phiB,
                                    std::string, etaB,
                                    std::string, gammas,
                                    std::string, sink,
                                    unsigned int, sourceTime,
                                    unsigned int, noiseA,
                                    unsigned int, noiseB,
                                    double, flowTime,
                                    std::string, output,
                                    unsigned int, detSrc,
                                    std::string, momProjector,
                                    std::string, sourceMom);
};

template <typename FImpl1, typename FImpl2>
class TStochasticMeson: public Module<StochasticMesonPar>
{
public:
    FERM_TYPE_ALIASES(FImpl1, 1);
    FERM_TYPE_ALIASES(FImpl2, 2);
    BASIC_TYPE_ALIASES(ScalarImplCR, Scalar);
    SINK_TYPE_ALIASES(Scalar);
    using MPManyType = MomentumProjectMany<ComplexField1, ComplexField1>;
    typedef Lattice<iScalar<iScalar<iMatrix<typename ComplexField1::vector_type, Nc>>>> ColourMatrixField1;
    typedef typename ColourMatrixField1::scalar_object                                  ColourMatrix1;
    
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
    class MomResult: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(MomResult,
                                        Gamma::Algebra, gamma_snk,
                                        Gamma::Algebra, gamma_src,
                                        double, flow_time,
                                        unsigned int, source_time,
                                        unsigned int, noise_a,
                                        unsigned int, noise_b,
                                        std::vector<std::vector<int>>, momentum,
                                        std::vector<std::vector<Complex>>, corr);
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
private:
    // M with eta = 1_spin x M, checking that eta has this form
    template <typename PropField>
    void spinTrivialColour(ColourMatrixField1 &M, const PropField &eta,
                           PropField &tmp, const std::string &name);
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
    std::vector<std::string> input = {par().phiA, par().etaA, par().phiB, par().etaB};
    
    if (!par().sink.empty())
    {
        input.push_back(par().sink);
    }
    if (!par().momProjector.empty())
    {
        input.push_back(par().momProjector);
        input.push_back(par().momProjector + "_momList");
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
    if (!par().sink.empty() and !par().momProjector.empty())
    {
        HADRONS_ERROR(Argument, "sink and momProjector cannot both be set");
    }
    if (!par().sourceMom.empty() and par().sink.empty())
    {
        HADRONS_ERROR(Argument, "sourceMom is only used with sink");
    }
    envTmpLat(LatticeComplex, "c");
    envTmpLat(PropagatorField1, "sinkProd");
    envTmpLat(PropagatorField1, "prod");
    if (!par().momProjector.empty())
    {
        // one Y field per gamma pair, and the source colour matrix field N,
        // each set projected together
        std::map<Gamma::Algebra, std::vector<Gamma::Algebra>> gammaMap;

        envTmp(std::vector<ColourMatrixField1>, "traces", 1, parseGammaString(gammaMap),
               envGetGrid(ComplexField1));
        envTmp(std::vector<ColourMatrixField1>, "sourceCol", 1, 1, envGetGrid(ComplexField1));
        envTmpLat(ColourMatrixField1, "colA");
        envTmpLat(ColourMatrixField1, "colB");
    }
    else if (!par().sourceMom.empty())
    {
        envTmpLat(LatticeComplex, "coor");
        envTmpLat(LatticeComplex, "sourcePhase");
    }
    envCreate(HadronsSerializable, getName(), 1, 0);
}

// spin structure of the noise /////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2>
template <typename PropField>
void TStochasticMeson<FImpl1, FImpl2>::spinTrivialColour(ColourMatrixField1 &M,
                                                         const PropField &eta,
                                                         PropField &tmp,
                                                         const std::string &name)
{
    M   = peekSpin(eta, 0, 0);
    tmp = Zero();
    for (int s = 0; s < Ns; ++s)
    {
        pokeSpin(tmp, M, s, s);
    }
    tmp = tmp - eta;

    const RealD defect = norm2(tmp), scale = norm2(eta);

    if (!(defect <= 1.0e-20*scale))
    {
        HADRONS_ERROR(Argument, "noise field '" + name + "' is not of the form "
                      "1_spin x colour matrix (|defect|^2/|eta|^2 = "
                      + std::to_string(defect/scale) + "); momProjector needs "
                      "spin-colour diagonal noise");
    }
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2>
void TStochasticMeson<FImpl1, FImpl2>::execute(void)
{
    if (par().detSrc == 1)
    {
        LOG(Message) << "Allowing identical sources - these should be deterministic or you're making an error." << std::endl;
    }

    if ((par().etaA == par().etaB) && (par().detSrc == 0))
    {
        HADRONS_ERROR(Argument, "detSrc = 0 : etaA and etaB must be distinct noise fields");
    }

    LOG(Message) << "Computing stochastic meson contraction '" << getName() << "' using"
                 << " phiA='" << par().phiA << "', etaA='" << par().etaA << "'"
                 << " phiB='" << par().phiB << "', etaB='" << par().etaB << "'"
                 << std::endl;

    std::vector<Result>    result;
    std::vector<MomResult> momResult;
    const bool             projected = !par().momProjector.empty();
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
    auto &phiA = envGet(PropagatorField1, par().phiA);
    auto &etaA  = envGet(PropagatorField1, par().etaA);
    auto &phiB = envGet(PropagatorField2, par().phiB);
    auto &etaB  = envGet(PropagatorField2, par().etaB);

    envGetTmp(PropagatorField1, sinkProd);
    envGetTmp(PropagatorField1, prod);

    if (projected)
    {
        LOG(Message) << "(using momentum projector '" << par().momProjector
                     << "')" << std::endl;

        auto &momList = envGet(std::vector<std::vector<int>>, par().momProjector + "_momList");
        auto &mp      = envGet(MPManyType, par().momProjector);
        const unsigned int nmom = momList.size();

        envGetTmp(std::vector<ColourMatrixField1>, traces);
        envGetTmp(std::vector<ColourMatrixField1>, sourceCol);
        envGetTmp(ColourMatrixField1, colA);
        envGetTmp(ColourMatrixField1, colB);

        // source: N_p = sum_y exp(-i p.y) M_A(y)^dag M_B(y) on the source time
        // slice, from the conjugate of the projection of conj(M_A^dag M_B)
        std::vector<std::vector<ColourMatrix1>> projSource;
        std::vector<ColourMatrix1>              Np(nmom);

        startTimer("source");
        spinTrivialColour(colA, etaA, prod, par().etaA);
        spinTrivialColour(colB, etaB, prod, par().etaB);
        sourceCol[0] = adj(colA)*colB;
        sourceCol[0] = conjugate(sourceCol[0]);
        stopTimer("source");
        startTimer("momentumProject");
        mp.ProjectMany(sourceCol, projSource);
        stopTimer("momentumProject");
        for (unsigned int m = 0; m < nmom; ++m)
        {
            Np[m] = conjugate(projSource[0][par().sourceTime + nt*m]);
        }

        // sink: Y(x) = tr_spin[phiB^dag g5 gSnk phiA adj(gSrc) g5] per gamma pair
        momResult.resize(nGammas);
        unsigned int i = 0;
        for (auto &ss : gammaMap)
        {
            Gamma gSnk(ss.first);

            if (ss.second.empty())
            {
                continue;
            }
            startTimer("sink");
            sinkProd = adj(phiB)*g5*gSnk*phiA;
            stopTimer("sink");
            for (Gamma::Algebra &gammaSource : ss.second)
            {
                Gamma gSrc(gammaSource);

                startTimer("sink");
                prod      = sinkProd*(adj(gSrc)*g5);
                traces[i] = traceSpin(prod);
                stopTimer("sink");
                momResult[i].gamma_snk = gSnk.g;
                momResult[i].gamma_src = gSrc.g;
                i++;
            }
        }

        std::vector<std::vector<ColourMatrix1>> projSink;

        startTimer("momentumProject");
        mp.ProjectMany(traces, projSink);
        stopTimer("momentumProject");
        for (unsigned int g = 0; g < momResult.size(); ++g)
        {
            auto &r = momResult[g];

            r.source_time = par().sourceTime;
            r.noise_a     = par().noiseA;
            r.noise_b     = par().noiseB;
            r.flow_time   = par().flowTime;
            r.momentum    = momList;
            r.corr.resize(nmom, std::vector<Complex>(nt));
            for (unsigned int m = 0; m < nmom; ++m)
            for (int t = 0; t < nt; ++t)
            {
                r.corr[m][t] = TensorRemove(trace(projSink[g][t + nt*m]*Np[m]));
            }
        }
    }
    else
    {
        if (par().sink.empty())
        {
            HADRONS_ERROR(Definition, "no sink or momProjector provided");
        }

        LOG(Message) << "(using sink '" << par().sink << "')" << std::endl;

        std::string sinkNs = vm().getModuleNamespace(env().getObjectModule(par().sink));
        if (sinkNs != "MSink")
        {
            HADRONS_ERROR(Definition, "sink must be an MSink module");
        }

        SinkFnScalar &sink   = envGet(SinkFnScalar, par().sink);
        const bool   srcMom  = !par().sourceMom.empty();
        envGetTmp(LatticeComplex, c);

        // source phase exp(-i p.y), opposite to the MSink::Point phase
        if (srcMom)
        {
            envGetTmp(LatticeComplex, coor);
            envGetTmp(LatticeComplex, sourcePhase);
            const Complex     im(0., 1.);
            std::vector<Real> p = strToVec<Real>(par().sourceMom);

            if (p.size() != static_cast<unsigned int>(env().getNd() - 1))
            {
                HADRONS_ERROR(Size, "sourceMom must have " + std::to_string(env().getNd() - 1)
                              + " components");
            }
            sourcePhase = Zero();
            for (unsigned int mu = 0; mu < p.size(); ++mu)
            {
                LatticeCoordinate(coor, mu);
                sourcePhase = sourcePhase + (p[mu]/env().getDim(mu))*coor;
            }
            sourcePhase = exp(static_cast<Real>(-2.*M_PI)*im*sourcePhase);
            LOG(Message) << "(source momentum " << par().sourceMom << ")" << std::endl;
        }

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
                // P = sum_{y, y0=sourceTime} [exp(-i p.y)] etaA^dag adj(Gamma_source) gamma5 etaB.
                prod = adj(etaA) * adj(gSrc) * g5 * etaB;
                if (srcMom)
                {
                    envGetTmp(LatticeComplex, sourcePhase);
                    prod = sourcePhase*prod;
                }
                // This locates the flowed source operator at sourceTime; it does
                // not restrict the temporal smearing already contained in its fields.
                auto P = sliceSum(prod, Tp)[par().sourceTime];

                // Q = phiB^dag g5 Gamma_sink phiA.  Gamma5 hermiticity puts
                // both wall-noise factors at the source rather than at the sink.
                c = trace(adj(phiB) * g5 * gSnk * phiA * P);

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
    }

    startTimer("I/O");
    auto &out = envGet(HadronsSerializable, getName());
    if (projected)
    {
        saveResult(par().output, "stochastic_meson", momResult);
        out = momResult;
    }
    else
    {
        saveResult(par().output, "stochastic_meson", result);
        out = result;
    }
    stopTimer("I/O");

    LOG(Message) << "Stochastic meson contraction complete." << std::endl;
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_StochasticMeson_hpp_
