#ifndef Hadrons_MDistil_UnsmearedMeson_hpp_
#define Hadrons_MDistil_UnsmearedMeson_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/TimerArray.hpp>
#include <Hadrons/Serialization.hpp>
#include <Hadrons/Modules/MDistil/Base.hpp>

BEGIN_HADRONS_NAMESPACE

/************************************************************
 *                   UnsmearedMeson                         *
 * Computes the following diagram:                          * 
 *                                                          * 
 *                                                          * 
 *           _______________                                * 
 *          /               \                               * 
 *         /                 \                              * 
 *        /                   \                             * 
 *  M(rho1,rho2)            M(phi1,phi2)                    * 
 *        \                   /                             * 
 *         \                 /                              * 
 *          \_______________/                               * 
 *                                                          * 
 *                                                          * 
 *      tSrc                   tSnk                         * 
 *                                                          * 
 *                                                          * 
 ************************************************************/
BEGIN_MODULE_NAMESPACE(MDistil)

class UnsmearedMesonPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(UnsmearedMesonPar,
                                    std::string,               output,          // file stem for the out file
                                    std::string,               RhoRhoStem,      // file stem for the rho-rho MFs
                                    std::string,               RhoRhoField,     // M(rho,rho) meson field for the Kpi
                                    std::string,               vectorStemC,     // charm
                                    std::string,               vectorStemL,     // light
                                    int,                       readRaw3DField,  // read raw perambulator output (1) or converted files (0) 
                                    std::string,               rawTimeSources,  // time sources of the raw data
                                    std::string,               noisePol,        // noise policy of v1 - assert compatibility with M(rho,rho) 
                                    std::vector<unsigned int>, tSrcs,           // source times
                                    std::string,               gammas,          // list of gamma matrices to go at sink
                                    std::vector<std::string>,  moms,            // list of momenta 
                                   );
};

template <typename FImpl>
class TUnsmearedMeson: public Module<UnsmearedMesonPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    class Result: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(Result,
                                        std::string,          gammaSrc,
                                        std::string,          gammaSnk,
                                        std::string,          momSrc,
                                        std::string,          momSnk,
                                        unsigned int,         tSrc,
                                        std::vector<Complex>, corr);
    };
    // constructor
    TUnsmearedMeson(const std::string name);
    // destructor
    virtual ~TUnsmearedMeson(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(UnsmearedMeson, TUnsmearedMeson<FIMPL>, MDistil);

/******************************************************************************
 *                 TUnsmearedMeson implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TUnsmearedMeson<FImpl>::TUnsmearedMeson(const std::string name)
: Module<UnsmearedMesonPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TUnsmearedMeson<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().noisePol};
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TUnsmearedMeson<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()+"_cl",getName()+"_cc"};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TUnsmearedMeson<FImpl>::setup(void)
{
    GridCartesian * gridHD = envGetGrid(FermionField);
    GridCartesian * gridLD = envGetSliceGrid(FermionField,gridHD->Nd() -1);
    
    envTmp   (FermionField,      "fermion3dtmp1" ,1, gridLD);
    envTmp   (FermionField,      "fermion3dtmp2" ,1, gridLD);
    envTmp   (FermionField,      "fermion3dtmp3" ,1, gridLD);
    envTmp   (FermionField,      "fermion3dtmp4" ,1, gridLD);
    envTmp   (PropagatorField,   "prop3dtmp"     ,1, gridLD);
    envTmp   (PropagatorField,   "prop3dtmp1"    ,1, gridLD);
    envTmp   (ComplexField,      "MesPhi1"       ,1, gridLD);
    envTmp   (ComplexField,      "MesPhi2"       ,1, gridLD);
    envTmp   (ComplexField,      "MesPhi3"       ,1, gridLD);
    envTmpLat(ComplexField,      "coor");

    auto &dilNoise = envGet(DistillationNoise<FImpl>, par().noisePol);
    int nDL = dilNoise.dilutionSize(DistillationNoise<FImpl>::Index::l);        
    int nDS = dilNoise.dilutionSize(DistillationNoise<FImpl>::Index::s);     

    Grid::Coordinate coor  = gridHD->GlobalDimensions();
    coor[Tdir] = nDL * nDS;
    Grid::GridCartesian * gridDD = new GridCartesian(coor, gridLD->_simd_layout, gridLD->_processors, *gridHD);

    envTmp(FermionField,    "fermionDDtmp_light" ,1, gridDD);
    envTmp(FermionField,    "fermionDDtmp_charm" ,1, gridDD);

    envCreate(HadronsSerializable, getName()+"_ll", 1, 0);
    envCreate(HadronsSerializable, getName()+"_cl", 1, 0);
    envCreate(HadronsSerializable, getName()+"_cc", 1, 0);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TUnsmearedMeson<FImpl>::execute(void)
{

    if (par().vectorStemC.empty() && par().vectorStemL.empty()) {
        HADRONS_ERROR(Argument, "at least one of vectorStemC and vectorStemL must be filled"); 
    } else if (par().vectorStemC.empty() && !par().vectorStemL.empty()) {
        LOG(Message) << "vectorStemC is empty, only computing pion (ll) correlator(s)" << std::endl;
    } else if (!par().vectorStemC.empty() && par().vectorStemL.empty()) {
        LOG(Message) << "vectorStemL is empty, only computing etac (cc) correlator(s)" << std::endl;
    } else {
        LOG(Message) << "both vectorStems are given, computing pion (ll), etac (cc), and D meson (cl) correlator(s)" << std::endl;
    }

    // general grid setup
    GridCartesian * gridHD = envGetGrid(FermionField);
    GridCartesian * gridLD = envGetSliceGrid(FermionField,gridHD->Nd() -1);
    const int Ntlocal{gridHD->LocalDimensions()[Tdir]};
    const int Ntfirst{gridHD->LocalStarts()[Tdir]};
    int nT=env().getDim(Tdir);

    std::vector<std::string> momSrcs, momSnks;
    int nMoms = 0;

    if (par().moms.empty()) 
    {
        LOG(Message) << "You have not specified any momenta, so all possible combinations up to P^2 = 4 will be listed." << std::endl;
        for (int i = -2; i <= 2; i++)
        {
            for (int j = -2; j <= 2; j++)
            {
                for (int k = -2; k <= 2; k++)
                {
                    int mom1P2 = i*i + j*j + k*k;
                    if (mom1P2 <= 4)
                    {
                        std::string mom1 = std::to_string(i) + "_" + std::to_string(j) + "_" + std::to_string(k);
                        int o(-i), p(-j), q(-k);
                        int mom2P2 = o*o + p*p + q*q;
                        if ((std::abs(o) <= 2) && (std::abs(p) <= 2) && (std::abs(q) <= 2) && (mom2P2 <= 4))
                        {
                            std::string mom2 = std::to_string(o) + " " + std::to_string(p) + " " + std::to_string(q);
                            if (std::find(momSrcs.begin(), momSrcs.end(), mom1) == momSrcs.end())
                            {
                                momSrcs.push_back(mom1);
                            }
                            if (std::find(momSnks.begin(), momSnks.end(), mom2) == momSnks.end())
                            {
                                momSnks.push_back(mom2);
                            }
                            nMoms++;
                        }
                    }
                }
            }
        }
    }
    else
    {
        LOG(Message) << "Using " << par().moms.size() << " possible momenta provided." << std::endl;
        for (auto mom : par().moms)
        {
            std::vector<int> momI = strToVec<int>(mom);
            std::string mom1 = std::to_string(momI[0]) + "_" + std::to_string(momI[1]) + "_" + std::to_string(momI[2]);
            int o(-momI[0]), p(-momI[1]), q(-momI[2]);
            std::string mom2 = std::to_string(o) + " " + std::to_string(p) + " " + std::to_string(q);
            if (std::find(momSrcs.begin(), momSrcs.end(), mom1) == momSrcs.end())
            {
                momSrcs.push_back(mom1);
            }
            if (std::find(momSnks.begin(), momSnks.end(), mom2) == momSnks.end())
            {
                momSnks.push_back(mom2);
            }
            nMoms++;
        }
    }

    LOG(Message) << "WARNING: Assuming ordering s + ns*(l + nl*t) in DilutedNoise.hpp. This code will break when this changes!" << std::endl;

    std::vector<unsigned int> tSrcs = par().tSrcs;
    for(auto tSrc : tSrcs)
    {
        if(tSrc>=nT)
        {
            HADRONS_ERROR(Range, "all tSrcs must be smaller than nT");
        }
    }

    // read input Kpi fields
    std::string RhoRhoGamma = par().RhoRhoField;
    std::map<std::string, ContractionDistilMesonField<ComplexD,ComplexF>> RhoRhoMesonMFs;
    startTimer("MesonField IO");
    for (auto kmom : momSrcs)
    {
        std::string mfPath = par().RhoRhoStem + "rho-rho." + std::to_string(vm().getTrajectory()) + "/" + RhoRhoGamma + "_p" + kmom + ".h5";   
        TimerArray timer1;
        auto it = RhoRhoMesonMFs.find(RhoRhoGamma+"_p"+kmom);
        if (it == RhoRhoMesonMFs.end()) 
        {
            LOG(Message) << "reading " << mfPath << std::endl;
            std::vector<std::vector<int>> tVec;
            for(auto tSrc : tSrcs)
            {
                std::vector<int> srcDiag = {static_cast<int>(tSrc), static_cast<int>(tSrc)};
                if(std::find(tVec.begin(), tVec.end(), srcDiag) == tVec.end())
                {
                    tVec.push_back(srcDiag);
                }
            }
            RhoRhoMesonMFs.try_emplace(RhoRhoGamma+"_p"+kmom, ContractionDistilMesonField<ComplexD,ComplexF>(mfPath, nT, timer1, tVec));
        }
        else
        {
            LOG(Message) << "already read " << mfPath << std::endl;
        }
    }
    stopTimer("MesonField IO");

    // noise class -- assert they are identical and an "exact distillation" policy
    auto &dilNoise = envGet(DistillationNoise<FImpl>, par().noisePol);
    int nNoise = dilNoise.size(); 
    if(nNoise>1)
    {
        HADRONS_ERROR(Implementation, "UnsmearedMeson only implemented for exact distillation");
    }
    int nDL = dilNoise.dilutionSize(DistillationNoise<FImpl>::Index::l);        
    int nDS = dilNoise.dilutionSize(DistillationNoise<FImpl>::Index::s);        
    int nDT = dilNoise.dilutionSize(DistillationNoise<FImpl>::Index::t);        

    std::vector<Gamma::Algebra> gammas = strToVec<Gamma::Algebra>(par().gammas);

    std::vector<Result> clResults, ccResults, llResults;
    int resultSize = gammas.size()*tSrcs.size()*nMoms;
    clResults.resize(resultSize);
    ccResults.resize(resultSize);
    llResults.resize(resultSize);
    LOG(Message) << "Results objects have gammas (" << gammas.size() << ") * tSrcs (" << tSrcs.size() 
                 << ") * nMoms (" << nMoms << ") = " << resultSize << " size" << std::endl;
    for (unsigned int tSrci = 0; tSrci < tSrcs.size(); tSrci++)
    {
        unsigned int tSrc = tSrcs[tSrci];
        for (unsigned int i = 0; i < nMoms*gammas.size(); i++)
        {
            unsigned int ridx = tSrci*nMoms*gammas.size() + i;
            clResults[ridx].corr.resize(nT);
            ccResults[ridx].corr.resize(nT);
            llResults[ridx].corr.resize(nT);
        }
    }

    auto fillResultMetadata = [&]()
    {
        for (unsigned int tSrci = 0; tSrci < tSrcs.size(); tSrci++)
        {
            unsigned int tSrc = tSrcs[tSrci];
            unsigned int tdx = tSrci*nMoms*gammas.size();
            for (unsigned int ddx = 0; ddx < momSrcs.size(); ddx++)
            {
                std::string momSrc = momSrcs[ddx];
                std::string momSnk = momSnks[ddx];
                for (unsigned int sdx = 0; sdx < gammas.size(); sdx++)
                {
                    std::stringstream gamStr;
                    gamStr << gammas[sdx];

                    clResults[tdx].gammaSnk = gamStr.str();
                    clResults[tdx].gammaSrc = RhoRhoGamma;
                    clResults[tdx].momSnk   = momSnk;
                    clResults[tdx].momSrc   = momSrc;
                    clResults[tdx].tSrc     = tSrc;

                    ccResults[tdx].gammaSnk = gamStr.str();
                    ccResults[tdx].gammaSrc = RhoRhoGamma;
                    ccResults[tdx].momSnk   = momSnk;
                    ccResults[tdx].momSrc   = momSrc;
                    ccResults[tdx].tSrc     = tSrc;

                    llResults[tdx].gammaSnk = gamStr.str();
                    llResults[tdx].gammaSrc = RhoRhoGamma;
                    llResults[tdx].momSnk   = momSnk;
                    llResults[tdx].momSrc   = momSrc;
                    llResults[tdx].tSrc     = tSrc;

                    tdx++;
                }
            }
        }
    };
    fillResultMetadata();
    
    // Temporary objects
    envGetTmp(FermionField,       fermion3dtmp1);
    envGetTmp(FermionField,       fermion3dtmp2);
    envGetTmp(FermionField,       fermion3dtmp3);
    envGetTmp(FermionField,       fermion3dtmp4);
    envGetTmp(PropagatorField,    prop3dtmp);
    envGetTmp(PropagatorField,    prop3dtmp1);
    envGetTmp(ComplexField,       MesPhi1);
    envGetTmp(ComplexField,       MesPhi2);
    envGetTmp(ComplexField,       MesPhi3);

    //envGetTmp(ComplexField, coor);
    //envGetTmp(ComplexField, ph);
    //envGetTmp(ComplexField, ph3d);
        
    envGetTmp(FermionField,    fermionDDtmp_light);
    envGetTmp(FermionField,    fermionDDtmp_charm);

    std::map<std::string, DistilMesonFieldMatrix<ComplexD>> srcMatCache;
    auto getSrcMat = [&](const std::string &key, const unsigned int tSrc) -> const DistilMesonFieldMatrix<ComplexD>&
    {
        const std::string cacheKey = key + "_t" + std::to_string(tSrc);
        auto it = srcMatCache.find(cacheKey);
        if(it == srcMatCache.end())
        {
            it = srcMatCache.emplace(cacheKey, RhoRhoMesonMFs.at(key)(tSrc,tSrc,tSrc)).first;
        }
        return it->second;
    };

    const bool readRaw3DField = (par().readRaw3DField != 0);
    std::vector<int> rawTimeSources;
    if(readRaw3DField)
    {
        if(!par().rawTimeSources.empty())
        {
            rawTimeSources = strToVec<int>(par().rawTimeSources);
        }
        else
        {
            rawTimeSources.resize(nDT);
            std::iota(rawTimeSources.begin(), rawTimeSources.end(), 0);
        }
    }

    auto readConvertedDD = [&](FermionField &field, const std::string &stem, const unsigned int tSm,
                               const int tLoc)
    {
        std::string tFileName = stem;
        tFileName.append("_DD");
        tFileName.append("_tSm");
        tFileName.append(std::to_string(tSm));
        tFileName.append("_tLoc");
        tFileName.append(std::to_string(tLoc));
        DistillationVectorsIo::readComponent(field, tFileName, 1, nDL, nDS, nDT, 0, vm().getTrajectory());
    };

    auto readRawDD = [&](FermionField &field, const std::string &stem, const unsigned int tSm,
                         const int tLoc)
    {
        auto sourceIt = std::find(rawTimeSources.begin(), rawTimeSources.end(), static_cast<int>(tSm));
        if(sourceIt == rawTimeSources.end())
        {
            HADRONS_ERROR(Io, "raw 3D field stem does not contain requested source time " + std::to_string(tSm));
        }
        const int sourceOffset = sourceIt - rawTimeSources.begin();
        const int skip = sourceOffset * nDL * nDS;
        std::string filename = stem + "." + std::to_string(vm().getTrajectory()) + "/t" +
                               std::to_string(tLoc) + "_pkg.bin";

        field = Zero();
        ScidacReader reader;
        reader.open(filename);
        for(int id = 0; id < nDL * nDS; id++)
        {
            std::array<unsigned int, 3> index = dilNoise.dilutionCoordinates(id);
            const int dk = index[DistillationNoise<FImpl>::Index::l];
            const int ds = index[DistillationNoise<FImpl>::Index::s];
            const int dSolve = dilNoise.dilutionIndex(tSm, dk, ds);
            DistillationVectorsIo::pkgComponentReader(reader, fermion3dtmp3, 1, nDL, nDS, nDT, dSolve,
                                                      (id == 0) ? skip : 0);
            InsertSliceLocal(fermion3dtmp3, field, 0, id, Tdir);
        }
        reader.close();
    };

    auto readPhiDD = [&](FermionField &field, const std::string &stem, const unsigned int tSm,
                         const int tLoc)
    {
        if(readRaw3DField)
        {
            readRawDD(field, stem, tSm, tLoc);
        }
        else
        {
            readConvertedDD(field, stem, tSm, tLoc);
        }
    };

    int tSnk;
    std::string tFileName;
    for (unsigned int tSrci = 0; tSrci < tSrcs.size(); tSrci++)
    {
        unsigned int tSrc = tSrcs[tSrci];

        if(tSrc>=nT)
        {
            HADRONS_ERROR(Range, "tSrc must be smaller than nT");
        }

        for (int t = 0; t < Ntlocal; t++)
        {
            tSnk = t + Ntfirst;

            std::vector<TComplex> ccBuf, clBuf, llBuf;

            // read perambulators
            envGetTmp(FermionField,    fermionDDtmp_charm);
            if (!par().vectorStemC.empty()) {
                LOG(Message) << "Starting charm perambulator I/O for (tSrc,tSnk) = (" << tSrc << "," << tSnk << ")" << std::endl;
                startTimer("phi_c I/O");
                readPhiDD(fermionDDtmp_charm, par().vectorStemC, tSrc, tSnk);
                stopTimer("phi_c I/O");
            }

            envGetTmp(FermionField,    fermionDDtmp_light);
            if (!par().vectorStemL.empty()) {
                LOG(Message) << "Starting light perambulator I/O for (tSrc,tSnk) = (" << tSrc << "," << tSnk << ")" << std::endl;
                startTimer("phi_l I/O");
                readPhiDD(fermionDDtmp_light, par().vectorStemL, tSrc, tSnk);
                stopTimer("phi_l I/O");
            }

            unsigned int tdx = tSrci*nMoms*gammas.size();
            for (unsigned int ddx = 0; ddx < momSrcs.size(); ddx++) 
            {
                std::string momSrc = momSrcs[ddx];
                std::string momSnk = momSnks[ddx];

                const auto &MFmult = getSrcMat(RhoRhoGamma+"_p"+momSrc, tSrc);

                for (unsigned int sdx = 0; sdx < gammas.size(); sdx++)
                {
                    Gamma::Algebra gamma = gammas[sdx];
                    Gamma gam(gamma);

                    std::stringstream gamStr; gamStr << gamma;

                    MesPhi1 = Zero();
                    MesPhi2 = Zero();
                    MesPhi3 = Zero();
                    for (int id1=0; id1<nDL*nDS; id1++)
                    {
                        startTimer("ExtractSliceLocal");
                        if (!par().vectorStemC.empty()) {
                            ExtractSliceLocal(fermion3dtmp1, fermionDDtmp_charm, 0, id1, Tdir);
                        }
                        if (!par().vectorStemL.empty()) {
                            ExtractSliceLocal(fermion3dtmp4, fermionDDtmp_light, 0, id1, Tdir);
                        }
                        stopTimer("ExtractSliceLocal");
                        for (int id2=0; id2<nDL*nDS; id2++)
                        {
                            if (!par().vectorStemL.empty()) {
                                startTimer("ExtractSliceLocal");
                                ExtractSliceLocal(fermion3dtmp2, fermionDDtmp_light, 0, id2, Tdir);
                                stopTimer("ExtractSliceLocal");

                                startTimer("computation");
                                fermion3dtmp3 = gam*fermion3dtmp2;
                                prop3dtmp = outerProductC(fermion3dtmp4, fermion3dtmp3);
                                MesPhi3 += trace(prop3dtmp)*MFmult(id2,id1);

                                if (!par().vectorStemC.empty()) {
                                    prop3dtmp1 = outerProductC(fermion3dtmp1, fermion3dtmp3);
                                    MesPhi1 += trace(prop3dtmp1)*MFmult(id2,id1);
                                }
                                stopTimer("computation");
                            }

                            if (!par().vectorStemC.empty()) {
                                startTimer("ExtractSliceLocal");
                                ExtractSliceLocal(fermion3dtmp2, fermionDDtmp_charm, 0, id2, Tdir);
                                stopTimer("ExtractSliceLocal");

                                startTimer("computation");
                                fermion3dtmp3 = gam*fermion3dtmp2;
                                prop3dtmp = outerProductC(fermion3dtmp1, fermion3dtmp3);
                                MesPhi2 += trace(prop3dtmp)*MFmult(id2,id1);
                                stopTimer("computation");
                            }
                        }
                    }

                    if (!par().vectorStemL.empty()) {
                        startTimer("final contraction ll");
                        sliceSum(MesPhi3, llBuf, Tdir);

                        LOG(Message) << "Updating llResults for (tdx,tSnk) = (" << tdx << "," << tSnk << ")" << std::endl;
                        llResults[tdx].corr[tSnk] = TensorRemove(llBuf[0]);

                        stopTimer("final contraction ll");

                        if (!par().vectorStemC.empty()) {
                            startTimer("final contraction cl");
                            sliceSum(MesPhi1, clBuf, Tdir);

                            LOG(Message) << "Updating clResults for (tdx,tSnk) = (" << tdx << "," << tSnk << ")" << std::endl;
                            clResults[tdx].corr[tSnk] = TensorRemove(clBuf[0]);

                            stopTimer("final contraction cl");
                        }
                    }

                    if (!par().vectorStemC.empty()) {
                        startTimer("final contraction cc");
                        sliceSum(MesPhi2, ccBuf, Tdir);

                        LOG(Message) << "Updating ccResults for (tdx,tSnk) = (" << tdx << "," << tSnk << ")" << std::endl;
                        ccResults[tdx].corr[tSnk] = TensorRemove(ccBuf[0]);

                        stopTimer("final contraction cc");
                    }

                    tdx++;
                }
            }
        }
    }
    auto mergeResultCorr = [&](std::vector<Result> &results)
    {
        startTimer("result time merge");
        const bool sliceBoss = gridLD->IsBoss();
        for (auto &result : results)
        {
            for (auto &corr : result.corr)
            {
                if (!sliceBoss)
                {
                    corr = Complex(0.0, 0.0);
                }
            }
            if (!result.corr.empty())
            {
                gridHD->GlobalSumVector(result.corr.data(), static_cast<int>(result.corr.size()));
            }
        }
        stopTimer("result time merge");
    };
    startTimer("results gather/io");
    LOG(Message) << "Writing results to " << par().output << std::endl;
    if (!par().vectorStemL.empty()) {
        mergeResultCorr(llResults);

        saveResult(par().output+"_ll", "llMeson", llResults);
        auto &llOut = envGet(HadronsSerializable, getName()+"_ll");
        llOut = llResults;

        if (!par().vectorStemC.empty()) {
            mergeResultCorr(clResults);

            saveResult(par().output+"_cl", "clMeson", clResults);
            auto &clOut = envGet(HadronsSerializable, getName()+"_cl");
            clOut = clResults;
        }
    }

    if (!par().vectorStemC.empty()) {
        mergeResultCorr(ccResults);

        saveResult(par().output+"_cc", "ccMeson", ccResults);
        auto &ccOut = envGet(HadronsSerializable, getName()+"_cc");
        ccOut = ccResults;
    }
    stopTimer("results gather/io");
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MDistil_UnsmearedMeson_hpp_
