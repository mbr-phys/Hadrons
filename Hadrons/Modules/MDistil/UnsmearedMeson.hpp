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
                                    std::string,               output,        // file stem for the out file
                                    std::string,               RhoRhoStem,    // file stem for the rho-rho MFs
                                    std::string,               RhoRhoField,   // M(rho,rho) meson field for the Kpi
                                    std::string,               vectorStemC,   // charm
                                    std::string,               vectorStemL,   // light
                                    std::string,               noisePol,      // noise policy of v1 - assert compatibility with M(rho,rho) 
                                    std::vector<unsigned int>, tSrcs,         // source times
                                    std::string,               gammas,        // list of gamma matrices to go at sink
                                    std::vector<std::string>,  moms,          // list of momenta 
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
    //envTmpLat(ComplexField,      "ph");
    //envTmp   (ComplexField,      "ph3d"          ,1, gridLD);
    envTmpLat(ComplexField,      "coor");

    auto &dilNoise = envGet(DistillationNoise<FImpl>, par().noisePol);
    int nDL = dilNoise.dilutionSize(DistillationNoise<FImpl>::Index::l);        
    int nDS = dilNoise.dilutionSize(DistillationNoise<FImpl>::Index::s);     

    Grid::Coordinate coor  = gridHD->GlobalDimensions();
    coor[3] = nDL * nDS;
    Grid::GridCartesian * gridDD = Grid::SpaceTimeGrid::makeFourDimGrid(coor, GridDefaultSimd(Nd,vComplex::Nsimd()), GridDefaultMpi());

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
            RhoRhoMesonMFs.try_emplace(RhoRhoGamma+"_p"+kmom, ContractionDistilMesonField<ComplexD,ComplexF>(mfPath, nT, timer1, 0, nT-1, ""));
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

    std::vector<unsigned int> tSrcs = par().tSrcs;
    for(auto tSrc : tSrcs)
    {
        if(tSrc>=nT)
        {
            HADRONS_ERROR(Range, "all tSrcs must be smaller than nT");
        }
    }

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

            // read perambulator
            LOG(Message) << "Starting charm perambulator I/O for (tSrc,tSnk) = (" << tSrc << "," << tSnk << ")" << std::endl;
            envGetTmp(FermionField,    fermionDDtmp_charm);
            startTimer("phi_c I/O");
            tFileName = par().vectorStemC;
            tFileName.append("_DD");
            tFileName.append("_tSm");
            tFileName.append(std::to_string(tSrc));
            tFileName.append("_tLoc");
            tFileName.append(std::to_string(tSnk));
            DistillationVectorsIo::readComponent(fermionDDtmp_charm, tFileName, 1, nDL, nDS, nDT, 0, vm().getTrajectory());
            stopTimer("phi_c I/O");

            // read perambulator
            LOG(Message) << "Starting light perambulator I/O for (tSrc,tSnk) = (" << tSrc << "," << tSnk << ")" << std::endl;
            envGetTmp(FermionField,    fermionDDtmp_light);
            startTimer("phi_l I/O");
            tFileName = par().vectorStemL;
            tFileName.append("_DD");
            tFileName.append("_tSm");
            tFileName.append(std::to_string(tSrc));
            tFileName.append("_tLoc");
            tFileName.append(std::to_string(tSnk));
            DistillationVectorsIo::readComponent(fermionDDtmp_light, tFileName, 1, nDL, nDS, nDT, 0, vm().getTrajectory());
            stopTimer("phi_l I/O");

            unsigned int tdx = tSrci*nMoms*gammas.size();
            for (unsigned int ddx = 0; ddx < momSrcs.size(); ddx++) 
            {
                std::string momSrc = momSrcs[ddx];
                std::string momSnk = momSnks[ddx];
                // momentum phase e^{ipx} for Hw
                //Complex           i(0.0,1.0);
                //std::vector<Real> p;
                //p  = strToVec<Real>(momSnk);
                //ph = Zero();
                //for(unsigned int mu = 0; mu < env().getNd(); mu++)
                //{
                //    LatticeCoordinate(coor, mu);
                //    ph = ph + (p[mu]/env().getDim(mu))*coor;
                //}
                //ph = exp((Real)(2*M_PI)*i*ph);
                //// 3D phase e^{ipx}
                //ExtractSliceLocal(ph3d,ph,0,t,Tdir);  

                ContractionDistilMesonField<ComplexD,ComplexF> &SrcMF = RhoRhoMesonMFs.at(RhoRhoGamma+"_p"+momSrc);
                DistilMesonFieldMatrix<ComplexD> MFmult = SrcMF(tSrc,tSrc,tSrc);

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
                        ExtractSliceLocal(fermion3dtmp1, fermionDDtmp_charm, 0, id1, Tdir);
                        ExtractSliceLocal(fermion3dtmp4, fermionDDtmp_light, 0, id1, Tdir);
                        stopTimer("ExtractSliceLocal");
                        for (int id2=0; id2<nDL*nDS; id2++)
                        {
                            startTimer("ExtractSliceLocal");
                            ExtractSliceLocal(fermion3dtmp2, fermionDDtmp_light, 0, id2, Tdir);
                            stopTimer("ExtractSliceLocal");
                            startTimer("computation");
                            fermion3dtmp3 = gam*fermion3dtmp2;
                            prop3dtmp = outerProductC(fermion3dtmp1, fermion3dtmp3);
                            MesPhi1 += trace(prop3dtmp)*MFmult(id2,id1);
                            stopTimer("computation");

                            startTimer("computation");
                            fermion3dtmp3 = gam*fermion3dtmp2;
                            prop3dtmp = outerProductC(fermion3dtmp4, fermion3dtmp3);
                            MesPhi3 += trace(prop3dtmp)*MFmult(id2,id1);
                            stopTimer("computation");

                            startTimer("ExtractSliceLocal");
                            ExtractSliceLocal(fermion3dtmp2, fermionDDtmp_charm, 0, id2, Tdir);
                            stopTimer("ExtractSliceLocal");
                            startTimer("computation");
                            fermion3dtmp3 = gam*fermion3dtmp2;
                            prop3dtmp1 = outerProductC(fermion3dtmp1, fermion3dtmp3);
                            MesPhi2 += trace(prop3dtmp1)*MFmult(id2,id1);
                            stopTimer("computation");
                        }
                    }

                    startTimer("final contraction cl");
                    sliceSum(MesPhi1, clBuf, Tdir);

                    LOG(Message) << "Updating clResults for (tdx,tSnk) = (" << tdx << "," << tSnk << ")" << std::endl;
                    clResults[tdx].corr[t] = TensorRemove(clBuf[0]);

                    if (t == 0) // only edit metadata on first tSnk for each tSrc
                    {
                        LOG(Message) << "Updating metadata for (tdx,tSnk) = (" << tdx << "," << tSnk << ")" << std::endl;
                        clResults[tdx].gammaSnk = gamStr.str();
                        clResults[tdx].gammaSrc = RhoRhoGamma;
                        clResults[tdx].momSnk   = momSnk;
                        clResults[tdx].momSrc   = momSrc;
                        clResults[tdx].tSrc     = tSrc;
                    }
                    stopTimer("final contraction cl");

                    startTimer("final contraction cc");
                    sliceSum(MesPhi2, ccBuf, Tdir);

                    LOG(Message) << "Updating ccResults for (tdx,tSnk) = (" << tdx << "," << tSnk << ")" << std::endl;
                    ccResults[tdx].corr[t] = TensorRemove(ccBuf[0]);

                    if (t == 0) // only edit metadata on first tSnk for each tSrc
                    {
                        LOG(Message) << "Updating metadata for (tdx,tSnk) = (" << tdx << "," << tSnk << ")" << std::endl;
                        ccResults[tdx].gammaSnk = gamStr.str();
                        ccResults[tdx].gammaSrc = RhoRhoGamma;
                        ccResults[tdx].momSnk   = momSnk;
                        ccResults[tdx].momSrc   = momSrc;
                        ccResults[tdx].tSrc     = tSrc;
                    }
                    stopTimer("final contraction cc");

                    startTimer("final contraction ll");
                    sliceSum(MesPhi3, llBuf, Tdir);

                    LOG(Message) << "Updating llResults for (tdx,tSnk) = (" << tdx << "," << tSnk << ")" << std::endl;
                    llResults[tdx].corr[t] = TensorRemove(llBuf[0]);

                    if (t == 0) // only edit metadata on first tSnk for each tSrc
                    {
                        LOG(Message) << "Updating metadata for (tdx,tSnk) = (" << tdx << "," << tSnk << ")" << std::endl;
                        llResults[tdx].gammaSnk = gamStr.str();
                        llResults[tdx].gammaSrc = RhoRhoGamma;
                        llResults[tdx].momSnk   = momSnk;
                        llResults[tdx].momSrc   = momSrc;
                        llResults[tdx].tSrc     = tSrc;
                    }
                    stopTimer("final contraction ll");

                    tdx++;
                }
            }
        }
    }
    startTimer("results io");
    LOG(Message) << "Writing results to " << par().output << std::endl;
    saveResult(par().output+"_cl", "clMeson", clResults);
    auto &clOut = envGet(HadronsSerializable, getName()+"_cl");
    clOut = clResults;
    saveResult(par().output+"_cc", "ccMeson", ccResults);
    auto &ccOut = envGet(HadronsSerializable, getName()+"_cc");
    ccOut = ccResults;
    saveResult(par().output+"_ll", "llMeson", llResults);
    auto &llOut = envGet(HadronsSerializable, getName()+"_ll");
    llOut = llResults;
    stopTimer("results io");
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MDistil_UnsmearedMeson_hpp_
