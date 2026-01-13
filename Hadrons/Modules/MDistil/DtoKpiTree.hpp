#ifndef Hadrons_MDistil_DtoKpiTree_hpp_
#define Hadrons_MDistil_DtoKpiTree_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/TimerArray.hpp>
#include <Hadrons/Serialization.hpp>
#include <Hadrons/Modules/MDistil/Base.hpp>

BEGIN_HADRONS_NAMESPACE

/*************************************************************************
 *                         DtoKpiTree                                    *
 * Computes the following diagram:                                       * 
 *                                                                       * 
 *                                                                       * 
 *           ______ v1   v3 ______                                       * 
 *          /      g12 x g34      \                                      * 
 *         /        v5   v4        \                                     * 
 *        /          \      \______ M(rho3,rho4)                         * 
 *  M(rho1,rho2)      \                                                  * 
 *        \            \__________                                       * 
 *         \                      \                                      * 
 *          \______________________ M(phi2,rho5)                         * 
 *                                                                       * 
 *                                                                       * 
 *  D(t=tD)           H_W(t)           tKpi                              * 
 *                                                                       * 
 *                                                                       * 
 *************************************************************************/
BEGIN_MODULE_NAMESPACE(MDistil)

typedef std::pair<Gamma::Algebra, Gamma::Algebra> GammaPair;

class DtoKpiTreePar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(DtoKpiTreePar,
                                    std::string,                output,       // file stem for the out file
                                    std::string,                RhoRhoStem,    // file stem for the rho-rho MFs
                                    std::string,                RhoPhiStem,    // file stem for the rho-phi MFs
                                    std::string,                RhoRhoField,   // M(rho,rho) meson field for the Kpi
                                    std::string,                RhoPhiField,   // M(rho,phi) meson field for the Kpi
                                    std::string,                DMesonField,   // M(rho,rho) meson field for the D
                                    std::string,                vectorStemC,   // charm
                                    std::string,                vectorStemL,   // SU(3) light
                                    std::string,                noisePol,      // noise policy of v1 - assert compatibility with M(rho,rho) 
                                    std::vector<unsigned int>,  tDs,           // times of D meson
                                    std::vector<unsigned int>,  tKpis,         // times of Kpi mesons
                                    std::string,                gammas,        // list of space-separated pairs of gamma matrices (g12 g34)
                                    std::string,                mom,           // momentum injected into Hw
                                                                                );
};

template <typename FImpl>
class TDtoKpiTree: public Module<DtoKpiTreePar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    class Result: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(Result,
                                        std::string,          gammaD,
                                        std::string,          gammaKpi_rhorho,
                                        std::string,          gammaKpi_rhophi,
                                        std::string,          gammaHw,
                                        unsigned int,         tD,
                                        unsigned int,         tKpi,
                                        std::vector<Complex>, corr);
    };
    // constructor
    TDtoKpiTree(const std::string name);
    // destructor
    virtual ~TDtoKpiTree(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(DtoKpiTree, TDtoKpiTree<FIMPL>, MDistil);

/******************************************************************************
 *                 TDtoKpiTree implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TDtoKpiTree<FImpl>::TDtoKpiTree(const std::string name)
: Module<DtoKpiTreePar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TDtoKpiTree<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().noisePol};
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TDtoKpiTree<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TDtoKpiTree<FImpl>::setup(void)
{
    GridCartesian * gridHD = envGetGrid(FermionField);
    GridCartesian * gridLD = envGetSliceGrid(FermionField,gridHD->Nd() -1);
    
    envTmp   (FermionField,    "fermion3dtmp1" ,1, gridLD);
    envTmp   (FermionField,    "fermion3dtmp2" ,1, gridLD);
    envTmp   (FermionField,    "fermion3dtmp3" ,1, gridLD);
    envTmp   (PropagatorField, "prop3dtmp"     ,1, gridLD);
    envTmp   (ComplexField,    "MKpiPhi"       ,1, gridLD);
    envTmp   (ComplexField,    "MDPhi"         ,1, gridLD);
    envTmpLat(ComplexField,    "ph");
    envTmp   (ComplexField,    "ph3d"          ,1, gridLD);
    envTmpLat(ComplexField,    "coor");

    auto &dilNoise = envGet(DistillationNoise<FImpl>, par().noisePol);
    int nDL = dilNoise.dilutionSize(DistillationNoise<FImpl>::Index::l);        
    int nDS = dilNoise.dilutionSize(DistillationNoise<FImpl>::Index::s);     

    Grid::Coordinate coor  = gridHD->GlobalDimensions();
    coor[3] = nDL * nDS;
    Grid::GridCartesian * gridDD = Grid::SpaceTimeGrid::makeFourDimGrid(coor, GridDefaultSimd(Nd,vComplex::Nsimd()), GridDefaultMpi());

    envTmp(FermionField,    "fermionDDtmp_light" ,1, gridDD);
    envTmp(FermionField,    "fermionDDtmp_charm" ,1, gridDD);

    envCreate(HadronsSerializable, getName(), 1, 0);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TDtoKpiTree<FImpl>::execute(void)
{
    // general grid setup
    GridCartesian * gridHD = envGetGrid(FermionField);
    GridCartesian * gridLD = envGetSliceGrid(FermionField,gridHD->Nd() -1);
    const int Ntlocal{gridHD->LocalDimensions()[Tdir]};
    const int Ntfirst{gridHD->LocalStarts()[Tdir]};
    int nT=env().getDim(Tdir);

    // read input D-meson field
    std::string mfPath = par().RhoRhoStem + "rho-rho." + std::to_string(vm().getTrajectory()) + "/" + par().DMesonField;   
    LOG(Message) << "reading " << mfPath << std::endl;
    TimerArray timer;
    // TODO: consider making a loader to directly get every e.g. 4th time slice instead of sequentially 
    ContractionDistilMesonField<ComplexD,ComplexF> DMesonMF(mfPath, nT, timer, 0, nT-1, "");
    std::string DGamma = par().DMesonField.substr(0,par().DMesonField.find('.'));

    // read input Kpi rhorho field
    mfPath = par().RhoRhoStem + "rho-rho." + std::to_string(vm().getTrajectory()) + "/" + par().RhoRhoField;   
    LOG(Message) << "reading " << mfPath << std::endl;
    TimerArray timer1;
    // TODO: consider making a loader to directly get every e.g. 4th time slice instead of sequentially 
    ContractionDistilMesonField<ComplexD,ComplexF> RhoRhoMF(mfPath, nT, timer1, 0, nT-1, "");
    std::string RhoRhoGamma = par().RhoRhoField.substr(0,par().RhoRhoField.find('.'));

    // read input Kpi rhophi field
    mfPath = par().RhoPhiStem + "rho-phi." + std::to_string(vm().getTrajectory()) + "/" + par().RhoPhiField;   
    LOG(Message) << "reading " << mfPath << std::endl;
    TimerArray timer2;
    // TODO: consider making a loader to directly get every e.g. 4th time slice instead of sequentially 
    ContractionDistilMesonField<ComplexD,ComplexF> RhoPhiMF(mfPath, nT, timer2);
    std::string RhoPhiGamma = par().RhoPhiField.substr(0,par().RhoPhiField.find('.'));

    // noise class -- assert they are identical and an "exact distillation" policy
    auto &dilNoise = envGet(DistillationNoise<FImpl>, par().noisePol);
    int nNoise = dilNoise.size(); 
    if(nNoise>1)
    {
        HADRONS_ERROR(Implementation, "DtoKpiTree only implemented for exact distillation");
    }
    int nDL = dilNoise.dilutionSize(DistillationNoise<FImpl>::Index::l);        
    int nDS = dilNoise.dilutionSize(DistillationNoise<FImpl>::Index::s);        
    int nDT = dilNoise.dilutionSize(DistillationNoise<FImpl>::Index::t);        

    std::vector<unsigned int> tKpis = par().tKpis; 
    for(auto tKpi : tKpis)
    {
        if(tKpi>=nT)
        {
            HADRONS_ERROR(Range, "all tKpis must be smaller than nT");
        }
    }

    std::vector<unsigned int> tDs = par().tDs;

    std::vector<GammaPair> gammas = strToVec<GammaPair>(par().gammas);

    std::vector<Result> results;
    results.resize(gammas.size()*tDs.size()*tKpis.size());
    
    // Temporary objects
    envGetTmp(FermionField,    fermion3dtmp1);
    envGetTmp(FermionField,    fermion3dtmp2);
    envGetTmp(FermionField,    fermion3dtmp3);
    envGetTmp(PropagatorField, prop3dtmp);
    envGetTmp(ComplexField,    MKpiPhi);
    envGetTmp(ComplexField,    MDPhi);

    // momentum phase e^{ipx} for Hw
    Complex           i(0.0,1.0);
    std::vector<Real> p;
    p  = strToVec<Real>(par().mom);
    envGetTmp(ComplexField, coor);
    envGetTmp(ComplexField, ph);
    envGetTmp(ComplexField, ph3d);
    ph = Zero();
    for(unsigned int mu = 0; mu < env().getNd(); mu++)
    {
        LatticeCoordinate(coor, mu);
        ph = ph + (p[mu]/env().getDim(mu))*coor;
    }
    ph = exp((Real)(2*M_PI)*i*ph);
        
    envGetTmp(FermionField,    fermionDDtmp_light);
    envGetTmp(FermionField,    fermionDDtmp_charm);

    unsigned int rdx = 0;
    for (unsigned int tDi = 0; tDi < tDs.size(); tDi++)
    {
        unsigned int tD = tDs[tDi];

        if(tD>=nT)
        {
            HADRONS_ERROR(Range, "tD must be smaller than nT");
        }
        // determine timeslices tH which are between tD and tKpi (shorter distance)
        for(auto tKpi : tKpis)
        {
            std::vector<unsigned int> tHs;
            int tDMinusTKpi = (tD - tKpi + nT) % nT;
            int tKpiMinusTD = (tKpi - tD + nT) % nT;
            if(tDMinusTKpi < tKpiMinusTD)
            {
                for(int iTH = 1; iTH < tDMinusTKpi; iTH++)
                {
                    int tH_tmp = (tKpi + iTH + nT) % nT;
                    tHs.push_back(tH_tmp);      
                }     
            }
            else
            {
                for(int iTH = 1; iTH < tKpiMinusTD; iTH++)
                {
                    int tH_tmp = (tD + iTH + nT) % nT;
                    tHs.push_back(tH_tmp);      
                }     
            }

            for (unsigned int i = 0; i < gammas.size(); i++)
            {
                unsigned int ridx = rdx*gammas.size() + i;
                results[ridx].corr.resize(tHs.size());
            }

            LOG(Message) << "WARNING: Assuming ordering s + ns*(l + nl*t) in DilutedNoise.hpp. This code will break when this changes!" << std::endl;
            int tH;
            std::vector<TComplex>  buf;
            std::string tFileName;

            for (int t = 0; t < Ntlocal; t++)
            {
                tH = t + Ntfirst;
                auto it = std::find(tHs.begin(), tHs.end(), tH);
                if (it == tHs.end())
                {
                    LOG(Message) << "Only computing three-point function between tD = " << tD << " and tKpi = " << tKpi << ", skipping tH = " << tH << std::endl;
                    continue;
                }   
                LOG(Message) << "Starting perambulator I/O for (tD,tKpi,tH) = (" << tD << "," << tKpi << "," << tH << ")" << std::endl;
                int tHi = std::distance(tHs.begin(), it);

                // read perambulators
                envGetTmp(FermionField,    fermionDDtmp_charm);
                startTimer("phi_c I/O");
                tFileName = par().vectorStemC;
                tFileName.append("_DD");
                tFileName.append("_tSm");
                tFileName.append(std::to_string(tD));
                tFileName.append("_tLoc");
                tFileName.append(std::to_string(tH));
                DistillationVectorsIo::readComponent(fermionDDtmp_charm, tFileName, 1, nDL, nDS, nDT, 0, vm().getTrajectory());
                stopTimer("phi_c I/O");

                envGetTmp(FermionField,    fermionDDtmp_light);
                startTimer("phi_l I/O");
                tFileName = par().vectorStemL;
                tFileName.append("_DD");
                tFileName.append("_tSm");
                tFileName.append(std::to_string(tKpi));
                tFileName.append("_tLoc");
                tFileName.append(std::to_string(tH));
                DistillationVectorsIo::readComponent(fermionDDtmp_light, tFileName, 1, nDL, nDS, nDT, 0, vm().getTrajectory());
                stopTimer("phi_l I/O");

                // 3D phase e^{ipx}
                ExtractSliceLocal(ph3d,ph,0,t,Tdir);  

                unsigned int sdx = 0;
                for (unsigned int sdx = 0; sdx < gammas.size(); sdx++)
                {
                    unsigned int tdx = rdx*gammas.size() + sdx;

                    Gamma::Algebra gam12 = gammas[sdx].first, gam34 = gammas[sdx].second;
                    Gamma g12(gam12), g34(gam34);

                    //   contract 2xphi_l with Kpi(rho,rho)
                    // & contract phi_l, phi_c, DMesonMF
                    // TODO: this needs to be smarter for colour-suppressed diagram variant
                    for (int id1=0; id1<nDL*nDS; id1++)
                    {
                        ExtractSliceLocal(fermion3dtmp1, fermionDDtmp_light, 0, id1, Tdir);
                        startTimer("computation contractPhis");
                        for (int id2=0; id2<nDL*nDS; id2++)
                        {
                            ExtractSliceLocal(fermion3dtmp2, fermionDDtmp_light, 0, id2, Tdir);
                            fermion3dtmp3 = g12*fermion3dtmp2;
                            fermion3dtmp2 = fermion3dtmp3*RhoRhoMF(tKpi,tKpi,tKpi)(id1,id2);
                            prop3dtmp = outerProduct(fermion3dtmp1, fermion3dtmp2);
                            // this object is sum_{spin,colour,d1,d2} (DMeson[d1,d2] * vector1[d1] * gamma12 * vector2[d2]) on timeslice tH
                            MKpiPhi += trace(prop3dtmp);

                            ExtractSliceLocal(fermion3dtmp2, fermionDDtmp_charm, 0, id2, Tdir);
                            fermion3dtmp3 = g34*fermion3dtmp1;
                            fermion3dtmp1 = fermion3dtmp3*RhoPhiMF(tKpi,tKpi,tD)(id1,id2)*DMesonMF(tD,tD,tD)(id1,id2);
                            prop3dtmp = outerProduct(fermion3dtmp2, fermion3dtmp1);
                            MDPhi += trace(prop3dtmp);
                        }
                        stopTimer("computation contractPhis");
                    }
                    startTimer("final contraction");
                    MKpiPhi = MDPhi*MKpiPhi*ph3d;
                    sliceSum(MKpiPhi, buf, Tdir);

                    LOG(Message) << "Updating results for (rdx,sdx,tdx,tH) = (" << rdx << "," << sdx << "," << tdx << "," << tH << ")" << std::endl;
                    results[tdx].corr[tHi] = TensorRemove(buf[0]);

                    stopTimer("final contraction");

                    if (tHi == 0) // only edit metadata on first tH for each (tD,tKpi)
                    {
                        LOG(Message) << "Updating metadata for (rdx,sdx,tdx,tH) = (" << rdx << "," << sdx << "," << tdx << "," << tH << ")" << std::endl;
                        std::stringstream gHw;
                        gHw << "(" << gam12 << " " << gam34 << ")";
                        results[tdx].gammaHw         = gHw.str();
                        results[tdx].gammaD          = DGamma;
                        results[tdx].gammaKpi_rhorho = RhoRhoGamma;
                        results[tdx].gammaKpi_rhophi = RhoPhiGamma;
                        results[tdx].tD              = tD;
                        results[tdx].tKpi            = tKpi;
                    }
                }
            }
            rdx++;
        }
    }
    LOG(Message) << "Writing results to " << par().output << std::endl;
    saveResult(par().output, "DtoKpiTree", results);
    auto &out = envGet(HadronsSerializable, getName());
    out = results;
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MDistil_DtoKpiTree_hpp_
