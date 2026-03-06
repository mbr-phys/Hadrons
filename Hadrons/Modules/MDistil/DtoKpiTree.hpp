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
 *  Sresults : colour singlet bilinears in Hamiltonian                   *
 *             -> two colour traces                                      *
 *  Rresults : colour rearranged bilinears in Hamiltonian                *
 *             -> one colour trace                                       *
 *                                                                       * 
 *************************************************************************/
BEGIN_MODULE_NAMESPACE(MDistil)

typedef std::pair<Gamma::Algebra, Gamma::Algebra> GammaPair;

class DtoKpiTreePar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(DtoKpiTreePar,
                                    std::string,               output,        // file stem for the out file
                                    std::string,               RhoRhoStem,    // file stem for the rho-rho MFs
                                    std::string,               RhoPhiStem,    // file stem for the rho-phi MFs
                                    std::string,               RhoRhoField,   // M(rho,rho) meson field for the Kpi
                                    std::string,               RhoPhiField,   // M(rho,phi) meson field for the Kpi
                                    std::string,               DMesonField,   // M(rho,rho) meson field for the D
                                    std::string,               vectorStemC,   // charm
                                    std::string,               vectorStemL,   // SU(3) light
                                    std::string,               noisePol,      // noise policy of v1 - assert compatibility with M(rho,rho) 
                                    std::vector<unsigned int>, tDs,           // times of D meson
                                    std::vector<unsigned int>, tKpis,         // times of Kpi mesons
                                    std::string,               gammas,        // list of space-separated pairs of gamma matrices (g12 g34)
                                    std::string,               momHw,         // momentum injected into Hw
                                    std::vector<std::string>,  momsD,         // list of momenta of incoming D meson
                                    std::vector<std::string>,  momsKpi,       // list of possible momenta of final states (if empty, will assume all up to P^2=4
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
                                        std::string,          momD,
                                        std::string,          momKpi_rhorho,
                                        std::string,          momKpi_rhophi,
                                        std::string,          momHw,
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
    std::vector<std::string> out = {getName()+"_tree",getName()+"_colour"};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TDtoKpiTree<FImpl>::setup(void)
{
    GridCartesian * gridHD = envGetGrid(FermionField);
    GridCartesian * gridLD = envGetSliceGrid(FermionField,gridHD->Nd() -1);
    
    envTmp   (FermionField,      "fermion3dtmp1" ,1, gridLD);
    envTmp   (FermionField,      "fermion3dtmp2" ,1, gridLD);
    envTmp   (FermionField,      "fermion3dtmp3" ,1, gridLD);
    envTmp   (PropagatorField,   "prop3dtmp"     ,1, gridLD);
    envTmp   (PropagatorField,   "prop3dtmp1"     ,1, gridLD);
    envTmp   (ColourMatrixField, "MKpiPhi"       ,1, gridLD);
    envTmp   (ColourMatrixField, "MDPhi"         ,1, gridLD);
    envTmp   (ComplexField,      "MColour"       ,1, gridLD);
    envTmpLat(ComplexField,      "ph");
    envTmp   (ComplexField,      "ph3d"          ,1, gridLD);
    envTmpLat(ComplexField,      "coor");

    auto &dilNoise = envGet(DistillationNoise<FImpl>, par().noisePol);
    int nDL = dilNoise.dilutionSize(DistillationNoise<FImpl>::Index::l);        
    int nDS = dilNoise.dilutionSize(DistillationNoise<FImpl>::Index::s);     

    Grid::Coordinate coor  = gridHD->GlobalDimensions();
    coor[3] = nDL * nDS;
    Grid::GridCartesian * gridDD = Grid::SpaceTimeGrid::makeFourDimGrid(coor, GridDefaultSimd(Nd,vComplex::Nsimd()), GridDefaultMpi());

    envTmp(FermionField,    "fermionDDtmp_light" ,1, gridDD);
    envTmp(FermionField,    "fermionDDtmp_charm" ,1, gridDD);

    envCreate(HadronsSerializable, getName()+"_singlet", 1, 0);
    envCreate(HadronsSerializable, getName()+"_rearranged", 1, 0);
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

    std::vector<std::string> Dmoms, Kpilist;
    std::map<std::string, std::vector<std::vector<std::string>>> Kpimoms;
    int nMoms = 0;

    if (par().momsKpi.empty()) 
    {
        LOG(Message) << "You have not specified any K-pi momenta, so all possible combinations up to P^2 = 4 will be listed." << std::endl;
    }
    else
    {
        LOG(Message) << "Using " << par().momsKpi.size() << " possible K-pi momenta provided." << std::endl;
    }
    LOG(Message) << "WARNING: Assuming ordering s + ns*(l + nl*t) in DilutedNoise.hpp. This code will break when this changes!" << std::endl;

    for (auto dmom : par().momsD)
    {
        std::vector<int> dmomI = strToVec<int>(dmom);
        std::string dstr = std::to_string(dmomI[0]) + "_" + std::to_string(dmomI[1]) + "_" + std::to_string(dmomI[2]);
        Dmoms.push_back(dstr);
        // if no momsKpi given, loop over all possible combinations up to P^2 = 4
        std::vector<std::vector<std::string>> Kpims;
        if (par().momsKpi.empty())
        {
            for (int i = -2; i <= 2; i++)
            {
                for (int j = -2; j <= 2; j++)
                {
                    for (int k = -2; k <= 2; k++)
                    {
                        int Kpi1P2 = i*i + j*j + k*k;
                        if (Kpi1P2 <= 4)
                        {
                            std::string Kpimom1 = std::to_string(i) + "_" + std::to_string(j) + "_" + std::to_string(k);
                            int o(dmomI[0]-i), p(dmomI[1]-j), q(dmomI[2]-k);
                            int Kpi2P2 = o*o + p*p + q*q;
                            if ((std::abs(o) <= 2) && (std::abs(p) <= 2) && (std::abs(q) <= 2) && (Kpi2P2 <= 4))
                            {
                                std::string Kpimom2 = std::to_string(o) + "_" + std::to_string(p) + "_" + std::to_string(q);
                                std::vector<std::string> Kpis = {Kpimom1, Kpimom2};
                                Kpims.push_back(Kpis);
                                if (std::find(Kpilist.begin(), Kpilist.end(), Kpimom1) == Kpilist.end())
                                {
                                    Kpilist.push_back(Kpimom1);
                                }
                                if (std::find(Kpilist.begin(), Kpilist.end(), Kpimom2) == Kpilist.end())
                                {
                                    Kpilist.push_back(Kpimom2);
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
            for (auto kmom : par().momsKpi) 
            {
                std::vector<int> kmomI = strToVec<int>(kmom);
                std::string Kpimom1 = std::to_string(kmomI[0]) + "_" + std::to_string(kmomI[1]) + "_" + std::to_string(kmomI[2]);
                int Kpi1P2 = kmomI[0]*kmomI[0] + kmomI[1]*kmomI[1] + kmomI[2]*kmomI[2];
                int o(dmomI[0]-kmomI[0]), p(dmomI[1]-kmomI[1]), q(dmomI[2]-kmomI[2]);
                int Kpi2P2 = o*o + p*p + q*q;
                if ((std::abs(o) <= 2) && (std::abs(p) <= 2) && (std::abs(q) <= 2) && (Kpi2P2 <= 4))
                {
                    std::string Kpimom2 = std::to_string(o) + "_" + std::to_string(p) + "_" + std::to_string(q);
                    std::vector<std::string> Kpis = {Kpimom1, Kpimom2};
                    Kpims.push_back(Kpis);
                    if (std::find(Kpilist.begin(), Kpilist.end(), Kpimom1) == Kpilist.end())
                    {
                        Kpilist.push_back(Kpimom1);
                    }
                    if (std::find(Kpilist.begin(), Kpilist.end(), Kpimom2) == Kpilist.end())
                    {
                        Kpilist.push_back(Kpimom2);
                    }
                    nMoms++;
                }
            }
        }
        Kpimoms.try_emplace(dstr, Kpims);
        LOG(Message) << "Will calculate " << Kpims.size() << " K-pi final state combinations for D meson momentum = " + dmom << std::endl;
    }

    // read input D-meson field
    std::string DGamma = par().DMesonField;
    std::map<std::string, ContractionDistilMesonField<ComplexD,ComplexF>> RhoRhoMesonMFs;
    startTimer("MesonField IO");
    for (auto dmom : Dmoms)
    {
        std::string mfPath = par().RhoRhoStem + "rho-rho." + std::to_string(vm().getTrajectory()) + "/" + DGamma + "_p" + dmom + ".h5";   
        LOG(Message) << "reading " << mfPath << std::endl;
        TimerArray timer;
        // TODO: consider making a loader to directly get every e.g. 4th time slice instead of sequentially 
        RhoRhoMesonMFs.try_emplace(DGamma+"_p"+dmom, ContractionDistilMesonField<ComplexD,ComplexF>(mfPath, nT, timer, 0, nT-1, ""));
    }

    // read input Kpi fields
    std::string RhoRhoGamma = par().RhoRhoField;
    std::string RhoPhiGamma = par().RhoPhiField;
    std::map<std::string, ContractionDistilMesonField<ComplexD,ComplexF>> RhoPhiMesonMFs;
    for (auto kmom : Kpilist)
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
    for(auto tD : tDs)
    {
        if(tD>=nT)
        {
            HADRONS_ERROR(Range, "all tDs must be smaller than nT");
        }
    }

    std::vector<GammaPair> gammas = strToVec<GammaPair>(par().gammas);

    std::vector<Result> Sresults, Rresults;
    int resultSize = gammas.size()*tDs.size()*tKpis.size()*nMoms;
    Sresults.resize(resultSize);
    Rresults.resize(resultSize);
    LOG(Message) << "Results objects have gammas (" << gammas.size() << ") * tDs (" << tDs.size() 
                 << ") * tKpis (" << tKpis.size() << ") * nMoms (" << nMoms << ") = " << resultSize << " size" << std::endl;
    int counter = 0;
    for (unsigned int tDi = 0; tDi < tDs.size(); tDi++)
    {
        unsigned int tD = tDs[tDi];
        for (unsigned int tKpii = 0; tKpii < tKpis.size(); tKpii++)
        {
            unsigned int tKpi = tKpis[tKpii];
            unsigned int size = (tKpi-tD <= nT/2) ? ((tKpi - tD + nT)%nT - 1) : ((tD - tKpi + nT)%nT - 1);
            for (unsigned int i = 0; i < nMoms*gammas.size(); i++)
            {
                unsigned int ridx = counter*nMoms*gammas.size() + i;
                //LOG(Message) << "Resizing (counter,ridx) = (" << counter << "," << ridx << ")" << std::endl;
                Sresults[ridx].corr.resize(size);
                Rresults[ridx].corr.resize(size);
            }
            counter++;
        }
    }
    
    // Temporary objects
    envGetTmp(FermionField,       fermion3dtmp1);
    envGetTmp(FermionField,       fermion3dtmp2);
    envGetTmp(FermionField,       fermion3dtmp3);
    envGetTmp(PropagatorField,    prop3dtmp);
    envGetTmp(PropagatorField,    prop3dtmp1);
    envGetTmp(ColourMatrixField,  MKpiPhi);
    envGetTmp(ColourMatrixField,  MDPhi);
    envGetTmp(ComplexField,       MColour);

    // momentum phase e^{ipx} for Hw
    Complex           i(0.0,1.0);
    std::vector<Real> p;
    p  = strToVec<Real>(par().momHw);
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

    for (unsigned int tDi = 0; tDi < tDs.size(); tDi++)
    {
        unsigned int tD = tDs[tDi];

        if(tD>=nT)
        {
            HADRONS_ERROR(Range, "tD must be smaller than nT");
        }

        int tH;
        std::string tFileName;
        for (int t = 0; t < Ntlocal; t++)
        {
            tH = t + Ntfirst;
            if (tH == tD)
            {
                LOG(Message) << "Not including contact terms, skipping tD = " << tD << " and tH = " << tH << std::endl;
                continue;
            }

            // 3D phase e^{ipx}
            ExtractSliceLocal(ph3d,ph,0,t,Tdir);  

            // read perambulator
            LOG(Message) << "Starting charm perambulator I/O for (tD,tH) = (" << tD << "," << tH << ")" << std::endl;
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

            for(unsigned int tKpii = 0; tKpii < tKpis.size(); tKpii++)
            {
                unsigned int tKpi = tKpis[tKpii];
                if (tH == tKpi)
                {
                    LOG(Message) << "Not including contact terms, skipping tKpi = " << tKpi << " and tH = " << tH << std::endl;
                    continue;
                }

                std::vector<unsigned int> tHs;
                int tDMinusTKpi = (tD - tKpi + nT) % nT;
                int tKpiMinusTD = (tKpi - tD + nT) % nT;
                int tHMinusTKpi = (tH - tKpi + nT) % nT;
                int tKpiMinusTH = (tKpi - tH + nT) % nT;
                int tHMinusTD   = (tH - tD + nT) % nT;
                int tDMinusTH   = (tD - tH + nT) % nT;
                //LOG(Message) << "(tD,tKpi,tH) = (" << tD << "," << tKpi << "," << tH << ")" << std::endl;
                //LOG(Message) << "    tDMinusTKpi = " << tDMinusTKpi << std::endl;
                //LOG(Message) << "    tKpiMinusTD = " << tKpiMinusTD << std::endl;
                //LOG(Message) << "    tHMinusTKpi = " << tHMinusTKpi << std::endl;
                //LOG(Message) << "    tKpiMinusTH = " << tKpiMinusTH << std::endl;
                //LOG(Message) << "      tHMinusTD = " << tHMinusTD << std::endl;
                //LOG(Message) << "      tDMinusTH = " << tDMinusTH << std::endl;
                int tHi;
                // double check indexing here
                if ((tDMinusTKpi < tKpiMinusTD) && (tDMinusTH < tDMinusTKpi)) 
                {
                    tHi = tHMinusTKpi-1; 
                    LOG(Message) << "--> computing backwards signal" << std::endl;
                }
                else if ((tKpiMinusTD <= tDMinusTKpi) && (tHMinusTD < tKpiMinusTD))
                {
                    tHi = tHMinusTD-1;
                    LOG(Message) << "--> computing forwards signal" << std::endl;
                }
                else
                {
                    LOG(Message) << "Only computing three-point functions between tD = " << tD << " and tKpi = " << tKpi << ", skipping tH = " << tH << std::endl;
                    continue;
                }

                std::vector<TComplex> Sbuf, Rbuf;

                // read perambulator
                LOG(Message) << "Starting light perambulator I/O for (tKpi,tH) = (" << tKpi << "," << tH << ")" << std::endl;
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

                unsigned int rdx = tDi*tKpis.size() + tKpii;
                unsigned int tdx = rdx*nMoms*gammas.size();
                for (unsigned int ddx = 0; ddx < Dmoms.size(); ddx++) 
                {
                    std::string dmom = Dmoms[ddx];
                    std::vector<std::vector<std::string>> Kmoms = Kpimoms.at(dmom);

                    ContractionDistilMesonField<ComplexD,ComplexF> &DMesonMF = RhoRhoMesonMFs.at(DGamma+"_p"+dmom);

                    for (unsigned int kdx = 0; kdx < Kmoms.size(); kdx++)
                    {
                        std::string Kmom1 = Kmoms[kdx][0], Kmom2 = Kmoms[kdx][1];

                        ContractionDistilMesonField<ComplexD,ComplexF> &RhoRhoMF = RhoRhoMesonMFs.at(RhoRhoGamma+"_p"+Kmom1);

                        startTimer("MesonField IO");
                        std::string mfPath = par().RhoPhiStem + "rho-phi." + std::to_string(vm().getTrajectory()) + "/" + RhoPhiGamma + "_p" + Kmom2 + ".h5";   
                        LOG(Message) << "reading " << mfPath << std::endl;
                        TimerArray timer2;
                        std::vector<std::vector<int>> tDtKpi = {{(int)tKpi,(int)tD}}; // only read one entry from h5 file
                        ContractionDistilMesonField<ComplexD,ComplexF> RhoPhiMF(mfPath, nT, timer2, tDtKpi);
                        stopTimer("MesonField IO");

                        // pre-contract distillation index id3 between two MFs
                        startTimer("MF mult");
                        DistilMesonFieldMatrix<ComplexD> MFmult;
                        A2AContraction::mul(MFmult, RhoPhiMF(tKpi,tKpi,tD), DMesonMF(tD,tD,tD));
                        stopTimer("MF mult");
                        for (unsigned int sdx = 0; sdx < gammas.size(); sdx++)
                        {
                            Gamma::Algebra gam12 = gammas[sdx].first, gam34 = gammas[sdx].second;
                            Gamma g12(gam12), g34(gam34);

                            MKpiPhi = Zero();
                            MDPhi   = Zero();
                            MColour = Zero();
                            //   contract 2xphi_l with Kpi(rho,rho)
                            // & contract phi_l, phi_c, DMesonMF
                            for (int id1=0; id1<nDL*nDS; id1++)
                            {
                                startTimer("ExtractSliceLocal");
                                ExtractSliceLocal(fermion3dtmp1, fermionDDtmp_light, 0, id1, Tdir);
                                stopTimer("ExtractSliceLocal");
                                for (int id2=0; id2<nDL*nDS; id2++)
                                {
                                    startTimer("ExtractSliceLocal");
                                    ExtractSliceLocal(fermion3dtmp2, fermionDDtmp_light, 0, id2, Tdir);
                                    stopTimer("ExtractSliceLocal");
                                    startTimer("computation contractPhis Tree");
                                    fermion3dtmp3 = g34*fermion3dtmp2;
                                    fermion3dtmp2 = fermion3dtmp3*RhoRhoMF(tKpi,tKpi,tKpi)(id2,id1);
                                    prop3dtmp = outerProduct(fermion3dtmp1, fermion3dtmp2);
                                    // sum_{spin,d1,d2} (vector4[d1] * gamma34 * vector3[d2] * PMF[d2,d1]) 
                                    MKpiPhi += traceSpin(prop3dtmp);
                                    stopTimer("computation contractPhis Tree");

                                    startTimer("ExtractSliceLocal");
                                    ExtractSliceLocal(fermion3dtmp2, fermionDDtmp_charm, 0, id2, Tdir);
                                    stopTimer("ExtractSliceLocal");
                                    startTimer("computation contractPhis Tree");
                                    fermion3dtmp3 = g12*fermion3dtmp1;
                                    prop3dtmp = outerProduct(fermion3dtmp2, fermion3dtmp3);
                                    // sum_{spin,d1,d2,d3} (DMeson[d3,d2] * vector1[d2] * gamma12 * vector2[d1] * PMF[d1,d3]) 
                                    MDPhi += traceSpin(prop3dtmp*MFmult(id1,id2));
                                    stopTimer("computation contractPhis Tree");
                                }
                            }
                            startTimer("final contraction singlet");
                            MColour = traceColour(MDPhi)*ph3d*traceColour(MKpiPhi);
                            sliceSum(MColour, Sbuf, Tdir);

                            LOG(Message) << "Updating Sresults for (tdx,tH) = (" << tdx << "," << tH << ")" << std::endl;
                            Sresults[tdx].corr[tHi] = TensorRemove(Sbuf[0]);

                            if (tHi == 0) // only edit metadata on first tH for each (tD,tKpi)
                            {
                                LOG(Message) << "Updating metadata for (tdx,tH) = (" << tdx << "," << tH << ")" << std::endl;
                                std::stringstream gHw;
                                gHw << "(" << gam12 << " " << gam34 << ")";
                                Sresults[tdx].gammaHw         = gHw.str();
                                Sresults[tdx].gammaD          = DGamma;
                                Sresults[tdx].gammaKpi_rhorho = RhoRhoGamma;
                                Sresults[tdx].gammaKpi_rhophi = RhoPhiGamma;
                                Sresults[tdx].momD            = dmom;
                                Sresults[tdx].momKpi_rhorho   = Kmom1;
                                Sresults[tdx].momKpi_rhophi   = Kmom2;
                                Sresults[tdx].momHw           = par().momHw;
                                Sresults[tdx].tD              = tD;
                                Sresults[tdx].tKpi            = tKpi;
                            }
                            stopTimer("final contraction singlet");

                            startTimer("final contraction rearranged");
                            MColour = traceColour(MDPhi*ph3d*MKpiPhi);
                            sliceSum(MColour, Rbuf, Tdir);

                            LOG(Message) << "Updating Rresults for (tdx,tH) = (" << tdx << "," << tH << ")" << std::endl;
                            Rresults[tdx].corr[tHi] = TensorRemove(Rbuf[0]);

                            if (tHi == 0) // only edit metadata on first tH for each (tD,tKpi)
                            {
                                LOG(Message) << "Updating metadata for (tdx,tH) = (" << tdx << "," << tH << ")" << std::endl;
                                std::stringstream gHw;
                                gHw << "(" << gam12 << " " << gam34 << ")";
                                Rresults[tdx].gammaHw         = gHw.str();
                                Rresults[tdx].gammaD          = DGamma;
                                Rresults[tdx].gammaKpi_rhorho = RhoRhoGamma;
                                Rresults[tdx].gammaKpi_rhophi = RhoPhiGamma;
                                Rresults[tdx].momD            = dmom;
                                Rresults[tdx].momKpi_rhorho   = Kmom1;
                                Rresults[tdx].momKpi_rhophi   = Kmom2;
                                Rresults[tdx].momHw           = par().momHw;
                                Rresults[tdx].tD              = tD;
                                Rresults[tdx].tKpi            = tKpi;
                            }
                            stopTimer("final contraction rearranged");

                            tdx++;
                        }
                    }
                }
            }
        }
    }
    startTimer("results io");
    LOG(Message) << "Writing results to " << par().output << std::endl;
    saveResult(par().output+"_singlet", "DtoKpiTreeSinglet", Sresults);
    auto &Sout = envGet(HadronsSerializable, getName()+"_singlet");
    Sout = Sresults;
    saveResult(par().output+"_rearranged", "DtoKpiTreeRearranged", Rresults);
    auto &Rout = envGet(HadronsSerializable, getName()+"_rearranged");
    Rout = Rresults;
    stopTimer("results io");
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MDistil_DtoKpiTree_hpp_
