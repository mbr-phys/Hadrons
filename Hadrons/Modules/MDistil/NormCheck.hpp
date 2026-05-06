#ifndef Hadrons_MDistil_NormCheck_hpp_
#define Hadrons_MDistil_NormCheck_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/TimerArray.hpp>
#include <Hadrons/Serialization.hpp>
#include <Hadrons/Modules/MDistil/Base.hpp>

BEGIN_HADRONS_NAMESPACE

/************************************************************
 *                   NormCheck                         *
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

class NormCheckPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(NormCheckPar,
                                    std::string,               vectorStemC,   // charm
                                    std::string,               vectorStemL,   // light
                                    std::string,               noisePol,      // noise policy of v1 - assert compatibility with M(rho,rho) 
                                    std::vector<unsigned int>, tSrcs,         // source times
                                   );
};

template <typename FImpl>
class TNormCheck: public Module<NormCheckPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    // constructor
    TNormCheck(const std::string name);
    // destructor
    virtual ~TNormCheck(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(NormCheck, TNormCheck<FIMPL>, MDistil);

/******************************************************************************
 *                 TNormCheck implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TNormCheck<FImpl>::TNormCheck(const std::string name)
: Module<NormCheckPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TNormCheck<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().noisePol};
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TNormCheck<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TNormCheck<FImpl>::setup(void)
{
    GridCartesian * gridHD = envGetGrid(FermionField);
    GridCartesian * gridLD = envGetSliceGrid(FermionField,gridHD->Nd() -1);
    
    envTmp   (FermionField,      "fermion3dtmp1" ,1, gridLD);

    auto &dilNoise = envGet(DistillationNoise<FImpl>, par().noisePol);
    int nDL = dilNoise.dilutionSize(DistillationNoise<FImpl>::Index::l);        
    int nDS = dilNoise.dilutionSize(DistillationNoise<FImpl>::Index::s);     

    Grid::Coordinate coor  = gridHD->GlobalDimensions();
    coor[3] = nDL * nDS;
    Grid::GridCartesian * gridDD = Grid::SpaceTimeGrid::makeFourDimGrid(coor, GridDefaultSimd(Nd,vComplex::Nsimd()), GridDefaultMpi());

    envTmp(FermionField,    "fermionDDtmp_light" ,1, gridDD);
    envTmp(FermionField,    "fermionDDtmp_charm" ,1, gridDD);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TNormCheck<FImpl>::execute(void)
{
    // general grid setup
    GridCartesian * gridHD = envGetGrid(FermionField);
    GridCartesian * gridLD = envGetSliceGrid(FermionField,gridHD->Nd() -1);
    const int Ntlocal{gridHD->LocalDimensions()[Tdir]};
    const int Ntfirst{gridHD->LocalStarts()[Tdir]};
    int nT=env().getDim(Tdir);

    LOG(Message) << "WARNING: Assuming ordering s + ns*(l + nl*t) in DilutedNoise.hpp. This code will break when this changes!" << std::endl;

    // noise class -- assert they are identical and an "exact distillation" policy
    auto &dilNoise = envGet(DistillationNoise<FImpl>, par().noisePol);
    int nNoise = dilNoise.size(); 
    if(nNoise>1)
    {
        HADRONS_ERROR(Implementation, "NormCheck only implemented for exact distillation");
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

    // Temporary objects
    envGetTmp(FermionField,       fermion3dtmp1);

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

            std::vector<TComplex> ccBuf, clBuf;

            // read perambulator
            LOG(Message) << "Starting charm perambulator I/O for (tSrc,tSnk) = (" << tSrc << "," << tSnk << ")" << std::endl;
            envGetTmp(FermionField,    fermionDDtmp_charm);
            //startTimer("phi_c I/O");
            tFileName = par().vectorStemC;
            tFileName.append("_DD");
            tFileName.append("_tSm");
            tFileName.append(std::to_string(tSrc));
            tFileName.append("_tLoc");
            tFileName.append(std::to_string(tSnk));
            DistillationVectorsIo::readComponent(fermionDDtmp_charm, tFileName, 1, nDL, nDS, nDT, 0, vm().getTrajectory());
            //stopTimer("phi_c I/O");

            LOG(Message) << "    ++++++++++++++++++++++++++++++ " << std::endl;
            LOG(Message) << "Checking charm perambulator I/O for (tSrc,tSnk) = (" << tSrc << "," << tSnk << ")" << std::endl;
            for (int id1=0; id1<nDL*nDS; id1++)
            {
                ExtractSliceLocal(fermion3dtmp1, fermionDDtmp_charm, 0, id1, Tdir);
                RealD nrm = norm2(fermion3dtmp1);
                LOG(Message) << "Distillation index = " << id1 << ", Norm = " << nrm << std::endl;
            }
            LOG(Message) << "    ++++++++++++++++++++++++++++++ " << std::endl;

            // read perambulator
            LOG(Message) << "Starting light perambulator I/O for (tSrc,tSnk) = (" << tSrc << "," << tSnk << ")" << std::endl;
            envGetTmp(FermionField,    fermionDDtmp_light);
            //startTimer("phi_l I/O");
            tFileName = par().vectorStemL;
            tFileName.append("_DD");
            tFileName.append("_tSm");
            tFileName.append(std::to_string(tSrc));
            tFileName.append("_tLoc");
            tFileName.append(std::to_string(tSnk));
            DistillationVectorsIo::readComponent(fermionDDtmp_light, tFileName, 1, nDL, nDS, nDT, 0, vm().getTrajectory());
            //stopTimer("phi_l I/O");

            LOG(Message) << "    ++++++++++++++++++++++++++++++ " << std::endl;
            LOG(Message) << "Checking light perambulator I/O for (tSrc,tSnk) = (" << tSrc << "," << tSnk << ")" << std::endl;
            for (int id1=0; id1<nDL*nDS; id1++)
            {
                ExtractSliceLocal(fermion3dtmp1, fermionDDtmp_light, 0, id1, Tdir);
                RealD nrm = norm2(fermion3dtmp1);
                LOG(Message) << "Distillation index = " << id1 << ", Norm = " << nrm << std::endl;
            }
            LOG(Message) << "    ++++++++++++++++++++++++++++++ " << std::endl;

        }
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MDistil_NormCheck_hpp_
