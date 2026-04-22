#ifndef Hadrons_MDistil_DMeson4QuarkField_hpp_
#define Hadrons_MDistil_DMeson4QuarkField_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/TimerArray.hpp>
#include <Hadrons/Modules/MDistil/Base.hpp>

BEGIN_HADRONS_NAMESPACE

/********************************************************************************
 *                         DMeson4QuarkField                                    *
 * Computes the following sub-diagram:                                          *
 *                                                                              *
 *                                                                              *
 *           ____            ___                                                *
 *          /    \          /                                                   *
 *         /      \        /                                                    *
 *        /         v1  v3                                                      *
 *  M(rho1,rho2)                                                                *
 *        \         v2  v4                                                      *
 *         \      /        \                                                    *
 *          \____/          \___                                                *
 *                                                                              *
 *                                                                              *
 *  D(t=tD)           H_W(t)       tKpi                                         *
 *                                                                              *
 *                                                                              *
 *******************************************************************************/
BEGIN_MODULE_NAMESPACE(MDistil)

class DMeson4QuarkFieldPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(DMeson4QuarkFieldPar,
                                    std::string,                outPath,       // file stem for the out file
                                    std::string,                DMesonStem,    // file stem for the D
                                    std::string,                DMesonField,   // M(rho,rho) meson field for the D
                                    std::string,                vectorStemC,   // charm
                                    std::string,                vectorStemL,   // SU(3) light
                                    std::string,                noisePol,      // noise policy of v1 - assert compatibility with M(rho,rho) 
                                    unsigned int,               batchIO,       // parallel I/O yes/no?
                                    unsigned int,               fewerTH,       // compute all tH or just the ones between tD and tKpi?
                                    unsigned int,               tD,            // time of D meson
                                    unsigned int,               DMesSize,      // max time slice for DMeson rhos
                                    std::vector<unsigned int>,  tKpi,          // time of yet uncontracted part
                                    Gamma::Algebra,             gamma12,       // between vector1 and vector2
                                    Gamma::Algebra,             gamma34,       // between vector3 and vector4
                                    std::string,                mom,           // momentum injected into Hw
                                    unsigned int,               blockSize,     // tunable parameters
                                    unsigned int,               cacheSize);
};

template <typename FImpl>
class TDMeson4QuarkField: public Module<DMeson4QuarkFieldPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    class Result: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(Result,
                                        Gamma::Algebra, gamma_snk,
                                        Gamma::Algebra, gamma_src,
                                        std::vector<Complex>, corr);
    };
    // constructor
    TDMeson4QuarkField(const std::string name);
    // destructor
    virtual ~TDMeson4QuarkField(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(DMeson4QuarkField, TDMeson4QuarkField<FIMPL>, MDistil);

/******************************************************************************
 *                 TDMeson4QuarkField implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TDMeson4QuarkField<FImpl>::TDMeson4QuarkField(const std::string name)
: Module<DMeson4QuarkFieldPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TDMeson4QuarkField<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().noisePol};
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TDMeson4QuarkField<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TDMeson4QuarkField<FImpl>::setup(void)
{
    GridCartesian * gridHD = envGetGrid(FermionField);
    GridCartesian * gridLD = envGetSliceGrid(FermionField,gridHD->Nd() -1);
    
    envTmpLat(FermionField,      "fermion4dtmp");
    envTmp   (FermionField,      "fermion3dtmp1" ,1, gridLD);
    envTmp   (FermionField,      "fermion3dtmp2" ,1, gridLD);
    envTmp   (FermionField,      "fermion3dtmp3" ,1, gridLD);
    envTmp   (PropagatorField,   "prop3dtmp"     ,1, gridLD);
    envTmp   (ColourMatrixField, "MPhiPhi"       ,1, gridLD);
    envTmp   (ComplexField,      "cplx3dtmp"     ,1, gridLD);
    envTmpLat(ComplexField,      "ph");
    envTmp   (ComplexField,      "ph3d"          ,1, gridLD);
    envTmpLat(ComplexField,      "coor");


    auto &dilNoise = envGet(DistillationNoise<FImpl>, par().noisePol);
    int nDL = dilNoise.dilutionSize(DistillationNoise<FImpl>::Index::l);        
    int nDS = dilNoise.dilutionSize(DistillationNoise<FImpl>::Index::s);     

    envTmp    (std::vector<FermionField>, "vec_light", 1, nDL * nDS, gridLD);
    // maybe only initialise this optionally if batchIO == true ?
    envTmp    (std::vector<FermionField>, "vec_charm", 1, nDL * nDS, gridLD);
    
    envTmp(Vector<HADRONS_DISTIL_IO_TYPE>, "Sblock_buf", 1, nDL * nDS * nDL * nDS);
    envTmp(Vector<HADRONS_DISTIL_TYPE>,    "Scache_buf", 1, nDL * nDS * nDL * nDS);
    envTmp(Vector<HADRONS_DISTIL_IO_TYPE>, "Rblock_buf", 1, nDL * nDS * nDL * nDS);
    envTmp(Vector<HADRONS_DISTIL_TYPE>,    "Rcache_buf", 1, nDL * nDS * nDL * nDS);

    Grid::Coordinate coor  = gridHD->GlobalDimensions();
    coor[3] = nDL * nDS;
    Grid::GridCartesian * gridDD = Grid::SpaceTimeGrid::makeFourDimGrid( coor, GridDefaultSimd(Nd,vComplex::Nsimd()), GridDefaultMpi() );
    envTmp   (FermionField,    "fermionDDtmp_light" ,1, gridDD);
    envTmp   (FermionField,    "fermionDDtmp_charm" ,1, gridDD);



    
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TDMeson4QuarkField<FImpl>::execute(void)
{
    // general grid setup
    GridCartesian * gridHD = envGetGrid(FermionField);
    GridCartesian * gridLD = envGetSliceGrid(FermionField,gridHD->Nd() -1);
    const int Ntlocal{gridHD->LocalDimensions()[Tdir]};
    const int Ntfirst{gridHD->LocalStarts()[Tdir]};
    int nT=env().getDim(Tdir);
    unsigned int tD = par().tD; 
    
    // block and cache to store the output in
    envGetTmp(Vector<HADRONS_DISTIL_IO_TYPE>, Sblock_buf);
    envGetTmp(Vector<HADRONS_DISTIL_TYPE>, Scache_buf);
    envGetTmp(Vector<HADRONS_DISTIL_IO_TYPE>, Rblock_buf);
    envGetTmp(Vector<HADRONS_DISTIL_TYPE>, Rcache_buf);
    
    LOG(Message) << "Loaded block and cache" << std::endl;
    // read input D-meson field
    std::string mfPath = par().DMesonStem + "rho-rho." + std::to_string(vm().getTrajectory()) + "/" + par().DMesonField;   
    LOG(Message) << "reading " << mfPath << std::endl;
    TimerArray timer;
    std::vector<std::vector<int>> Dvec = {{(int)tD, (int)tD}};
    ContractionDistilMesonField<ComplexD,ComplexF> DMeson(mfPath, nT, timer, Dvec);
     
    // noise class -- assert they are identical and an "exact distillation" policy
    auto &dilNoise = envGet(DistillationNoise<FImpl>, par().noisePol);
    int nNoise = dilNoise.size(); 
    if(nNoise>1)
    {
        HADRONS_ERROR(Implementation, "DMeson4QuarkField only implemented for exact distillation");
    }
    int nDL = dilNoise.dilutionSize(DistillationNoise<FImpl>::Index::l);        
    int nDS = dilNoise.dilutionSize(DistillationNoise<FImpl>::Index::s);        
    int nDT = dilNoise.dilutionSize(DistillationNoise<FImpl>::Index::t);        
    // other input parameters
    std::vector<unsigned int> tKpi = par().tKpi; 
    if(tD>=nT)
    {
        HADRONS_ERROR(Range, "tD must be smaller than nT");
    }
    for(auto tKp : tKpi)
    {
        if(tKp>=nT)
        {
            HADRONS_ERROR(Range, "tKpi must be smaller than nT");
        }
    }
    Gamma                  g12(par().gamma12);
    Gamma                  g34(par().gamma34);
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

    // Temporary objects
    envGetTmp(FermionField,      fermion4dtmp);
    envGetTmp(FermionField,      fermion3dtmp1);
    envGetTmp(FermionField,      fermion3dtmp2);
    envGetTmp(FermionField,      fermion3dtmp3);
    envGetTmp(PropagatorField,   prop3dtmp);
    envGetTmp(ColourMatrixField, MPhiPhi);
    envGetTmp(ComplexField,      cplx3dtmp);
    
    int fewerTH = par().fewerTH;
    if(fewerTH && Ntlocal < nT)
    {
        LOG(Message) << "WARNING: Option 'fewerTH' is only optimised for trivial mpi layout in time-direction. You are likely wasting resources in this job!" << std::endl;
    }
    int batchIO = par().batchIO;
    envGetTmp(std::vector<FermionField>,    vec_light);
    envGetTmp(std::vector<FermionField>,    vec_charm);
    envGetTmp(FermionField,    fermionDDtmp_light);
    envGetTmp(FermionField,    fermionDDtmp_charm);



    // initialise file and metadata
    DistilMesonFieldMetadata<FImpl> md;
    for (auto pmu: p)
    {
        md.Momentum.push_back(pmu);
    }
    std::stringstream ss2;
    ss2 << par().gamma12 << "_" << par().gamma34;
    md.Operator          = ss2.str();
    md.Nt                = nT;   
    md.Nvec              = nDL;     //nvec=nDL for exact
    md.NoisePair         = {0,0};
    md.MesonFieldType    = "MD-Hw";
    md.RelativeSide      = "none";
    md.NoiseHashLeft     = "0";
    md.NoiseHashRight    = "0";
    //md.TimeDilutionLeft  = dilNoise.getMap();//[Index::t];
    //md.TimeDilutionRight = dilNoise.getMap();//[Index::t];
    //md.LapDilutionLeft   = index1[DistillationNoise<FImpl>::Index::l];
    //md.LapDilutionRight  = index1[DistillationNoise<FImpl>::Index::l];
    //md.SpinDilutionLeft  = index1[DistillationNoise<FImpl>::Index::s];
    //md.SpinDilutionRight = index1[DistillationNoise<FImpl>::Index::s];
   

    std::vector<DistilMatrixIo<HADRONS_DISTIL_IO_TYPE>> matrix_ioS(tKpi.size());
    std::vector<DistilMatrixIo<HADRONS_DISTIL_IO_TYPE>> matrix_ioR(tKpi.size());
    startTimer("file creation");
    // file name of output

    std::string DGamma = par().DMesonField.substr(0,par().DMesonField.find('.'));

    int iKpi=0;
    for(auto tKp : tKpi)
    {
        std::string outPathS = par().outPath, outPathR = par().outPath; 
        std::stringstream ss;
        ss << DGamma << "__" << par().gamma12 << "_" << par().gamma34 << "_p";
        for (unsigned int mu = 0; mu < p.size(); ++mu)
                ss << p[mu] << ((mu == p.size() - 1) ? "" : "_");
        //ss << ".h5";   
        std::string singletstr = ss.str()+"_singlet.h5", rearrstr = ss.str()+"_rearranged.h5";
        outPathS += "/D-Hw.tD" + std::to_string(tD) +".tKpi"+ std::to_string(tKp) + "." + std::to_string(vm().getTrajectory()) + "/" + singletstr;
        outPathR += "/D-Hw.tD" + std::to_string(tD) +".tKpi"+ std::to_string(tKp) + "." + std::to_string(vm().getTrajectory()) + "/" + rearrstr;
        
        makeFileDir(outPathS, gridHD);
        makeFileDir(outPathR, gridHD);
        unsigned int myRank = gridHD->ThisRank(); 
        DistilMatrixIo<HADRONS_DISTIL_IO_TYPE> mIOS(outPathS, DISTIL_MATRIX_NAME, nT, nDL * nDS, nDL * nDS);
        DistilMatrixIo<HADRONS_DISTIL_IO_TYPE> mIOR(outPathR, DISTIL_MATRIX_NAME, nT, nDL * nDS, nDL * nDS);
        if(myRank==0)
        {
            mIOS.initFile(md);
            mIOR.initFile(md);
        }
        gridHD->Barrier();
        matrix_ioS[iKpi] = mIOS;
        matrix_ioR[iKpi] = mIOR;
        iKpi++;
    }
    stopTimer("file creation");
    // determine timeslices tH which are between tD and tKpi (shorter distance)
    std::vector<std::vector<unsigned int>> tHs;
    std::vector<unsigned int> tHs_flat;
    for(auto tKp : tKpi)
    {
        std::vector<unsigned int> tH_iKpi;
        int tDMinusTKpi = (tD - tKp + nT) % nT;
        int tKpiMinusTD = (tKp - tD + nT) % nT;
        if(tDMinusTKpi < tKpiMinusTD)
        {
            for(int iTH = 1; iTH < tDMinusTKpi; iTH++)
            {
                int tH_tmp = (tKp + iTH + nT) % nT;
                tH_iKpi.push_back(tH_tmp);      
                if( std::find(tHs_flat.begin(), tHs_flat.end(), tH_tmp) == tHs_flat.end())
                {
                    tHs_flat.push_back(tH_tmp);     
                } 
            }     
            tHs.push_back(tH_iKpi);
        }
        else
        {
            for(int iTH = 1; iTH < tKpiMinusTD; iTH++)
            {
                int tH_tmp = (tD + iTH + nT) % nT;
                tH_iKpi.push_back(tH_tmp);      
                if( std::find(tHs_flat.begin(), tHs_flat.end(), tH_tmp) == tHs_flat.end())
                {
                    tHs_flat.push_back(tH_tmp);     
                } 
            }     
            tHs.push_back(tH_iKpi);
        }
    }
    
    LOG(Message) << "WARNING: Assuming ordering s + ns*(l + nl*t) in DilutedNoise.hpp. This code will break when this changes!" << std::endl;
    // variables used in the loop structure
    int dk1,ds1,dk2,ds2,dSolve1,dSolve2,tH;
    std::array<unsigned int, 3> index1,index2;
    std::vector<TComplex>  Sbuf, Rbuf;
    std::string tFileName;

    const uint i_rank =  gridHD->ThisRank();
    const uint N_ranks = gridHD->RankCount(); 
 
 
    // loop over tH
    for (int t = 0; t < Ntlocal; t++ )
    {
        tH = t + Ntfirst;
        if(fewerTH && (std::find(tHs_flat.begin(), tHs_flat.end(), tH) == tHs_flat.end()) )
        {
            LOG(Message) << "Only computing three-point function for timeslices between tD and tKpi, skipping tH = " << tH << std::endl;
            continue;
        } 
        MPhiPhi=Zero();    
        // 3D phase e^{ipx}
        ExtractSliceLocal(ph3d,ph,0,t,Tdir);  
        // initialise vector 2
        if(batchIO)
        {
            startTimer("MPhiPhi I/O: light");
            LOG(Message) << "Starting batch I/O" << std::endl;
            /*tFileName = par().vectorStemL;
            tFileName.append("_tSm");
            tFileName.append(std::to_string(tD));
            tFileName.append("_tLoc");
            tFileName.append(std::to_string(tH));
            DistillationVectorsIo::read(vec_light, tFileName, 1, nDL, nDS, nDT, false, vm().getTrajectory());*/
            tFileName = par().vectorStemL;
            tFileName.append("_DD");
            tFileName.append("_tSm");
            tFileName.append(std::to_string(tD));
            tFileName.append("_tLoc");
            tFileName.append(std::to_string(tH));
            DistillationVectorsIo::readComponent(fermionDDtmp_light, tFileName, 1, nDL, nDS, nDT, 0, vm().getTrajectory());
            stopTimer("MPhiPhi I/O: light");
            startTimer("MPhiPhi I/O: charm");
            LOG(Message) << "Starting batch I/O" << std::endl;
            /*tFileName = par().vectorStemC;
            tFileName.append("_tSm");
            tFileName.append(std::to_string(tD));
            tFileName.append("_tLoc");
            tFileName.append(std::to_string(tH));
            DistillationVectorsIo::read(vec_charm, tFileName, 1, nDL, nDS, nDT, false, vm().getTrajectory());*/
            tFileName = par().vectorStemC;
            tFileName.append("_DD");
            tFileName.append("_tSm");
            tFileName.append(std::to_string(tD));
            tFileName.append("_tLoc");
            tFileName.append(std::to_string(tH));
            DistillationVectorsIo::readComponent(fermionDDtmp_charm, tFileName, 1, nDL, nDS, nDT, 0, vm().getTrajectory());
            stopTimer("MPhiPhi I/O: charm");
        }
        else
        {
            startTimer("MPhiPhi I/O: light");
            for(int id2=0; id2<nDL * nDS; id2++)
            {
                index2 = dilNoise.dilutionCoordinates(id2);  
                dk2 = index2[DistillationNoise<FImpl>::Index::l];
                ds2 = index2[DistillationNoise<FImpl>::Index::s];
                dSolve2 = dilNoise.dilutionIndex(tD,dk2,ds2);
                tFileName = par().vectorStemL;
                tFileName.append("_t");
                tFileName.append(std::to_string(tH));
                DistillationVectorsIo::readComponent(fermion3dtmp2, tFileName, 1, nDL, nDS, nDT, dSolve2, vm().getTrajectory());
                // this is vector 2 on timeslice tH 
                vec_light[id2]=fermion3dtmp2;
            }
            stopTimer("MPhiPhi I/O: light");
        }
        for(int id1=0; id1<nDL * nDS; id1++)
        {
            // this line is where the ordering s + ns*(l + nl*t) is assumed
            index1 = dilNoise.dilutionCoordinates(id1);  
            dk1 = index1[DistillationNoise<FImpl>::Index::l];
            ds1 = index1[DistillationNoise<FImpl>::Index::s];
            dSolve1 = dilNoise.dilutionIndex(tD,dk1,ds1);
            if(batchIO)
            {
                //fermion3dtmp1 = vec_charm[id1];
                ExtractSliceLocal(fermion3dtmp1, fermionDDtmp_charm,0,id1,Tdir);
            }
            else
            {
                startTimer("MPhiPhi I/O: charm");
                // this is vector 1 on timeslice tH 
                tFileName = par().vectorStemC;
                tFileName.append("_t");
                tFileName.append(std::to_string(tH));
                DistillationVectorsIo::readComponent(fermion3dtmp1, tFileName, 1, nDL, nDS, nDT, dSolve1, vm().getTrajectory());
                stopTimer("MPhiPhi I/O: charm");
            }
            startTimer("computation MPhiPhi");
            for(int id2=0; id2<nDL * nDS; id2++)
            {
                //fermion3dtmp3 = g12*vec_light[id2];
                ExtractSliceLocal(fermion3dtmp2, fermionDDtmp_light,0,id2,Tdir);
                fermion3dtmp3 = g12*fermion3dtmp2;
                fermion3dtmp2 = fermion3dtmp3*DMeson(tD,tD,tD)(id2,id1);
                prop3dtmp = outerProductC(fermion3dtmp1,fermion3dtmp2);
                // this object is sum_{spin,d1,d2} (DMeson[d1,d2] * vector1[d1] * gamma12 * vector2[d2]) on timeslice tH
                MPhiPhi += traceSpin(prop3dtmp);
            }
            stopTimer("computation MPhiPhi");
        }
        iKpi=0;
        for(auto tKp : tKpi)
        {
            if(fewerTH && (std::find(tHs[iKpi].begin(), tHs[iKpi].end(), tH) == tHs[iKpi].end()) )
            {
                LOG(Message) << "Only computing three-point function for timeslices between tD and tKpi, skipping tH = " << tH << " for tKpi = " << tKp << " and tD = " << tD << std::endl;
                iKpi++;
                continue;
            } 
            startTimer("K-Pi I/O: light");
            if(batchIO)
            {
                LOG(Message) << "Starting batch I/O" << std::endl;
                /*std::string tFileName = par().vectorStemL;
                tFileName.append("_tSm"); 
                tFileName.append(std::to_string(tKp));  
                tFileName.append("_tLoc");
                tFileName.append(std::to_string(tH));
                DistillationVectorsIo::read(vec_light, tFileName, 1, nDL, nDS, nDT, false, vm().getTrajectory());*/
                tFileName = par().vectorStemL;
                tFileName.append("_DD");
                tFileName.append("_tSm");
                tFileName.append(std::to_string(tKp));
                tFileName.append("_tLoc");
                tFileName.append(std::to_string(tH));
                DistillationVectorsIo::readComponent(fermionDDtmp_light, tFileName, 1, nDL, nDS, nDT, 0, vm().getTrajectory());
            } 
            else
            {
                for(int id2=0; id2<nDL * nDS; id2++)
                {
                    index2 = dilNoise.dilutionCoordinates(id2);  
                    dk2 = index2[DistillationNoise<FImpl>::Index::l];
                    ds2 = index2[DistillationNoise<FImpl>::Index::s];
                    dSolve2 = dilNoise.dilutionIndex(tKp,dk2,ds2);
                    // this is vector 2 on timeslice tH 
                    std::string tFileName = par().vectorStemL;
                    tFileName.append("_t");
                    tFileName.append(std::to_string(tH));
                    DistillationVectorsIo::readComponent(fermion3dtmp2, tFileName, 1, nDL, nDS, nDT, dSolve2, vm().getTrajectory());
                    vec_light[id2]=fermion3dtmp2;
                }
            }
            stopTimer("K-Pi I/O: light");
            /*************************************************
            FE: checked here that a sliceSum over MPhiPhi
            reproduces exactly a contraction of meson fields
            tr[ M(rho,rho; tD,tD,tD) * M(phi,phi; tD,tD,tD) ]
            *************************************************/
            startTimer("computation D-4quark");
            DistilMatrixSetIo<ComplexF> Sblock(Sblock_buf.data(), 1 , 1, nDL * nDS, nDL * nDS);
            DistilMatrixSetIo<ComplexF> Rblock(Rblock_buf.data(), 1 , 1, nDL * nDS, nDL * nDS);
            for(int id1=0; id1<nDL * nDS; id1++)
            {
                //fermion3dtmp1 = vec_light[id1];
                ExtractSliceLocal(fermion3dtmp1, fermionDDtmp_light,0,id1,Tdir);
                for(int id2=0; id2<nDL * nDS; id2++)
                {
                    // no caching for the moment - but keep this here in case anyone wants to optimise this code at some stage
                    DistilMatrixSetCache<ComplexD> Scache(Scache_buf.data(), 1, 1, 1, 1, 1);
                    DistilMatrixSetCache<ComplexD> Rcache(Rcache_buf.data(), 1, 1, 1, 1, 1);
                    //fermion3dtmp3 = g34*vec_light[id2];          
                    ExtractSliceLocal(fermion3dtmp2, fermionDDtmp_light,0,id2,Tdir);
                    fermion3dtmp3 = g34*fermion3dtmp2;
                    fermion3dtmp2 = fermion3dtmp3;
                    prop3dtmp = outerProductC(fermion3dtmp1,fermion3dtmp2);

                    // colour-singlet -> two colour traces
                    cplx3dtmp = trace(prop3dtmp)*traceColour(MPhiPhi)*ph3d;
                    sliceSum(cplx3dtmp,Sbuf,Tdir);
                    Scache(0,0,0,0,0)=TensorRemove(Sbuf[0]);
                    Sblock(0,0,id1,id2) = Scache(0,0,0,0,0);                

                    // colour-rearranged -> one colour trace
                    cplx3dtmp = traceColour(traceSpin(prop3dtmp)*MPhiPhi)*ph3d;
                    sliceSum(cplx3dtmp,Rbuf,Tdir);                
                    Rcache(0,0,0,0,0)=TensorRemove(Rbuf[0]);
                    Rblock(0,0,id1,id2) = Rcache(0,0,0,0,0);                
                }
            }
            stopTimer("computation D-4quark");
            startTimer("serial write I/O");
            LOG(Message) << "Starting serial IO for tH = " << tH << std::endl;
            DistilMatrixSetTimeSliceIo<ComplexF> Sblock_relative(Sblock_buf.data(), 1, nDL * nDS, nDL * nDS);
            DistilMatrixSetTimeSliceIo<ComplexF> Rblock_relative(Rblock_buf.data(), 1, nDL * nDS, nDL * nDS);
            std::string dataset_name = std::to_string(tKp)+"-"+std::to_string(tKp);
            gridHD->Barrier();
            for(int iIO=0; iIO<N_ranks; iIO++)
            {
                if(iIO==i_rank)
                {
                    LOG(Message) << "Writing from rank " << i_rank << std::endl;
                    matrix_ioS[iKpi].saveBlock(Sblock_relative, 0, 0, 0, dataset_name, 0, nDL * nDS, std::to_string(tH));
                    matrix_ioR[iKpi].saveBlock(Rblock_relative, 0, 0, 0, dataset_name, 0, nDL * nDS, std::to_string(tH));
                }
                gridHD->Barrier();
            }
            stopTimer("serial write I/O");
            iKpi++;
        }
    }
    
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MDistil_DMeson4QuarkField_hpp_
