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
 *  The output F(id1,id2) has open distillation indices at tKpi only; colour    *
 *  is fully contracted at the H_W vertex. The colour labels refer to the       *
 *  contraction pairing used here, (12) = {charm, light from D} and             *
 *  (34) = {two light lines to tKpi}:                                           *
 *  singlet    : colour traced within (12) and within (34)                      *
 *               -> two colour traces                                           *
 *  rearranged : colour runs across (12) and (34)                               *
 *               -> one colour trace                                            *
 *  Mapping to Q1/Q2 (Fierz) requires the full (V-A)x(V-A) set of gamma pairs.  *
 *                                                                              *
 *  Perambulators are read from the raw "<stem>.<traj>/t<tH>_pkg.bin" files.    *
 *                                                                              *
 *******************************************************************************/
BEGIN_MODULE_NAMESPACE(MDistil)

typedef std::pair<Gamma::Algebra, Gamma::Algebra> GammaPair;

class DMeson4QuarkFieldPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(DMeson4QuarkFieldPar,
                                    std::string,                outPath,        // file stem for the out files
                                    std::string,                RhoRhoStem,     // file stem for the rho-rho MFs
                                    std::string,                DMesonGamma,    // gamma of the M(rho,rho) meson field for the D
                                    std::vector<std::string>,   momsD,          // list of momenta of the D meson field
                                    std::string,                vectorStemC,    // charm
                                    std::string,                vectorStemL,    // SU(3) light
                                    std::string,                rawTimeSources, // time sources of the raw data
                                    std::string,                noisePol,       // noise policy of v1 - assert compatibility with M(rho,rho)
                                    unsigned int,               fewerTH,        // compute all tH or just the ones between tD and tKpi?
                                    std::vector<unsigned int>,  tDs,            // times of D meson
                                    std::vector<unsigned int>,  tKpis,          // times of yet uncontracted part
                                    std::string,                gammas,         // list of space-separated pairs of gamma matrices (g12 g34)
                                    std::string,                mom);           // momentum injected into Hw
};

template <typename FImpl>
class TDMeson4QuarkField: public Module<DMeson4QuarkFieldPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
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

    envTmp   (FermionField,      "fermion3dtmp"  ,1, gridLD);
    envTmp   (FermionField,      "chiS"          ,1, gridLD);
    envTmp   (FermionField,      "chiR"          ,1, gridLD);
    envTmp   (PropagatorField,   "prop3dtmp"     ,1, gridLD);
    envTmp   (ColourMatrixField, "MPhiPhi"       ,1, gridLD);
    envTmpLat(ComplexField,      "ph");
    envTmp   (ComplexField,      "ph3d"          ,1, gridLD);
    envTmpLat(ComplexField,      "coor");

    auto &dilNoise = envGet(DistillationNoise<FImpl>, par().noisePol);
    int nDL = dilNoise.dilutionSize(DistillationNoise<FImpl>::Index::l);
    int nDS = dilNoise.dilutionSize(DistillationNoise<FImpl>::Index::s);

    envTmp(std::vector<FermionField>, "charmSlices", 1, nDL * nDS, gridLD);
    envTmp(std::vector<FermionField>, "lightSlices", 1, nDL * nDS, gridLD);
    // holds psi on the D side and gamma34 * phi_l on the K-pi side
    envTmp(std::vector<FermionField>, "workSlices",  1, nDL * nDS, gridLD);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TDMeson4QuarkField<FImpl>::execute(void)
{
    if (par().vectorStemC.empty() || par().vectorStemL.empty())
    {
        HADRONS_ERROR(Argument, "DMeson4QuarkField requires both vectorStemC and vectorStemL");
    }

    // general grid setup
    GridCartesian * gridHD = envGetGrid(FermionField);
    GridCartesian * gridLD = envGetSliceGrid(FermionField,gridHD->Nd() -1);
    const int Ntlocal{gridHD->LocalDimensions()[Tdir]};
    const int Ntfirst{gridHD->LocalStarts()[Tdir]};
    const int nT = env().getDim(Tdir);
    const std::string traj = std::to_string(vm().getTrajectory());

    // noise class -- assert "exact distillation" policy
    auto &dilNoise = envGet(DistillationNoise<FImpl>, par().noisePol);
    int nNoise = dilNoise.size();
    if(nNoise>1)
    {
        HADRONS_ERROR(Implementation, "DMeson4QuarkField only implemented for exact distillation");
    }
    const int nDL = dilNoise.dilutionSize(DistillationNoise<FImpl>::Index::l);
    const int nDS = dilNoise.dilutionSize(DistillationNoise<FImpl>::Index::s);
    const int nDT = dilNoise.dilutionSize(DistillationNoise<FImpl>::Index::t);
    const int nDil = nDL * nDS;
    LOG(Message) << "WARNING: Assuming ordering s + ns*(l + nl*t) in DilutedNoise.hpp. This code will break when this changes!" << std::endl;

    // time slices
    const std::vector<unsigned int> tDs = par().tDs, tKpis = par().tKpis;
    if (tDs.empty() || tKpis.empty())
    {
        HADRONS_ERROR(Argument, "DMeson4QuarkField requires at least one tD and one tKpi");
    }
    for(auto tD : tDs)
    {
        if(tD>=nT)
        {
            HADRONS_ERROR(Range, "all tDs must be smaller than nT");
        }
    }
    for(auto tKpi : tKpis)
    {
        if(tKpi>=nT)
        {
            HADRONS_ERROR(Range, "all tKpis must be smaller than nT");
        }
    }
    const bool fewerTH = (par().fewerTH != 0);
    if(fewerTH && Ntlocal < nT)
    {
        LOG(Message) << "WARNING: Option 'fewerTH' is only optimised for trivial mpi layout in time-direction. You are likely wasting resources in this job!" << std::endl;
    }

    // time sources contained in the raw perambulator files
    std::vector<int> rawTimeSources;
    if(!par().rawTimeSources.empty())
    {
        rawTimeSources = strToVec<int>(par().rawTimeSources);
    }
    else
    {
        rawTimeSources.resize(nDT);
        std::iota(rawTimeSources.begin(), rawTimeSources.end(), 0);
    }
    auto getRawSourceOffset = [&](const unsigned int tSm)
    {
        auto sourceIt = std::find(rawTimeSources.begin(), rawTimeSources.end(), static_cast<int>(tSm));
        if(sourceIt == rawTimeSources.end())
        {
            HADRONS_ERROR(Io, "raw 3D field stem does not contain requested source time " + std::to_string(tSm));
        }
        return static_cast<int>(sourceIt - rawTimeSources.begin());
    };
    // read sources in file order, so that a single forward pass through each file suffices
    auto sourceOrder = [&](const std::vector<unsigned int> &times)
    {
        std::vector<unsigned int> order(times.size());
        std::iota(order.begin(), order.end(), 0);
        std::stable_sort(order.begin(), order.end(), [&](const unsigned int a, const unsigned int b)
        {
            return getRawSourceOffset(times[a]) < getRawSourceOffset(times[b]);
        });
        return order;
    };
    const std::vector<unsigned int> tDOrder = sourceOrder(tDs), tKpiOrder = sourceOrder(tKpis);

    // D meson momenta
    std::vector<std::string> Dmoms;
    for (auto dmom : par().momsD)
    {
        std::vector<int> dmomI = strToVec<int>(dmom);
        if (dmomI.size() != 3)
        {
            HADRONS_ERROR(Size, "each D meson momentum must have 3 components");
        }
        Dmoms.push_back(std::to_string(dmomI[0]) + "_" + std::to_string(dmomI[1]) + "_" + std::to_string(dmomI[2]));
    }
    if (Dmoms.empty())
    {
        HADRONS_ERROR(Argument, "DMeson4QuarkField requires at least one D meson momentum");
    }

    // gamma pairs, and the distinct gamma12 and gamma34 among them
    const std::vector<GammaPair> gammas = strToVec<GammaPair>(par().gammas);
    if (gammas.empty())
    {
        HADRONS_ERROR(Argument, "DMeson4QuarkField requires at least one gamma pair");
    }
    std::vector<Gamma::Algebra> g12s, g34s;
    std::vector<unsigned int>   g12Index, g34Index;
    auto distinctIndex = [](std::vector<Gamma::Algebra> &list, const Gamma::Algebra g)
    {
        auto it = std::find(list.begin(), list.end(), g);
        if (it == list.end())
        {
            list.push_back(g);
            return static_cast<unsigned int>(list.size() - 1);
        }
        return static_cast<unsigned int>(it - list.begin());
    };
    for (auto &gp : gammas)
    {
        g12Index.push_back(distinctIndex(g12s, gp.first));
        g34Index.push_back(distinctIndex(g34s, gp.second));
    }

    const unsigned int nTD = tDs.size(), nTKpi = tKpis.size(), nPD = Dmoms.size(), nGP = gammas.size();

    // read the D-meson fields, only the diagonal time slices tD needed by this job
    std::vector<std::vector<DistilMesonFieldMatrix<ComplexD>>> DMat(nPD);
    std::vector<std::vector<int>> tVec;
    for (auto tD : tDs)
    {
        tVec.push_back({static_cast<int>(tD), static_cast<int>(tD)});
    }
    startTimer("MesonField IO");
    for (unsigned int pDi = 0; pDi < nPD; pDi++)
    {
        std::string mfPath = par().RhoRhoStem + "rho-rho." + traj + "/" + par().DMesonGamma + "_p" + Dmoms[pDi] + ".h5";
        LOG(Message) << "reading " << mfPath << " for " << tVec.size() << " diagonal time slices" << std::endl;
        TimerArray timer;
        ContractionDistilMesonField<ComplexD,ComplexF> DMeson(mfPath, nT, timer, tVec);
        for (auto tD : tDs)
        {
            DMat[pDi].push_back(DMeson(tD,tD,tD));
            if (DMat[pDi].back().rows() != nDil || DMat[pDi].back().cols() != nDil)
            {
                HADRONS_ERROR(Size, "D meson field dimensions do not match the dilution size");
            }
        }
    }
    stopTimer("MesonField IO");

    // momentum phase e^{ipx} for Hw
    Complex           i(0.0,1.0);
    std::vector<Real> p;
    p  = strToVec<Real>(par().mom);
    if (p.size() != env().getNd())
    {
        HADRONS_ERROR(Size, "Hw momentum must have " + std::to_string(env().getNd()) + " components");
    }
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

    // output files, one per (tD, tKpi, D momentum, gamma pair) and colour structure
    auto blockIndex = [&](const unsigned int tDi, const unsigned int pDi, const unsigned int gpi,
                          const unsigned int tKpii)
    {
        return ((tDi*nPD + pDi)*nGP + gpi)*nTKpi + tKpii;
    };
    const unsigned int nBlocks = nTD*nPD*nGP*nTKpi;
    std::vector<DistilMatrixIo<HADRONS_DISTIL_IO_TYPE>> matrix_ioS(nBlocks), matrix_ioR(nBlocks);
    startTimer("file creation");
    for (unsigned int tDi = 0; tDi < nTD; tDi++)
    for (unsigned int pDi = 0; pDi < nPD; pDi++)
    for (unsigned int gpi = 0; gpi < nGP; gpi++)
    {
        DistilMesonFieldMetadata<FImpl> md;
        for (auto pmu: p)
        {
            md.Momentum.push_back(pmu);
        }
        std::stringstream ssOp;
        ssOp << gammas[gpi].first << "_" << gammas[gpi].second;
        md.Operator          = ssOp.str();
        md.Nt                = nT;
        md.Nvec              = nDL;     //nvec=nDL for exact
        md.NoisePair         = {0,0};
        md.MesonFieldType    = "MD-Hw";
        md.RelativeSide      = "none";
        md.NoiseHashLeft     = "0";
        md.NoiseHashRight    = "0";

        std::stringstream ss;
        ss << par().DMesonGamma << "_p" << Dmoms[pDi] << "__" << gammas[gpi].first << "_" << gammas[gpi].second << "_p";
        for (unsigned int mu = 0; mu < p.size(); ++mu)
                ss << p[mu] << ((mu == p.size() - 1) ? "" : "_");
        for (unsigned int tKpii = 0; tKpii < nTKpi; tKpii++)
        {
            std::string dir = par().outPath + "/D-Hw.tD" + std::to_string(tDs[tDi]) + ".tKpi" + std::to_string(tKpis[tKpii]) + "." + traj + "/";
            std::string outPathS = dir + ss.str() + "_singlet.h5", outPathR = dir + ss.str() + "_rearranged.h5";
            makeFileDir(outPathS, gridHD);
            makeFileDir(outPathR, gridHD);
            const unsigned int b = blockIndex(tDi, pDi, gpi, tKpii);
            matrix_ioS[b] = DistilMatrixIo<HADRONS_DISTIL_IO_TYPE>(outPathS, DISTIL_MATRIX_NAME, nT, nDil, nDil);
            matrix_ioR[b] = DistilMatrixIo<HADRONS_DISTIL_IO_TYPE>(outPathR, DISTIL_MATRIX_NAME, nT, nDil, nDil);
            if(gridHD->IsBoss())
            {
                matrix_ioS[b].initFile(md);
                matrix_ioR[b].initFile(md);
            }
        }
    }
    gridHD->Barrier();
    stopTimer("file creation");

    // is tH needed for the pair (tD,tKpi)? with fewerTH only the timeslices strictly
    // between tD and tKpi (shorter distance) are computed
    auto needsTH = [&](const int tD, const int tKpi, const int tH)
    {
        if (!fewerTH)
        {
            return true;
        }
        const int tDMinusTKpi = (tD - tKpi + nT) % nT;
        const int tKpiMinusTD = (tKpi - tD + nT) % nT;
        if (tDMinusTKpi < tKpiMinusTD)
        {
            const int tHMinusTKpi = (tH - tKpi + nT) % nT;
            return (tHMinusTKpi >= 1) && (tHMinusTKpi < tDMinusTKpi);
        }
        const int tHMinusTD = (tH - tD + nT) % nT;
        return (tHMinusTD >= 1) && (tHMinusTD < tKpiMinusTD);
    };

    // raw perambulator reader; records are ordered [source offset][dilution index]
    auto readRawDDSlices = [&](std::vector<FermionField> &slices, ScidacReader &reader,
                               const std::string &filename, int &currentSourceOffset,
                               const unsigned int tSm)
    {
        const int sourceOffset = getRawSourceOffset(tSm);
        if ((currentSourceOffset < 0) || (sourceOffset < currentSourceOffset))
        {
            if (currentSourceOffset >= 0)
            {
                reader.close();
            }
            reader.open(filename);
            currentSourceOffset = 0;
        }

        for (int i = 0; i < (sourceOffset - currentSourceOffset) * nDil; ++i)
        {
            reader.skipScidacFieldRecord();
        }

        for (int id = 0; id < nDil; ++id)
        {
            std::array<unsigned int, 3> index = dilNoise.dilutionCoordinates(id);
            const int dk = index[DistillationNoise<FImpl>::Index::l];
            const int ds = index[DistillationNoise<FImpl>::Index::s];
            const int dSolve = dilNoise.dilutionIndex(tSm, dk, ds);
            DistillationVectorsIo::Record record;

            reader.readScidacFieldRecord(slices[id], record);
            if (record.index != static_cast<unsigned int>(dSolve))
            {
                HADRONS_ERROR(Io, "vector index mismatch: record.index = " + std::to_string(record.index) +
                                       " and componentIndex = " + std::to_string(dSolve));
            }
            if (record.nNoise != 1 || record.nDL != nDL || record.nDS != nDS || record.nDT != nDT)
            {
                HADRONS_ERROR(Io, "dilution parameter mismatch");
            }
        }
        currentSourceOffset = sourceOffset + 1;
    };

    // Temporary objects
    envGetTmp(FermionField,              fermion3dtmp);
    envGetTmp(FermionField,              chiS);
    envGetTmp(FermionField,              chiR);
    envGetTmp(PropagatorField,           prop3dtmp);
    envGetTmp(ColourMatrixField,         MPhiPhi);
    envGetTmp(std::vector<FermionField>, charmSlices);
    envGetTmp(std::vector<FermionField>, lightSlices);
    envGetTmp(std::vector<FermionField>, workSlices);

    // D-side objects per (tD, D momentum, gamma12) on the current tH, with the Hw phase folded in:
    //   MTph = e^{ipx} MPhiPhi^T, wph = e^{ipx} tr_c[MPhiPhi]
    auto dIndex = [&](const unsigned int tDi, const unsigned int pDi, const unsigned int g12i)
    {
        return (tDi*nPD + pDi)*g12s.size() + g12i;
    };
    std::vector<ColourMatrixField> MTph(nTD*nPD*g12s.size(), ColourMatrixField(gridLD));
    std::vector<ComplexField>      wph(nTD*nPD*g12s.size(),  ComplexField(gridLD));

    // results of the current tH: rank-local partial sums, reduced over the slice before writing
    const size_t blockSize = static_cast<size_t>(nDil)*nDil;
    std::vector<ComplexD>                Sbuf(nBlocks*blockSize), Rbuf(nBlocks*blockSize);
    std::vector<HADRONS_DISTIL_IO_TYPE>  ioBuf(blockSize);
    std::vector<char>                    blockDone(nBlocks);
    const int nProcT  = gridHD->ProcessorGrid()[Tdir];
    const int myProcT = gridHD->ThisProcessorCoor()[Tdir];

    unsigned int charmReads = 0, lightReads = 0, skippedTH = 0;

    // loop over tH
    for (int t = 0; t < Ntlocal; t++)
    {
        const int tH = t + Ntfirst;

        // which tD and tKpi are needed on this tH
        std::vector<char> tDNeeded(nTD, 0), tKpiNeeded(nTKpi, 0);
        for (unsigned int tDi = 0; tDi < nTD; tDi++)
        for (unsigned int tKpii = 0; tKpii < nTKpi; tKpii++)
        {
            if (needsTH(tDs[tDi], tKpis[tKpii], tH))
            {
                tDNeeded[tDi]     = 1;
                tKpiNeeded[tKpii] = 1;
            }
        }
        std::fill(blockDone.begin(), blockDone.end(), 0);
        const bool anyWork = std::any_of(tDNeeded.begin(), tDNeeded.end(), [](const char c){ return c != 0; });

        if (!anyWork)
        {
            skippedTH++;
            LOG(Message) << "Only computing three-point function for timeslices between tD and tKpi, skipping tH = " << tH << std::endl;
        }
        else
        {
            // 3D phase e^{ipx}
            ExtractSliceLocal(ph3d,ph,0,t,Tdir);

            const std::string charmFilename = par().vectorStemC + "." + traj + "/t" + std::to_string(tH) + "_pkg.bin";
            const std::string lightFilename = par().vectorStemL + "." + traj + "/t" + std::to_string(tH) + "_pkg.bin";
            ScidacReader charmReader, lightReader;
            int charmSourceOffset = -1, lightSourceOffset = -1;

            // D side: sum_{spin,d1,d2} (DMeson[d2,d1] * vector1[d1] * gamma12 * vector2[d2]) on timeslice tH
            for (auto tDi : tDOrder)
            {
                if (!tDNeeded[tDi])
                {
                    continue;
                }
                const unsigned int tD = tDs[tDi];
                LOG(Message) << "Starting perambulator I/O for (tD,tH) = (" << tD << "," << tH << ")" << std::endl;
                startTimer("phi_c I/O");
                readRawDDSlices(charmSlices, charmReader, charmFilename, charmSourceOffset, tD);
                stopTimer("phi_c I/O");
                charmReads++;
                startTimer("phi_l I/O");
                readRawDDSlices(lightSlices, lightReader, lightFilename, lightSourceOffset, tD);
                stopTimer("phi_l I/O");
                lightReads++;

                startTimer("computation MPhiPhi");
                for (unsigned int pDi = 0; pDi < nPD; pDi++)
                {
                    const auto &D = DMat[pDi][tDi];
                    // psi[d1] = sum_{d2} vector2[d2] * DMeson[d2,d1]
                    for (int id1 = 0; id1 < nDil; id1++)
                    {
                        workSlices[id1] = Zero();
                        for (int id2 = 0; id2 < nDil; id2++)
                        {
                            workSlices[id1] += lightSlices[id2]*D(id2,id1);
                        }
                    }
                    for (unsigned int g12i = 0; g12i < g12s.size(); g12i++)
                    {
                        Gamma g12(g12s[g12i]);
                        MPhiPhi = Zero();
                        for (int id1 = 0; id1 < nDil; id1++)
                        {
                            fermion3dtmp = g12*workSlices[id1];
                            // outerProductC(l, r) = conj(l)*r
                            prop3dtmp = outerProductC(charmSlices[id1], fermion3dtmp);
                            MPhiPhi += traceSpin(prop3dtmp);
                        }
                        const unsigned int di = dIndex(tDi, pDi, g12i);
                        MTph[di] = ph3d*transpose(MPhiPhi);
                        wph[di]  = ph3d*traceColour(MPhiPhi);
                    }
                }
                stopTimer("computation MPhiPhi");
            }
            charmReader.close();

            /*************************************************
            FE: checked here that a sliceSum over MPhiPhi
            reproduces exactly a contraction of meson fields
            tr[ M(rho,rho; tD,tD,tD) * M(phi,phi; tD,tD,tD) ]
            *************************************************/

            // K-pi side, with prop = outerProductC(vector3[d1], gamma34 * vector4[d2]):
            //   singlet    : sum_x tr[prop] tr_c[MPhiPhi] e^{ipx} = <vector3[d1] | wph * gamma34 * vector4[d2]>
            //   rearranged : sum_x tr_c[tr_s[prop] MPhiPhi] e^{ipx} = <vector3[d1] | MTph * gamma34 * vector4[d2]>
            for (auto tKpii : tKpiOrder)
            {
                if (!tKpiNeeded[tKpii])
                {
                    continue;
                }
                const unsigned int tKpi = tKpis[tKpii];
                LOG(Message) << "Starting light perambulator I/O for (tKpi,tH) = (" << tKpi << "," << tH << ")" << std::endl;
                startTimer("phi_l I/O");
                readRawDDSlices(lightSlices, lightReader, lightFilename, lightSourceOffset, tKpi);
                stopTimer("phi_l I/O");
                lightReads++;

                startTimer("computation D-4quark");
                for (unsigned int g34i = 0; g34i < g34s.size(); g34i++)
                {
                    Gamma g34(g34s[g34i]);
                    for (int id = 0; id < nDil; id++)
                    {
                        workSlices[id] = g34*lightSlices[id];
                    }
                    for (unsigned int gpi = 0; gpi < nGP; gpi++)
                    {
                        if (g34Index[gpi] != g34i)
                        {
                            continue;
                        }
                        for (unsigned int tDi = 0; tDi < nTD; tDi++)
                        {
                            if (!needsTH(tDs[tDi], tKpi, tH))
                            {
                                continue;
                            }
                            for (unsigned int pDi = 0; pDi < nPD; pDi++)
                            {
                                const unsigned int di = dIndex(tDi, pDi, g12Index[gpi]);
                                const unsigned int b  = blockIndex(tDi, pDi, gpi, tKpii);
                                ComplexD *S = Sbuf.data() + b*blockSize;
                                ComplexD *R = Rbuf.data() + b*blockSize;
                                for (int id2 = 0; id2 < nDil; id2++)
                                {
                                    chiS = wph[di]*workSlices[id2];
                                    chiR = MTph[di]*workSlices[id2];
                                    for (int id1 = 0; id1 < nDil; id1++)
                                    {
                                        S[id1*nDil + id2] = rankInnerProduct(lightSlices[id1], chiS);
                                        R[id1*nDil + id2] = rankInnerProduct(lightSlices[id1], chiR);
                                    }
                                }
                                blockDone[b] = 1;
                            }
                        }
                    }
                }
                stopTimer("computation D-4quark");
            }
            lightReader.close();

            // sum over the spatial ranks of this timeslice (all of them follow the same control flow)
            startTimer("slice reduction");
            for (unsigned int b = 0; b < nBlocks; b++)
            {
                if (blockDone[b])
                {
                    gridLD->GlobalSumVector(Sbuf.data() + b*blockSize, static_cast<int>(blockSize));
                    gridLD->GlobalSumVector(Rbuf.data() + b*blockSize, static_cast<int>(blockSize));
                }
            }
            stopTimer("slice reduction");
        }

        // serial write: one writer per temporal partition, partitions in turn
        startTimer("serial write I/O");
        for (int iProcT = 0; iProcT < nProcT; iProcT++)
        {
            if ((iProcT == myProcT) && gridLD->IsBoss())
            {
                for (unsigned int b = 0; b < nBlocks; b++)
                {
                    if (!blockDone[b])
                    {
                        continue;
                    }
                    const unsigned int tKpi = tKpis[b % nTKpi];
                    const std::string dataset_name = std::to_string(tKpi)+"-"+std::to_string(tKpi);
                    for (size_t k = 0; k < blockSize; k++)
                    {
                        const ComplexD &z = Sbuf[b*blockSize + k];
                        ioBuf[k] = HADRONS_DISTIL_IO_TYPE(z.real(), z.imag());
                    }
                    matrix_ioS[b].saveBlock(ioBuf.data(), 0, 0, nDil, nDil, std::to_string(tH), dataset_name, nDil);
                    for (size_t k = 0; k < blockSize; k++)
                    {
                        const ComplexD &z = Rbuf[b*blockSize + k];
                        ioBuf[k] = HADRONS_DISTIL_IO_TYPE(z.real(), z.imag());
                    }
                    matrix_ioR[b].saveBlock(ioBuf.data(), 0, 0, nDil, nDil, std::to_string(tH), dataset_name, nDil);
                }
            }
            gridHD->Barrier();
        }
        stopTimer("serial write I/O");
    }
    LOG(Message) << "Perambulator I/O summary on this temporal partition: charm reads = " << charmReads
                 << ", light reads = " << lightReads << ", skipped tH = " << skippedTH << std::endl;
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MDistil_DMeson4QuarkField_hpp_
