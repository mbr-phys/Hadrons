#ifndef Hadrons_MDistil_SmearedMeson_hpp_
#define Hadrons_MDistil_SmearedMeson_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Serialization.hpp>
#include <Hadrons/NamedTensor.hpp>
#include <Hadrons/Modules/MDistil/Base.hpp>

BEGIN_HADRONS_NAMESPACE
BEGIN_MODULE_NAMESPACE(MDistil)

// Exact-distillation meson two-point functions from native, on-disk
// PerambIndexTensor partitions and pre-existing rho-rho meson-field files.
class SmearedMesonPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(SmearedMesonPar,
                                    std::string,               output,
                                    std::string,               RhoRhoStem,
                                    std::string,               RhoRhoField,
                                    std::string,               perambStemC,
                                    std::string,               perambStemL,
                                    std::string,               noisePol,
                                    std::vector<unsigned int>, tSrcs,
                                    std::string,               gammas,
                                    std::vector<std::string>,  moms);
};

template <typename FImpl>
class TSmearedMeson: public Module<SmearedMesonPar>
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

    TSmearedMeson(const std::string name);
    virtual ~TSmearedMeson(void) {};
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    virtual void setup(void);
    virtual void execute(void);
};

MODULE_REGISTER_TMP(SmearedMeson, TSmearedMeson<FIMPL>, MDistil);

template <typename FImpl>
TSmearedMeson<FImpl>::TSmearedMeson(const std::string name)
: Module<SmearedMesonPar>(name)
{}

template <typename FImpl>
std::vector<std::string> TSmearedMeson<FImpl>::getInput(void)
{
    return {par().noisePol};
}

template <typename FImpl>
std::vector<std::string> TSmearedMeson<FImpl>::getOutput(void)
{
    return {getName()+"_ll", getName()+"_cl", getName()+"_cc"};
}

template <typename FImpl>
void TSmearedMeson<FImpl>::setup(void)
{
    envCreate(HadronsSerializable, getName()+"_ll", 1, 0);
    envCreate(HadronsSerializable, getName()+"_cl", 1, 0);
    envCreate(HadronsSerializable, getName()+"_cc", 1, 0);
}

template <typename FImpl>
void TSmearedMeson<FImpl>::execute(void)
{
    if (par().perambStemL.empty() && par().perambStemC.empty())
    {
        HADRONS_ERROR(Argument, "at least one of perambStemL and perambStemC must be set");
    }
    if (par().RhoRhoStem.empty() || par().RhoRhoField.empty())
    {
        HADRONS_ERROR(Argument, "RhoRhoStem and RhoRhoField must be set");
    }
    if (par().moms.empty())
    {
        HADRONS_ERROR(Argument, "moms must list the source momenta to contract");
    }

    GridCartesian *grid = envGetGrid(FermionField);
    GridCartesian *gridSlice = envGetSliceGrid(FermionField, grid->Nd() - 1);
    const unsigned int nT = env().getDim(Tdir);
    const unsigned int tFirst = grid->LocalStarts()[Tdir];
    const unsigned int tLocal = grid->LocalDimensions()[Tdir];

    LOG(Message) << "gammas = '" << par().gammas << "'" << std::endl;
    std::vector<Gamma::Algebra> gammas = strToVec<Gamma::Algebra>(par().gammas);

    if (gammas.empty()) 
    {
        HADRONS_ERROR(Argument, "gammas must contain at least one sink gamma");
    }

    if (!par().tSrcs.empty() && *std::max_element(par().tSrcs.begin(), par().tSrcs.end()) >= nT)
    {
        HADRONS_ERROR(Range, "all tSrcs must be smaller than nT");
    }

    // The noise policy supplies the native file dimensions.  Exact
    // distillation means every source LapH-spin degree of freedom is a column
    // of tau: one noise vector, nDL == nVec, and full spin dilution.
    auto &dilNoise = envGet(DistillationNoise<FImpl>, par().noisePol);
    const unsigned int nVec = dilNoise.getNl();
    const unsigned int nDL = dilNoise.dilutionSize(DistillationNoise<FImpl>::Index::l);
    const unsigned int nDS = dilNoise.dilutionSize(DistillationNoise<FImpl>::Index::s);
    const unsigned int nNoise = dilNoise.size();
    if (nDL != nVec || nNoise != 1 || nDS != Ns)
    {
        HADRONS_ERROR(Implementation, "SmearedMeson currently requires exact, fully diluted noise");
    }
    const bool hasL = !par().perambStemL.empty();
    const bool hasC = !par().perambStemC.empty();
    const unsigned int nLS = nVec*Ns;
    using Matrix = Eigen::Matrix<Complex, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

    std::map<unsigned int, PerambIndexTensor> perambLCache, perambCCache;
    auto loadPeramb = [&](const std::string &stem, std::map<unsigned int, PerambIndexTensor> &cache,
                          const unsigned int tSrc) -> PerambIndexTensor &
    {
        auto it = cache.find(tSrc);
        if (it == cache.end())
        {
            it = cache.try_emplace(tSrc, nT, nVec, nDL, nNoise, nDS).first;
            const std::string filename = stem + "." + std::to_string(vm().getTrajectory()) +
                                         "/iDT_" + std::to_string(tSrc) + "." +
                                         std::to_string(vm().getTrajectory());
            LOG(Message) << "reading " << filename << ".h5" << std::endl;
            it->second.read(filename);
            if (it->second.MetaData.timeDilutionIndex != static_cast<int>(tSrc))
            {
                HADRONS_ERROR(Io, "perambulator time-dilution index does not match its filename");
            }
        }
        return it->second;
    };

    auto makeTau = [&](const std::string &stem, std::map<unsigned int, PerambIndexTensor> &cache,
                       const unsigned int tSink, const unsigned int tSrc,
                       const Gamma *leftGamma = nullptr)
    {
        auto &p = loadPeramb(stem, cache, tSrc);
        Matrix tau(nLS, nLS);
        for (unsigned int iVec = 0; iVec < nVec; ++iVec)
        for (unsigned int iSpin = 0; iSpin < Ns; ++iSpin)
        for (unsigned int jVec = 0; jVec < nVec; ++jVec)
        for (unsigned int jSpin = 0; jSpin < Ns; ++jSpin)
        {
            auto spin = p.tensor(tSink, iVec, jVec, 0, jSpin)();
            if (leftGamma)
            {
                spin = (*leftGamma)*spin;
            }
            tau(iVec*Ns + iSpin, jVec*Ns + jSpin) =
                spin(iSpin)();
        }
        return tau;
    };

    auto fieldPath = [&](const std::string &gamma, const std::string &mom)
    {
        return par().RhoRhoStem + "rho-rho." + std::to_string(vm().getTrajectory()) + "/" + gamma + "_p" + mom + ".h5";
    };
    // The pre-existing rho-rho meson field is the source operator.  Cache the
    // one requested diagonal block separately for every source time.
    std::map<std::string, ContractionDistilMesonField<ComplexD, ComplexF>> sourceRhoRho;
    auto loadRhoRho = [&](std::map<std::string, ContractionDistilMesonField<ComplexD, ComplexF>> &cache,
                          const std::string &key, const std::string &gamma, const std::string &mom,
                          std::vector<std::vector<int>> &times)
        -> ContractionDistilMesonField<ComplexD, ComplexF> &
    {
        auto it = cache.find(key);
        if (it == cache.end())
        {
            TimerArray timer;
            it = cache.emplace(key, ContractionDistilMesonField<ComplexD, ComplexF>(fieldPath(gamma, mom), nT, timer, times)).first;
        }
        return it->second;
    };

    std::vector<Result> ll, cl, cc;
    const unsigned int nResult = par().tSrcs.size()*par().moms.size()*gammas.size();
    if (hasL) ll.resize(nResult);
    if (hasL && hasC) cl.resize(nResult);
    if (hasC) cc.resize(nResult);

    unsigned int resultIndex = 0;
    for (auto tSrc: par().tSrcs)
    for (const auto &momSrc: par().moms)
    {
        const auto pSrc = strToVec<int>(momSrc);
        if (pSrc.size() != static_cast<size_t>(grid->Nd() - 1))
        {
            HADRONS_ERROR(Size, "each momentum must have Nd-1 components");
        }
        std::string srcMom, snkMom;
        for (unsigned int mu = 0; mu < pSrc.size(); ++mu)
        {
            srcMom += std::to_string(pSrc[mu]) + (mu + 1 == pSrc.size() ? "" : "_");
            snkMom += std::to_string(-pSrc[mu]) + (mu + 1 == pSrc.size() ? "" : "_");
        }
        std::vector<std::vector<int>> srcTVec = {{static_cast<int>(tSrc), static_cast<int>(tSrc)}};
        const std::string srcKey = par().RhoRhoField + "_p" + srcMom + "_t" + std::to_string(tSrc);
        auto &srcField = loadRhoRho(sourceRhoRho, srcKey, par().RhoRhoField, srcMom, srcTVec);
        const Matrix src = srcField(tSrc, tSrc, tSrc);
        if (src.rows() != nLS || src.cols() != nLS)
        {
            HADRONS_ERROR(Size, "rho-rho source matrix is incompatible with perambulator LapH-spin dimension");
        }

        for (const auto gamma: gammas)
        {
            std::stringstream gammaName;
            gammaName << gamma;
            Gamma gam(gamma);
            auto initResult = [&](Result &r)
            {
                r.gammaSrc = par().RhoRhoField;
                r.gammaSnk = gammaName.str();
                r.momSrc = srcMom;
                r.momSnk = snkMom;
                r.tSrc = tSrc;
                r.corr.assign(nT, Complex(0., 0.));
            };
            if (hasL) initResult(ll[resultIndex]);
            if (hasL && hasC) initResult(cl[resultIndex]);
            if (hasC) initResult(cc[resultIndex]);
            for (unsigned int tSink = tFirst; tSink < tFirst + tLocal; ++tSink)
            {
                if (hasL)
                {
                    const Matrix first = makeTau(par().perambStemL, perambLCache, tSink, tSrc);
                    const Matrix second = makeTau(par().perambStemL, perambLCache, tSink, tSrc, &gam);
                    ll[resultIndex].corr[tSink] = (first.adjoint()*second*src).trace();
                }
                if (hasL && hasC)
                {
                    const Matrix first = makeTau(par().perambStemC, perambCCache, tSink, tSrc);
                    const Matrix second = makeTau(par().perambStemL, perambLCache, tSink, tSrc, &gam);
                    cl[resultIndex].corr[tSink] = (first.adjoint()*second*src).trace();
                }
                if (hasC)
                {
                    const Matrix first = makeTau(par().perambStemC, perambCCache, tSink, tSrc);
                    const Matrix second = makeTau(par().perambStemC, perambCCache, tSink, tSrc, &gam);
                    cc[resultIndex].corr[tSink] = (first.adjoint()*second*src).trace();
                }
            }

            ++resultIndex;
        }
    }

    auto reduceAndSave = [&](std::vector<Result> &results, const std::string &suffix, const std::string &tag)
    {
        for (auto &r: results)
        {
            // The values were computed on temporal owner ranks only.
            if (!gridSlice->IsBoss()) std::fill(r.corr.begin(), r.corr.end(), Complex(0., 0.));
            grid->GlobalSumVector(r.corr.data(), static_cast<int>(r.corr.size()));
        }
        saveResult(par().output + suffix, tag, results);
    };
    if (hasL)
    {
        reduceAndSave(ll, "_ll", "llSmearedMeson");
        envGet(HadronsSerializable, getName()+"_ll") = ll;
    }
    if (hasL && hasC)
    {
        reduceAndSave(cl, "_cl", "clSmearedMeson");
        envGet(HadronsSerializable, getName()+"_cl") = cl;
    }
    if (hasC)
    {
        reduceAndSave(cc, "_cc", "ccSmearedMeson");
        envGet(HadronsSerializable, getName()+"_cc") = cc;
    }
}

END_MODULE_NAMESPACE
END_HADRONS_NAMESPACE

#endif
