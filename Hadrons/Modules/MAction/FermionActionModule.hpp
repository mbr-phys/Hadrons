#ifndef Hadrons_MAction_FermionActionModule_hpp_
#define Hadrons_MAction_FermionActionModule_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/Serialization.hpp>
#include <Hadrons/GaugeBoundary.hpp>

BEGIN_HADRONS_NAMESPACE
BEGIN_MODULE_NAMESPACE(MAction)

template <typename Par>
class FermionActionModule : public Module<Par>
{
public:
    // constructor
    FermionActionModule(const std::string name) : Module<Par>(name){};

protected:
    template <typename ImplPar>
    void parseBoundary(ImplPar &implParams);
    // add the gauge field's boundary metadata to the inputs if it exists, so
    // that it is still alive when checkGaugeBoundary runs in setup
    void addGaugeBoundaryInput(std::vector<std::string> &in);
private:
    template <typename ImplPar>
    void checkGaugeBoundary(const ImplPar &implParams);
};

template <typename Par>
template <typename ImplPar>
void FermionActionModule<Par>::parseBoundary(ImplPar &implParams)
{
    if (!this->par().boundary.empty())
    {
        implParams.boundary_phases = strToVec<Complex>(this->par().boundary);
    }
    if (!this->par().twist.empty())
    {
        implParams.twist_n_2pi_L = strToVec<Real>(this->par().twist);
        if (implParams.twist_n_2pi_L.size() <= 1)
        {
            if (envHasType(std::vector<double>, this->par().twist))
            {
                auto &t = envGet(std::vector<double>, this->par().twist);
                implParams.twist_n_2pi_L = t;
            }
            else if (envHasType(HadronsSerializable, this->par().twist))
            {
                auto &t = envGet(HadronsSerializable, this->par().twist);
                implParams.twist_n_2pi_L = t.template get<std::vector<double>>();
            }
            else
            {
                HADRONS_ERROR(Definition, "cannot interpret twist input '" + this->par().twist + "'");
            }
        }
    }
    LOG(Message) << "Fermion boundary conditions: " << implParams.boundary_phases
                 << std::endl;
    LOG(Message) << "Twists: " << implParams.twist_n_2pi_L
                 << std::endl;
    if (implParams.boundary_phases.size() != this->env().getNd())
    {
        HADRONS_ERROR(Size, "Wrong number of boundary phase");
    }
    if (implParams.twist_n_2pi_L.size() != this->env().getNd())
    {
        HADRONS_ERROR(Size, "Wrong number of twist");
    }
    checkGaugeBoundary(implParams);
}

template <typename Par>
void FermionActionModule<Par>::addGaugeBoundaryInput(std::vector<std::string> &in)
{
    const std::string metadata = GaugeBoundary::metadataName(this->par().gauge);

    if (this->env().hasObject(metadata) &&
        (this->env().getObjectModule(metadata) >= 0))
    {
        in.push_back(metadata);
    }
}

// The temporal fermion boundary must match the gauge field's declared boundary:
// an open gauge field needs a fixed (zero-phase) temporal fermion boundary,
// otherwise the operator hops across the unused link T -> 0 (error). A fixed
// fermion boundary on a periodic or undeclared gauge field is allowed but
// unusual (warning). Single- and double-precision actions are checked alike
// since the gauge precision cast forwards the metadata.
template <typename Par>
template <typename ImplPar>
void FermionActionModule<Par>::checkGaugeBoundary(const ImplPar &implParams)
{
    const std::string gauge    = this->par().gauge;
    const std::string metadata = GaugeBoundary::metadataName(gauge);
    const bool        fixed    =
        (std::abs(implParams.boundary_phases[this->env().getNd() - 1]) == 0.);

    if (this->env().hasCreatedObject(metadata))
    {
        auto &info = envGet(GaugeBoundaryInfo, metadata);

        if ((info.type == GaugeBoundary::open) && !fixed)
        {
            HADRONS_ERROR(Argument, "gauge field '" + gauge + "' has open temporal "
                          "boundaries but the temporal fermion boundary phase "
                          "is not 0");
        }
        else if ((info.type == GaugeBoundary::periodic) && fixed)
        {
            LOG(Warning) << "fixed temporal fermion boundary on gauge field '"
                         << gauge << "', which is declared periodic" << std::endl;
        }
        else if (!GaugeBoundary::isKnown(info) && fixed)
        {
            LOG(Warning) << "fixed temporal fermion boundary on gauge field '"
                         << gauge << "', whose boundary is unknown" << std::endl;
        }
    }
    else if (fixed)
    {
        LOG(Warning) << "fixed temporal fermion boundary on gauge field '"
                     << gauge << "', which has no boundary metadata" << std::endl;
    }
}

END_MODULE_NAMESPACE
END_HADRONS_NAMESPACE

#endif // Hadrons_MAction_FermionActionModule_hpp_