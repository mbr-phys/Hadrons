/*
 * FermionFlow_OBC.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2022
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Matthew Black    <matthewkblack@protonmail.com>
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
#ifndef Hadrons_MGradientFlow_FermionFlow_OBC_hpp_
#define Hadrons_MGradientFlow_FermionFlow_OBC_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Serialization.hpp>
#include <Hadrons/GaugeBoundary.hpp>
#include <Hadrons/TimerArray.hpp>
#include <Hadrons/Modules/MGradientFlow/Utils_OBC.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                    Open-Boundary Fermion Field Gradient Flow                             *
 *  Flows multiple fields (FermionField and/or PropagatorField) together,    *
 *  sharing gauge field evolution and RK stage computation.                   *
 *                                                                            *
 *  Parameters:                                                               *
 *  - props: list of field names to flow                                      *
 *  - propTypes: optional list of field types ("FermionField" or              *
 *             "PropagatorField"). If empty, defaultType is used for all.    *
 *  - defaultType: default field type when propTypes is empty                 *
 *  - bc: fermion temporal BC: +1 periodic, -1 antiperiodic, 0 fixed *
 *  - cG: open gauge-boundary coefficient (used only when bc == 0)       *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MGradientFlow)

class FermionFlow_OBCPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(FermionFlow_OBCPar,
                                    std::string, output,
                                    std::vector<std::string>, props,
                                    std::vector<std::string>, propTypes,
                                    std::string, defaultType,
                                    std::vector<std::string>, outProps,
                                    std::string, gauge,
                                    int, bc,
                                    double, cG,
                                    int, steps,
                                    double, step_size,
                                    int, meas_interval);
};

template <typename FImpl, typename GImpl, typename FlowAction>
class TFermionFlow_OBC: public Module<FermionFlow_OBCPar>
{
public:
    BASIC_TYPE_ALIASES(FImpl,);
    GAUGE_TYPE_ALIASES(GImpl,);
    typedef typename FImpl::FermionField FermionField;
    typedef OBC::Evolution<FlowAction, GImpl, FImpl> EvolutionType;
public:
    // constructor
    TFermionFlow_OBC(const std::string name);
    // destructor
    virtual ~TFermionFlow_OBC(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(FermionFlow_OBC,ARG(TFermionFlow_OBC<FIMPL, GIMPL, OBC::WilsonAction<GIMPL>>),MGradientFlow);

/******************************************************************************
 *                     TFermionFlow_OBC implementation                            *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl, typename GImpl, typename FlowAction>
TFermionFlow_OBC<FImpl,GImpl,FlowAction>::TFermionFlow_OBC(const std::string name)
: Module<FermionFlow_OBCPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl, typename GImpl, typename FlowAction>
std::vector<std::string> TFermionFlow_OBC<FImpl,GImpl,FlowAction>::getInput(void)
{
    std::vector<std::string> in = {par().gauge};
    const std::string gaugeMetadata =
        GaugeBoundary::metadataName(par().gauge);
    if (env().hasObject(gaugeMetadata) &&
        (env().getObjectModule(gaugeMetadata) >= 0)) {
        in.push_back(gaugeMetadata);
    }
    for (std::string q : par().props) {
        in.push_back(q);
    }

    return in;
}

template <typename FImpl, typename GImpl, typename FlowAction>
std::vector<std::string> TFermionFlow_OBC<FImpl,GImpl,FlowAction>::getOutput(void)
{
    if (par().steps < 0) {
        HADRONS_ERROR(Argument, "steps must be non-negative");
    }
    if (par().meas_interval <= 0) {
        HADRONS_ERROR(Argument, "meas_interval must be positive");
    }

    std::vector<std::string> out = {
        getName(), getName()+"_U",
        GaugeBoundary::metadataName(getName()+"_U")
    };

    // output flowed fields at measurement intervals
    for (int i = 1; i <= par().steps; i++) 
    {
        if ((i % par().meas_interval == 0) || (i == par().steps)) {
            double ft = par().step_size * i;
            std::stringstream ftt; ftt << std::fixed << std::setprecision(2) << ft;
            if (par().outProps.empty()) {
                for (std::string q : par().props) {
                    out.push_back(q+"_t"+ftt.str());
                }
            } else {
                for (std::string q : par().outProps) {
                    out.push_back(q);
                }
            }
        }
    }

    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl, typename GImpl, typename FlowAction>
void TFermionFlow_OBC<FImpl,GImpl,FlowAction>::setup(void)
{
    if (par().bc != -1 && par().bc != 0 && par().bc != 1) {
        HADRONS_ERROR(Argument, "bc must be +1 (periodic), -1 (antiperiodic), or 0 (open)");
    }
    if (par().steps < 0) {
        HADRONS_ERROR(Argument, "steps must be non-negative");
    }
    if (par().meas_interval <= 0) {
        HADRONS_ERROR(Argument, "meas_interval must be positive");
    }
    if (!std::isfinite(par().step_size) || par().step_size <= 0.0) {
        HADRONS_ERROR(Argument, "step_size must be finite and positive");
    }
    if (!std::isfinite(par().cG) || par().cG < 0.0) {
        HADRONS_ERROR(Argument,
                      "cG must be explicitly set to a finite non-negative value");
    }

    envCreateLat(GaugeField, getName()+"_U");
    const std::string effectiveGaugeBc =
        GaugeBoundary::gaugeTypeFromFermionBc(par().bc);
    const double effectiveCG = (par().bc == 0) ? par().cG : 1.0;
    envCreate(GaugeBoundaryInfo,
              GaugeBoundary::metadataName(getName()+"_U"), 1,
              effectiveGaugeBc, Nd - 1, effectiveCG, true,
              true, "WilsonAction", par().step_size,
              par().steps * par().step_size, "field");

    // Validate props and propTypes
    if (!par().propTypes.empty() && (par().props.size() != par().propTypes.size())) {
        HADRONS_ERROR(Argument, "propTypes must be empty or have the same size as props");
    }
    if (!par().outProps.empty() && (par().props.size() != par().outProps.size())) {
        HADRONS_ERROR(Argument, "outProps must be empty or have the same size as props");
    }

    // Determine field types for each prop
    std::vector<std::string> fieldTypes;
    if (par().propTypes.empty()) {
        // Use defaultType for all fields
        if (par().defaultType.empty()) {
            HADRONS_ERROR(Argument, "defaultType must be specified when propTypes is empty");
        }
        fieldTypes.resize(par().props.size(), par().defaultType);
    } else {
        fieldTypes = par().propTypes;
    }

    // Validate field types and create temporaries
    for (size_t i = 0; i < par().props.size(); i++) {
        std::string q = par().props[i];
        std::string type = fieldTypes[i];
        if (type == "FermionField") {
            envTmpLat(FermionField, q+"_wf");
        } else if (type == "PropagatorField") {
            envTmpLat(PropagatorField, q+"_wf");
        } else {
            HADRONS_ERROR(Argument, "Unknown field type: " + type + " for field " + q);
        }
    }

    // create output fields
    for (int i = 1; i <= par().steps; i++) 
    {
        if (( i % par().meas_interval == 0) || (i == par().steps)) {
            double ft = par().step_size * i;
            std::stringstream ftt; ftt << std::fixed << std::setprecision(2) << ft;
            if (par().outProps.empty()) {
                for (size_t j = 0; j < par().props.size(); j++) {
                    std::string q = par().props[j];
                    std::string type = fieldTypes[j];
                    if (type == "FermionField") {
                        envCreateLat(FermionField, q+"_t"+ftt.str());
                    } else if (type == "PropagatorField") {
                        envCreateLat(PropagatorField, q+"_t"+ftt.str());
                    } else {
                        HADRONS_ERROR(Argument, "Unknown field type: " + type + " for field " + q);
                    }
                }
            } else {
                for (size_t j = 0; j < par().outProps.size(); j++) {
                    std::string q = par().outProps[j];
                    std::string type = fieldTypes[j];
                    if (type == "FermionField") {
                        envCreateLat(FermionField, q);
                    } else if (type == "PropagatorField") {
                        envCreateLat(PropagatorField, q);
                    } else {
                        HADRONS_ERROR(Argument, "Unknown field type: " + type + " for field " + q);
                    }
                }
            }
        }
    }
    envCreate(HadronsSerializable, getName(), 1, 0);
    const std::string gaugeBc = GaugeBoundary::gaugeTypeFromFermionBc(par().bc);
    envTmp(EvolutionType, "evolve", 1, envGetGrid(GaugeField), 3.0,
           par().step_size, -1.0, par().step_size, gaugeBc, par().cG, this);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl, typename GImpl, typename FlowAction>
void TFermionFlow_OBC<FImpl,GImpl,FlowAction>::execute(void)
{
    envGetTmp(EvolutionType, evolve);

    const std::string expectedGaugeBc =
        GaugeBoundary::gaugeTypeFromFermionBc(par().bc);
    const std::string inputMetadata =
        GaugeBoundary::metadataName(par().gauge);
    if (env().hasCreatedObject(inputMetadata)) {
        auto &info = envGet(GaugeBoundaryInfo, inputMetadata);
        const std::string mismatch = GaugeBoundary::compatibilityError(
            info, expectedGaugeBc, Nd - 1, par().cG);
        if (!mismatch.empty()) {
            HADRONS_ERROR(Argument,
                          "gauge boundary metadata mismatch for '" +
                          par().gauge + "': " + mismatch);
        }
        if (!GaugeBoundary::isKnown(info)) {
            LOG(Warning) << "Gauge boundary metadata for '" << par().gauge
                         << "' is unknown; proceeding with module parameters"
                         << std::endl;
        }
    } else {
        LOG(Warning) << "No gauge boundary metadata found for '" << par().gauge
                     << "'; proceeding with module parameters" << std::endl;
    }

    std::string type = evolve.SG.action_name();
    std::string ga = "GaugeAction";
    std::string::size_type i = type.find(ga);
    if (i != std::string::npos) {
        type.erase(i, ga.length());
    }

    std::string props = "";
    for (std::string q : par().props) props += q + " ";
    
    // Determine field types for logging
    std::vector<std::string> fieldTypes;
    if (par().propTypes.empty()) {
        fieldTypes.resize(par().props.size(), par().defaultType);
    } else {
        fieldTypes = par().propTypes;
    }

    const char *bcName = (par().bc > 0) ? "periodic" :
                         ((par().bc < 0) ? "antiperiodic" : "open");
    const std::string gaugeBc = GaugeBoundary::gaugeTypeFromFermionBc(par().bc);
    LOG(Message) << "Setting up " << type << " Field Flow on '" << par().gauge << "' Gauge Field and "
                 << props << "with " << par().props.size() << " field"
                 << ((par().props.size() > 1) ? "s" : "")
                 << " with " << bcName << " temporal fermion boundary conditions and "
                 << gaugeBc << " temporal gauge boundary conditions"
                 << ((par().bc == 0) ? " (cG=" + std::to_string(par().cG) + ")" : "")
                 << ", using " << par().steps << " step"
                 << ((par().steps > 1) ? "s." : ".") << std::endl;
    std::cout << evolve.SG.LogParameters();
    if (par().bc == 0) {
        LOG(Warning) << "FermionFlow_OBC gauge results retain the existing "
                     << "naive gauge diagnostics; the recorded action uses "
                     << "the OBC action and the temporal Polyakov loop is not "
                     << "a physical closed-loop OBC observable." << std::endl;
    }

    // Validate field types and separate into parallel vectors
    std::vector<std::string> fermionFieldNames;
    std::vector<std::string> propFieldNames;
    
    for (size_t i = 0; i < par().props.size(); i++) {
        std::string q = par().props[i];
        std::string type = fieldTypes[i];
        if (type == "FermionField") {
            fermionFieldNames.push_back(q);
        } else if (type == "PropagatorField") {
            propFieldNames.push_back(q);
        } else {
            HADRONS_ERROR(Argument, "Unknown field type: " + type);
        }
    }

    // Spatial boundaries remain periodic.  For bc == 0 the evolution utility
    // additionally projects the fermion field at both temporal endpoints.
    std::vector<int> gaugeStageBc = {1, 1, 1, par().bc};

    auto &out     = envGet(HadronsSerializable, getName());
    auto &Uresult = out.template hold<OBC::GaugeResult>();

    auto &U   = envGet(GaugeField, par().gauge);
    auto &Uwf = envGet(GaugeField, getName()+"_U");
    Uwf = U;

    // Initialize all flowed fields (both FermionField and PropagatorField)
    for (const auto& q : fermionFieldNames) {
        auto &qj = envGet(FermionField, q);
        FermionField &qjwf = *env().template getObject<FermionField>(getName()+"_tmp_"+q+"_wf");
        qjwf = qj;
    }
    for (const auto& q : propFieldNames) {
        auto &qj = envGet(PropagatorField, q);
        PropagatorField &qjwf = *env().template getObject<PropagatorField>(getName()+"_tmp_"+q+"_wf");
        qjwf = qj;
    }
    
    // apply flow equations
    double flowt = 0.0;
    evolve.gauge_status(Uwf,Uresult,flowt);
    
    for (unsigned int step = 1; step <= par().steps; step++) {
        flowt += evolve.epsilon;
        std::stringstream ftt; ftt << std::fixed << std::setprecision(2) << flowt;

        // evolve gauge field ONCE for all fields
        startTimer("gauge field flow time "+ftt.str());
        std::vector<GaugeField> &Wi = evolve.evolve_gaugeFF(Uwf,gaugeStageBc);
        stopTimer("gauge field flow time "+ftt.str());

        // measure gauge observables
        evolve.gauge_status(Uwf,Uresult,flowt);

        // Flow ALL FermionFields with shared gauge stages
        for (const auto& q : fermionFieldNames) {
            FermionField &qjwf = *env().template getObject<FermionField>(getName()+"_tmp_"+q+"_wf");
            startTimer("FermionField "+q+" flow time "+ftt.str());
            evolve.laplace_flow(Wi[0],Wi[1],Wi[2],qjwf, par().bc);
            stopTimer("FermionField "+q+" flow time "+ftt.str());
        }
        
        // Flow ALL PropagatorFields with shared gauge stages
        for (const auto& q : propFieldNames) {
            PropagatorField &qjwf = *env().template getObject<PropagatorField>(getName()+"_tmp_"+q+"_wf");
            startTimer("PropagatorField "+q+" flow time "+ftt.str());
            evolve.laplace_flow(Wi[0],Wi[1],Wi[2],qjwf, par().bc);
            stopTimer("PropagatorField "+q+" flow time "+ftt.str());
        }

        // Save outputs at measurement intervals
        if (( step % par().meas_interval == 0) || (step == par().steps)) {
            std::string suffix = "_t"+ftt.str();
            for (size_t j = 0; j < fermionFieldNames.size(); j++) {
                std::string q = fermionFieldNames[j];
                // Find global index in original props list
                auto it = std::find(par().props.begin(), par().props.end(), q);
                size_t globalIdx = std::distance(par().props.begin(), it);
                std::string qo = par().outProps.empty() ? q + suffix : par().outProps[globalIdx];
                auto &qji = envGet(FermionField, qo);
                FermionField &qjwf = *env().template getObject<FermionField>(getName()+"_tmp_"+q+"_wf");
                qji = qjwf;
            }
            for (size_t j = 0; j < propFieldNames.size(); j++) {
                std::string q = propFieldNames[j];
                // Find global index in original props list
                auto it = std::find(par().props.begin(), par().props.end(), q);
                size_t globalIdx = std::distance(par().props.begin(), it);
                std::string qo = par().outProps.empty() ? q + suffix : par().outProps[globalIdx];
                auto &qji = envGet(PropagatorField, qo);
                PropagatorField &qjwf = *env().template getObject<PropagatorField>(getName()+"_tmp_"+q+"_wf");
                qji = qjwf;
            }
        }
    }
    saveResult(par().output,"gauge_obs",Uresult);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MGradientFlow_FermionFlow_OBC_hpp_
