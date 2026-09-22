/*
 * FermionFlow.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
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
#ifndef Hadrons_MGradientFlow_FermionFlow_hpp_
#define Hadrons_MGradientFlow_FermionFlow_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Serialization.hpp>
#include <Hadrons/TimerArray.hpp>
#include <Hadrons/Modules/MGradientFlow/Utils.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                    Fermion Field Gradient Flow                             *
 *  Flows multiple fields (FermionField and/or PropagatorField) together,    *
 *  sharing gauge field evolution and RK stage computation.                   *
 *                                                                            *
 *  Parameters:                                                               *
 *  - props: list of field names to flow                                      *
 *  - propTypes: optional list of field types ("FermionField" or              *
 *             "PropagatorField"). If empty, defaultType is used for all.    *
 *  - defaultType: default field type when propTypes is empty                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MGradientFlow)

class FermionFlowPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(FermionFlowPar,
                                    std::string, output,
                                    std::vector<std::string>, props,
                                    std::vector<std::string>, propTypes,
                                    std::string, defaultType,
                                    std::vector<std::string>, outProps,
                                    std::string, gauge,
                                    int, bc,
                                    int, steps,
                                    double, step_size,
                                    int, meas_interval);
};

template <typename FImpl, typename GImpl, typename FlowAction>
class TFermionFlow: public Module<FermionFlowPar>
{
public:
    BASIC_TYPE_ALIASES(FImpl,);
    GAUGE_TYPE_ALIASES(GImpl,);
    typedef typename FImpl::FermionField FermionField;
    typedef Evolution<FlowAction, GImpl, FImpl> EvolutionType;
public:
    // constructor
    TFermionFlow(const std::string name);
    // destructor
    virtual ~TFermionFlow(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(FermionFlow,
                    ARG(TFermionFlow<FIMPL, GIMPL, WilsonGaugeAction<GIMPL>>),
                    MGradientFlow);

/******************************************************************************
 *                     TFermionFlow implementation                            *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl, typename GImpl, typename FlowAction>
TFermionFlow<FImpl,GImpl,FlowAction>::TFermionFlow(const std::string name)
: Module<FermionFlowPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl, typename GImpl, typename FlowAction>
std::vector<std::string> TFermionFlow<FImpl,GImpl,FlowAction>::getInput(void)
{
    std::vector<std::string> in = {par().gauge}; 
    for (std::string q : par().props) {
        in.push_back(q);
    }
    
    return in;
}

template <typename FImpl, typename GImpl, typename FlowAction>
std::vector<std::string> TFermionFlow<FImpl,GImpl,FlowAction>::getOutput(void)
{
    std::vector<std::string> out = {getName(), getName()+"_U"};

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
void TFermionFlow<FImpl,GImpl,FlowAction>::setup(void)
{
    envCreateLat(GaugeField, getName()+"_U");

    // Validate props and propTypes
    if (!par().propTypes.empty() && (par().props.size() != par().propTypes.size())) {
        HADRONS_ERROR(Argument, "propTypes must be empty or have the same size as props");
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
                for (std::string q : par().outProps) {
                    // Infer type from input props
                    auto it = std::find(par().props.begin(), par().props.end(), q);
                    if (it != par().props.end()) {
                        size_t idx = std::distance(par().props.begin(), it);
                        std::string type = fieldTypes[idx];
                        if (type == "FermionField") {
                            envCreateLat(FermionField, q);
                        } else if (type == "PropagatorField") {
                            envCreateLat(PropagatorField, q);
                        } else {
                            HADRONS_ERROR(Argument, "Unknown field type: " + type + " for field " + q);
                        }
                    } else {
                        HADRONS_ERROR(Argument, "outProp " + q + " not found in props");
                    }
                }
            }
        }
    }
    envCreate(HadronsSerializable, getName(), 1, 0);
    envTmp(EvolutionType, "evolve", 1, envGetGrid(GaugeField), 3.0, par().step_size, 
        -1.0, par().step_size);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl, typename GImpl, typename FlowAction>
void TFermionFlow<FImpl,GImpl,FlowAction>::execute(void)
{
    // action
    FlowAction SG = FlowAction(3.0);

    std::string type = SG.action_name();
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

    LOG(Message) << "Setting up " << type << " Field Flow on '" << par().gauge << "' Gauge Field and "
                 << props << "with " << par().props.size() << " field" 
                 << ((par().props.size() > 1) ? "s" : "")
                 << " with ppp" << ((par().bc < 0) ? "a" : "p") << " boundary conditions and "
                 << par().steps << " step" << ((par().steps > 1) ? "s." : ".") << std::endl;

    if ((par().outProps.size() != par().props.size()) && !par().outProps.empty()) {
        HADRONS_ERROR(Argument, "outProps should either be empty or be the same size as props");
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

    // set boundary conditions for gauge field
    std::vector<int> bc = {1,1,1};
    if (par().bc < 0) bc.push_back(-1);
    else bc.push_back(1);

    auto &out     = envGet(HadronsSerializable, getName());
    auto &Uresult = out.template hold<GaugeResult>();
    envGetTmp(EvolutionType, evolve);

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
        std::vector<GaugeField> &Wi = evolve.evolve_gaugeFF(Uwf,bc);
        stopTimer("gauge field flow time "+ftt.str());

        // measure gauge observables
        evolve.gauge_status(Uwf,Uresult,flowt);

        // Flow ALL FermionFields with shared gauge stages
        for (const auto& q : fermionFieldNames) {
            FermionField &qjwf = *env().template getObject<FermionField>(getName()+"_tmp_"+q+"_wf");
            startTimer("FermionField "+q+" flow time "+ftt.str());
            evolve.laplace_flow(Wi[0],Wi[1],Wi[2],qjwf);
            stopTimer("FermionField "+q+" flow time "+ftt.str());
        }
        
        // Flow ALL PropagatorFields with shared gauge stages
        for (const auto& q : propFieldNames) {
            PropagatorField &qjwf = *env().template getObject<PropagatorField>(getName()+"_tmp_"+q+"_wf");
            startTimer("PropagatorField "+q+" flow time "+ftt.str());
            evolve.laplace_flow(Wi[0],Wi[1],Wi[2],qjwf);
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

#endif // Hadrons_MGradientFlow_FermionFlow_hpp_
