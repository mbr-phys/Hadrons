/*
 * AdjointFermionFlow_OBC.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2026
 *
 * Author: Matthew Black <matthewkblack@protonmail.com>
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
#ifndef Hadrons_MGradientFlow_AdjointFermionFlow_OBC_hpp_
#define Hadrons_MGradientFlow_AdjointFermionFlow_OBC_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Serialization.hpp>
#include <Hadrons/GaugeBoundary.hpp>
#include <Hadrons/Modules/MGradientFlow/Utils_OBC.hpp>

BEGIN_HADRONS_NAMESPACE

/*
 * Open-boundary adjoint fermion flow module
 * ---------------------------
 * Implements backward fermion flow using reverse RK stepping.
 * 
 * Takes stochastic noise sources eta defined at flow time tau and flows them
 * backward to s=0, producing back-flowed sources xi(tau;0) that can be used
 * as RHS for ordinary Dirac solves.
 * 
 * Gauge field must be provided as input (pre-computed via forward flow).
 * Supports step-by-step workflow (one step per module instance).
 * 
 * Parameters:
 * - gauge: gauge field U_s at earlier flow time (for computing RK stages)
 * - sources: input fields xi at s+eps (to be flowed backward)
 * - outSources: output fields xi at s (after backward step)
 * - steps: number of steps (typically 1 for step-by-step workflow)
 * - step_size: eps (flow step size)
 * - bc: fermion temporal BC (+1 periodic, -1 antiperiodic, 0 fixed)
 * - cG: open gauge-boundary coefficient (used only when bc == 0)
 */

/******************************************************************************
 *                         TAdjointFermionFlow_OBC                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MGradientFlow)

class AdjointFermionFlow_OBCPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(AdjointFermionFlow_OBCPar,
                                    std::string, gauge,
                                    std::string, stage1,
                                    std::string, stage2,
                                    std::vector<std::string>, sources,
                                    std::vector<std::string>, sourceTypes,
                                    std::string, defaultType,
                                    std::vector<std::string>, outSources,
                                    int, bc,
                                    double, cG,
                                    int, steps,
                                    double, step_size);
};

template <typename FImpl, typename GImpl, typename FlowAction>
class TAdjointFermionFlow_OBC: public Module<AdjointFermionFlow_OBCPar>
{
public:
    BASIC_TYPE_ALIASES(FImpl,);
    GAUGE_TYPE_ALIASES(GImpl,);
    typedef typename FImpl::FermionField FermionField;
    typedef OBC::Evolution<FlowAction, GImpl, FImpl> EvolutionType;
public:
    // constructor
    TAdjointFermionFlow_OBC(const std::string name);
    // destructor
    virtual ~TAdjointFermionFlow_OBC(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(AdjointFermionFlow_OBC,ARG(TAdjointFermionFlow_OBC<FIMPL, GIMPL, OBC::WilsonAction<GIMPL>>),MGradientFlow);

/******************************************************************************
 *                    TAdjointFermionFlow_OBC implementation                      *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl, typename GImpl, typename FlowAction>
TAdjointFermionFlow_OBC<FImpl,GImpl,FlowAction>::TAdjointFermionFlow_OBC(const std::string name)
: Module<AdjointFermionFlow_OBCPar>(name)
{}

// dependencies ////////////////////////////////////////////////////////////////
template <typename FImpl, typename GImpl, typename FlowAction>
std::vector<std::string> TAdjointFermionFlow_OBC<FImpl,GImpl,FlowAction>::getInput(void)
{
    std::vector<std::string> in = {par().gauge};
    auto addMetadataDependency = [this, &in](const std::string &field) {
        const std::string metadata = GaugeBoundary::metadataName(field);
        if (env().hasObject(metadata) &&
            (env().getObjectModule(metadata) >= 0)) {
            in.push_back(metadata);
        }
    };
    addMetadataDependency(par().gauge);
    if (par().stage1.empty() != par().stage2.empty()) {
        HADRONS_ERROR(Argument, "stage1 and stage2 must either both be set or both be empty");
    }
    if (!par().stage1.empty()) {
        in.push_back(par().stage1);
        in.push_back(par().stage2);
        addMetadataDependency(par().stage1);
        addMetadataDependency(par().stage2);
    }
    for (const auto& src : par().sources) {
        in.push_back(src);
    }
    return in;
}

template <typename FImpl, typename GImpl, typename FlowAction>
std::vector<std::string> TAdjointFermionFlow_OBC<FImpl,GImpl,FlowAction>::getOutput(void)
{
    std::vector<std::string> out = par().outSources;
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl, typename GImpl, typename FlowAction>
void TAdjointFermionFlow_OBC<FImpl,GImpl,FlowAction>::setup(void)
{
    if (par().sources.empty() || par().outSources.empty()) {
        HADRONS_ERROR(Argument, "sources and outSources must both be non-empty");
    }
    if (par().sources.size() != par().outSources.size()) {
        HADRONS_ERROR(Argument, "sources and outSources must have the same size");
    }
    if (par().steps != 1) {
        HADRONS_ERROR(Argument, "AdjointFermionFlow currently supports exactly one step; provide the earlier gauge field for each module instance");
    }
    if (!std::isfinite(par().step_size) || par().step_size <= 0.0) {
        HADRONS_ERROR(Argument, "step_size must be finite and positive");
    }
    if (par().bc != -1 && par().bc != 0 && par().bc != 1) {
        HADRONS_ERROR(Argument, "bc must be +1 (periodic), -1 (antiperiodic), or 0 (open)");
    }
    if (!std::isfinite(par().cG) || par().cG < 0.0) {
        HADRONS_ERROR(Argument,
                      "cG must be explicitly set to a finite non-negative value");
    }
    if (par().stage1.empty() != par().stage2.empty()) {
        HADRONS_ERROR(Argument, "stage1 and stage2 must either both be set or both be empty");
    }

    std::vector<std::string> fieldTypes;
    if (par().sourceTypes.empty()) {
        if (par().defaultType.empty()) {
            HADRONS_ERROR(Argument, "defaultType must be specified when sourceTypes is empty");
        }
        fieldTypes.resize(par().sources.size(), par().defaultType);
    } else {
        if (par().sourceTypes.size() != par().sources.size()) {
            HADRONS_ERROR(Argument, "sourceTypes must be empty or have the same size as sources");
        }
        fieldTypes = par().sourceTypes;
    }

    // Create output fields
    for (size_t i = 0; i < par().outSources.size(); ++i) {
        if (fieldTypes[i] == "FermionField") {
            envCreateLat(FermionField, par().outSources[i]);
        } else if (fieldTypes[i] == "PropagatorField") {
            envCreateLat(PropagatorField, par().outSources[i]);
        } else {
            HADRONS_ERROR(Argument, "Unknown source field type: " + fieldTypes[i]);
        }
    }

    const std::string gaugeBc = GaugeBoundary::gaugeTypeFromFermionBc(par().bc);
    envTmp(EvolutionType, "evolve", 1, envGetGrid(GaugeField), 3.0,
           par().step_size, par().step_size, 0.0,
           gaugeBc, par().cG, this);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl, typename GImpl, typename FlowAction>
void TAdjointFermionFlow_OBC<FImpl,GImpl,FlowAction>::execute(void)
{
    envGetTmp(EvolutionType, evolve);

    const std::string expectedGaugeBc =
        GaugeBoundary::gaugeTypeFromFermionBc(par().bc);
    auto validateBoundaryMetadata =
        [this, &expectedGaugeBc](const std::string &field,
                                const std::string &role)
            -> GaugeBoundaryInfo * {
            const std::string metadata = GaugeBoundary::metadataName(field);
            if (!env().hasCreatedObject(metadata)) {
                LOG(Warning) << "No gauge boundary metadata found for " << role
                             << " '" << field
                             << "'; proceeding with module parameters"
                             << std::endl;
                return nullptr;
            }

            auto &info = envGet(GaugeBoundaryInfo, metadata);
            const std::string mismatch = GaugeBoundary::compatibilityError(
                info, expectedGaugeBc, par().cG);
            if (!mismatch.empty()) {
                HADRONS_ERROR(Argument,
                              "gauge boundary metadata mismatch for " + role +
                              " '" + field + "': " + mismatch);
            }
            if (!GaugeBoundary::isKnown(info)) {
                LOG(Warning) << "Gauge boundary metadata for " << role << " '"
                             << field
                             << "' is unknown; proceeding with module parameters"
                             << std::endl;
            }
            return &info;
        };

    GaugeBoundaryInfo *gaugeInfo =
        validateBoundaryMetadata(par().gauge, "gauge field");
    if (!par().stage1.empty()) {
        GaugeBoundaryInfo *stage1Info =
            validateBoundaryMetadata(par().stage1, "RK stage W1");
        GaugeBoundaryInfo *stage2Info =
            validateBoundaryMetadata(par().stage2, "RK stage W2");

        if ((stage1Info != nullptr) && (stage2Info != nullptr) &&
            stage1Info->hasFlowInfo && stage2Info->hasFlowInfo) {
            const std::string expectedAction = evolve.SG.action_name();
            auto close = [](const double a, const double b) {
                const double scale = 1.0 + std::max(std::abs(a), std::abs(b));
                return std::abs(a - b) <= 1.0e-12 * scale;
            };

            if ((stage1Info->flowAction != expectedAction) ||
                (stage2Info->flowAction != expectedAction)) {
                HADRONS_ERROR(Argument,
                              "saved RK stages use flow action '" +
                              stage1Info->flowAction + "'/'" +
                              stage2Info->flowAction +
                              "' but adjoint reconstruction expects '" +
                              expectedAction + "'");
            }
            if (!close(stage1Info->stepSize, par().step_size) ||
                !close(stage2Info->stepSize, par().step_size)) {
                HADRONS_ERROR(Argument,
                              "saved RK stage step size does not match "
                              "AdjointFermionFlow_OBC step_size");
            }
            if ((stage1Info->stage != "W1") ||
                (stage2Info->stage != "W2")) {
                HADRONS_ERROR(Argument,
                              "saved RK stage metadata does not identify W1/W2");
            }
            if (!close(stage1Info->flowTime, stage2Info->flowTime)) {
                HADRONS_ERROR(Argument,
                              "saved W1 and W2 stages have different start times");
            }
            if ((gaugeInfo != nullptr) && gaugeInfo->hasFlowInfo) {
                if ((gaugeInfo->flowAction != expectedAction) ||
                    !close(gaugeInfo->stepSize, par().step_size) ||
                    !close(gaugeInfo->flowTime, stage1Info->flowTime)) {
                    HADRONS_ERROR(Argument,
                                  "gauge field and saved RK stages have "
                                  "inconsistent flow provenance");
                }
            }
        } else {
            LOG(Warning) << "Saved RK stages do not both provide flow "
                         << "provenance metadata; only boundary compatibility "
                         << "was checked." << std::endl;
        }
    }

    LOG(Message) << "Starting adjoint fermion flow for '" << getName() << "'" << std::endl;
    LOG(Message) << "  gauge field: " << par().gauge << std::endl;
    if (!par().stage1.empty()) {
        LOG(Message) << "  saved RK stages: " << par().stage1 << " "
                     << par().stage2 << std::endl;
    }
    const std::string gaugeBc = GaugeBoundary::gaugeTypeFromFermionBc(par().bc);
    LOG(Message) << "  fermion temporal BC: " << par().bc << std::endl;
    LOG(Message) << "  gauge temporal BC: " << gaugeBc;
    if (par().bc == 0) {
        LOG(Message) << " (cG=" << par().cG << ")";
    }
    LOG(Message) << std::endl;
    LOG(Message) << "  steps: " << par().steps << ", step_size: " << par().step_size << std::endl;
    LOG(Message) << "  sources: ";
    for (const auto& src : par().sources) LOG(Message) << src << " ";
    LOG(Message) << std::endl;
    LOG(Message) << "  outSources: ";
    for (const auto& out : par().outSources) LOG(Message) << out << " ";
    LOG(Message) << std::endl;
    std::cout << evolve.SG.LogParameters();

    // Get gauge field
    auto &U = envGet(GaugeField, par().gauge);
    
    std::vector<std::string> fieldTypes;
    if (par().sourceTypes.empty()) {
        fieldTypes.resize(par().sources.size(), par().defaultType);
    } else {
        fieldTypes = par().sourceTypes;
    }

    // Get input and output fields, retaining the source/output pairing.
    std::vector<FermionField*> fermionIn(par().sources.size(), nullptr);
    std::vector<FermionField*> fermionOut(par().outSources.size(), nullptr);
    std::vector<PropagatorField*> propagatorIn(par().sources.size(), nullptr);
    std::vector<PropagatorField*> propagatorOut(par().outSources.size(), nullptr);
    for (size_t i = 0; i < par().sources.size(); i++) {
        if (fieldTypes[i] == "FermionField") {
            fermionIn[i] = &envGet(FermionField, par().sources[i]);
            fermionOut[i] = &envGet(FermionField, par().outSources[i]);
        } else if (fieldTypes[i] == "PropagatorField") {
            propagatorIn[i] = &envGet(PropagatorField, par().sources[i]);
            propagatorOut[i] = &envGet(PropagatorField, par().outSources[i]);
        } else {
            HADRONS_ERROR(Argument, "Unknown source field type: " + fieldTypes[i]);
        }
    }
    
    std::vector<int> bc = {1, 1, 1, par().bc};
    
    // Apply adjoint flow for each step.
    for (int step = 0; step < par().steps; step++) {
        std::stringstream step_ss;
        step_ss << std::fixed << std::setprecision(2) << ((step + 1) * par().step_size);
        LOG(Message) << "Adjoint flow step " << (step + 1) << "/" << par().steps 
                     << " (tau=" << step_ss.str() << ")" << std::endl;
        
        auto applyAdjoint = [&](GaugeField &W0, GaugeField &W1, GaugeField &W2) {
            for (size_t i = 0; i < par().sources.size(); i++) {
                if (fieldTypes[i] == "FermionField") {
                    *fermionOut[i] = *fermionIn[i];
                    evolve.adjoint_laplace_flow(W0, W1, W2, *fermionOut[i], par().bc);
                } else {
                    *propagatorOut[i] = *propagatorIn[i];
                    evolve.adjoint_laplace_flow(W0, W1, W2, *propagatorOut[i], par().bc);
                }
            }
        };

        if (par().stage1.empty()) {
            GaugeField U_copy = U;
            std::vector<GaugeField> &Wi = evolve.evolve_gaugeFF(U_copy, bc);
            applyAdjoint(Wi[0], Wi[1], Wi[2]);
        } else {
            // The saved gauge fields are raw gauge-flow stages. Apply the
            // fermion temporal boundary condition to private copies.
            GaugeField W0 = U;
            GaugeField W1 = envGet(GaugeField, par().stage1);
            GaugeField W2 = envGet(GaugeField, par().stage2);
            evolve.gauge_apply_boundary(W0, bc);
            evolve.gauge_apply_boundary(W1, bc);
            evolve.gauge_apply_boundary(W2, bc);
            applyAdjoint(W0, W1, W2);
        }
    }
    
    LOG(Message) << "Adjoint fermion flow complete for '" << getName() << "'" << std::endl;
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MGradientFlow_AdjointFermionFlow_OBC_hpp_
