/*
 * GaugeFlow_OBC.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2026
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
#ifndef Hadrons_MGradientFlow_GaugeFlow_OBC_hpp_
#define Hadrons_MGradientFlow_GaugeFlow_OBC_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/GaugeBoundary.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Serialization.hpp>
#include <Hadrons/Modules/MGradientFlow/Utils_OBC.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *             Gauge-field gradient flow with explicit gauge BC              *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MGradientFlow)

class GaugeFlow_OBCPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(GaugeFlow_OBCPar,
                                    std::string, output,
                                    std::string, gauge,
                                    int, steps,
                                    double, step_size,
                                    int, meas_interval,
                                    bool, save_history,
                                    bool, save_rk_stages,
                                    std::string, gaugeBC,
                                    double, cG,
                                    std::string, maxTau,
                                    std::string, c_plaq,
                                    std::string, c_rect);
};

template <typename GImpl, typename FlowAction>
class TGaugeFlow_OBC: public Module<GaugeFlow_OBCPar>
{
public:
    GAUGE_TYPE_ALIASES(GImpl,);
    typedef OBC::Evolution<FlowAction, GImpl, FIMPL> EvolutionType;
public:
    TGaugeFlow_OBC(const std::string name);
    virtual ~TGaugeFlow_OBC(void) {};
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    virtual void setup(void);
    virtual void execute(void);
};

MODULE_REGISTER_TMP(WilsonFlow_OBC,ARG(TGaugeFlow_OBC<GIMPL, OBC::WilsonAction<GIMPL>>),MGradientFlow);
MODULE_REGISTER_TMP(SymanzikFlow_OBC,ARG(TGaugeFlow_OBC<GIMPL, SymanzikGaugeAction<GIMPL>>),MGradientFlow);
MODULE_REGISTER_TMP(CustomFlow_OBC,ARG(TGaugeFlow_OBC<GIMPL, PlaqPlusRectangleAction<GIMPL>>),MGradientFlow);

/******************************************************************************
 *                     TGaugeFlow_OBC implementation                          *
 ******************************************************************************/
template <typename GImpl, typename FlowAction>
TGaugeFlow_OBC<GImpl, FlowAction>::TGaugeFlow_OBC(const std::string name)
: Module<GaugeFlow_OBCPar>(name)
{}

template <typename GImpl, typename FlowAction>
std::vector<std::string> TGaugeFlow_OBC<GImpl, FlowAction>::getInput(void)
{
    return {par().gauge};
}

template <typename GImpl, typename FlowAction>
std::vector<std::string> TGaugeFlow_OBC<GImpl, FlowAction>::getOutput(void)
{
    if (par().steps < 0) {
        HADRONS_ERROR(Argument, "steps must be non-negative");
    }
    if (par().meas_interval <= 0) {
        HADRONS_ERROR(Argument, "meas_interval must be positive");
    }

    std::vector<std::string> out = {
        getName(), getName() + "_U",
        GaugeBoundary::metadataName(getName() + "_U")
    };

    // Retain fixed-step gauge snapshots only when requested. Adaptive flow has
    // a runtime-dependent schedule and publishes only its final field.
    if (par().save_history && par().maxTau.empty()) {
        const std::string initial = getName() + "_U_t0.00";
        out.push_back(initial);
        out.push_back(GaugeBoundary::metadataName(initial));
        for (int i = 1; i <= par().steps; ++i) {
            if ((i % par().meas_interval == 0) || (i == par().steps)) {
                const double ft = par().step_size * i;
                std::stringstream ftt;
                ftt << std::fixed << std::setprecision(2) << ft;
                const std::string field = getName() + "_U_t" + ftt.str();
                out.push_back(field);
                out.push_back(GaugeBoundary::metadataName(field));
            }
        }
    }
    if (par().save_rk_stages && par().save_history && par().maxTau.empty()) {
        for (int i = 0; i < par().steps; ++i) {
            const double ft = par().step_size * i;
            std::stringstream ftt;
            ftt << std::fixed << std::setprecision(2) << ft;
            const std::string w1 = getName() + "_W1_t" + ftt.str();
            const std::string w2 = getName() + "_W2_t" + ftt.str();
            out.push_back(w1);
            out.push_back(GaugeBoundary::metadataName(w1));
            out.push_back(w2);
            out.push_back(GaugeBoundary::metadataName(w2));
        }
    }

    return out;
}

template <typename GImpl, typename FlowAction>
void TGaugeFlow_OBC<GImpl, FlowAction>::setup(void)
{
    if (par().gaugeBC.empty()) {
        HADRONS_ERROR(Argument,
                      "gaugeBC must be explicitly set to 'periodic' or 'open'");
    }
    if (!GaugeBoundary::isValidType(par().gaugeBC) ||
        (par().gaugeBC == GaugeBoundary::unknown)) {
        HADRONS_ERROR(Argument, "gaugeBC must be 'periodic' or 'open'");
    }
    if (!std::isfinite(par().cG) || par().cG < 0.0) {
        HADRONS_ERROR(Argument,
                      "cG must be explicitly set to a finite non-negative value");
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
    if (par().save_rk_stages &&
        (!par().save_history || !par().maxTau.empty() ||
         par().meas_interval != 1)) {
        HADRONS_ERROR(Argument,
                      "save_rk_stages requires fixed-step save_history "
                      "with meas_interval = 1");
    }

    double maxTau = -1.0;
    if (!par().maxTau.empty()) {
        try {
            maxTau = std::stod(par().maxTau);
        } catch (const std::exception &) {
            HADRONS_ERROR(Argument, "maxTau must be a finite positive number");
        }
        if (!std::isfinite(maxTau) || maxTau <= 0.0) {
            HADRONS_ERROR(Argument, "maxTau must be a finite positive number");
        }
    }

    double cPlaq = 0.0;
    double cRect = 0.0;
    if constexpr (std::is_same_v<FlowAction,
                                 PlaqPlusRectangleAction<GImpl>>) {
        if (par().c_plaq.empty() || par().c_rect.empty()) {
            HADRONS_ERROR(Argument,
                          "CustomFlow_OBC requires both c_plaq and c_rect");
        }
        try {
            cPlaq = std::stod(par().c_plaq);
            cRect = std::stod(par().c_rect);
        } catch (const std::exception &) {
            HADRONS_ERROR(Argument,
                          "c_plaq and c_rect must be finite numbers");
        }
        if (!std::isfinite(cPlaq) || !std::isfinite(cRect)) {
            HADRONS_ERROR(Argument,
                          "c_plaq and c_rect must be finite numbers");
        }
    }

    const double effectiveCG = (par().gaugeBC == GaugeBoundary::open)
                             ? par().cG : 1.0;
    std::string actionName;
    if constexpr (std::is_same_v<FlowAction, OBC::WilsonAction<GImpl>>) {
        actionName = "WilsonAction";
    } else if constexpr (std::is_same_v<FlowAction,
                                        SymanzikGaugeAction<GImpl>>) {
        actionName = "SymanzikGaugeAction";
    } else {
        actionName = "PlaqPlusRectangleAction";
    }
    auto createBoundaryInfo =
        [this, effectiveCG, &actionName](const std::string &field,
                                        const double flowTime,
                                        const std::string &stage) {
            envCreate(GaugeBoundaryInfo, GaugeBoundary::metadataName(field), 1,
                      par().gaugeBC, effectiveCG, true,
                      true, actionName, par().step_size, flowTime, stage);
        };

    envCreateLat(GaugeField, getName() + "_U");
    const double finalFlowTime = par().maxTau.empty()
                               ? par().steps * par().step_size : -1.0;
    createBoundaryInfo(getName() + "_U", finalFlowTime, "field");

    if (par().save_history && par().maxTau.empty()) {
        const std::string initial = getName() + "_U_t0.00";
        envCreateLat(GaugeField, initial);
        createBoundaryInfo(initial, 0.0, "field");
        for (int i = 1; i <= par().steps; ++i) {
            if ((i % par().meas_interval == 0) || (i == par().steps)) {
                const double ft = par().step_size * i;
                std::stringstream ftt;
                ftt << std::fixed << std::setprecision(2) << ft;
                const std::string field = getName() + "_U_t" + ftt.str();
                envCreateLat(GaugeField, field);
                createBoundaryInfo(field, ft, "field");
            }
        }
    }
    if (par().save_rk_stages && par().save_history && par().maxTau.empty()) {
        for (int i = 0; i < par().steps; ++i) {
            const double ft = par().step_size * i;
            std::stringstream ftt;
            ftt << std::fixed << std::setprecision(2) << ft;
            const std::string w1 = getName() + "_W1_t" + ftt.str();
            const std::string w2 = getName() + "_W2_t" + ftt.str();
            envCreateLat(GaugeField, w1);
            createBoundaryInfo(w1, ft, "W1");
            envCreateLat(GaugeField, w2);
            createBoundaryInfo(w2, ft, "W2");
        }
    }

    envCreate(HadronsSerializable, getName(), 1, 0);
    if constexpr (std::is_same_v<FlowAction,
                                 PlaqPlusRectangleAction<GImpl>>) {
        envTmp(EvolutionType, "evolve", 1, envGetGrid(GaugeField),
               cPlaq, cRect, par().step_size, maxTau, par().step_size,
               par().gaugeBC, par().cG, this);
    } else {
        envTmp(EvolutionType, "evolve", 1, envGetGrid(GaugeField), 3.0,
               par().step_size, maxTau, par().step_size,
               par().gaugeBC, par().cG, this);
    }
}

template <typename GImpl, typename FlowAction>
void TGaugeFlow_OBC<GImpl, FlowAction>::execute(void)
{
    envGetTmp(EvolutionType, evolve);
    auto &out    = envGet(HadronsSerializable, getName());
    auto &result = out.template hold<OBC::GaugeResult>();
    auto &U      = envGet(GaugeField, par().gauge);
    auto &Uwf    = envGet(GaugeField, getName() + "_U");

    double maxTau = -1.0;
    LOG(Message) << "Setting up gauge gradient flow on '" << par().gauge
                 << "' with " << par().gaugeBC
                 << " temporal gauge boundary conditions";
    if (par().gaugeBC == GaugeBoundary::open) {
        LOG(Message) << " (cG=" << par().cG << ")";
    }
    LOG(Message) << " and " << par().steps << " step"
                 << ((par().steps != 1) ? "s." : ".") << std::endl;
    if (par().gaugeBC == GaugeBoundary::open &&
        std::abs(par().cG - 1.0) > 1.0e-14) {
        LOG(Message) << "OBC plaquette diagnostics use the flow action boundary "
                     << "weight (cG=" << par().cG << "); this differs from "
                     << "openQCD's conventional weighted plaquette when cG != 1."
                     << std::endl;
    }
    if (!par().maxTau.empty()) {
        maxTau = std::stod(par().maxTau);
        LOG(Message) << "Using adaptive algorithm with maxTau = "
                     << maxTau << std::endl;
    }
    std::cout << evolve.SG.LogParameters();

    Uwf = U;
    if (par().save_history && par().maxTau.empty()) {
        auto &U0 = envGet(GaugeField, getName() + "_U_t0.00");
        U0 = Uwf;
    }

    double flowt = 0.0;
    LOG(Message) << "Step 0 (tau = " << flowt << ")" << std::endl;
    LOG(Message) << "Compute observables" << std::endl;
    startTimer("observables");
    evolve.gauge_status(Uwf, result, flowt);
    stopTimer("observables");

    if (par().steps != 0) {
        if (maxTau > 0.0) {
            unsigned int step = 0;
            do {
                ++step;
                flowt += evolve.epsilon;
                LOG(Message) << "Step " << step << " (tau = " << flowt
                             << ")" << std::endl;
                startTimer("evolution");
                evolve.evolve_gauge_adaptive(Uwf);
                stopTimer("evolution");
                if (step % par().meas_interval == 0) {
                    LOG(Message) << "Compute observables" << std::endl;
                    startTimer("observables");
                    evolve.gauge_status(Uwf, result, flowt);
                    stopTimer("observables");
                }
            } while (evolve.taus < maxTau);
        } else {
            for (unsigned int step = 1;
                 step <= static_cast<unsigned int>(par().steps); ++step) {
                flowt += evolve.epsilon;
                LOG(Message) << "Step " << step << " (tau = " << flowt
                             << ")" << std::endl;
                startTimer("evolution");
                std::vector<GaugeField> &Wi = evolve.evolve_gauge(Uwf);
                stopTimer("evolution");

                if (step % par().meas_interval == 0) {
                    LOG(Message) << "Compute observables" << std::endl;
                    startTimer("observables");
                    evolve.gauge_status(Uwf, result, flowt);
                    stopTimer("observables");
                }
                if (par().save_history &&
                    ((step % par().meas_interval == 0) ||
                     (step == static_cast<unsigned int>(par().steps)))) {
                    std::stringstream ftt;
                    ftt << std::fixed << std::setprecision(2) << flowt;
                    auto &Us = envGet(GaugeField,
                                      getName() + "_U_t" + ftt.str());
                    Us = Uwf;
                }
                if (par().save_rk_stages) {
                    const double startTime = par().step_size * (step - 1);
                    std::stringstream ftt;
                    ftt << std::fixed << std::setprecision(2) << startTime;
                    auto &W1 = envGet(GaugeField,
                                      getName() + "_W1_t" + ftt.str());
                    auto &W2 = envGet(GaugeField,
                                      getName() + "_W2_t" + ftt.str());
                    W1 = Wi[1];
                    W2 = Wi[2];
                }
            }
        }
    }

    auto &finalMetadata = envGet(
        GaugeBoundaryInfo, GaugeBoundary::metadataName(getName() + "_U"));
    finalMetadata.flowTime = (maxTau > 0.0) ? evolve.taus : flowt;

    saveResult(par().output, "gauge_obs", result);
}

END_MODULE_NAMESPACE
END_HADRONS_NAMESPACE

#endif // Hadrons_MGradientFlow_GaugeFlow_OBC_hpp_
