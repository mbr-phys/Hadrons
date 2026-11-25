/*
 * GaugeFlow.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
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
#ifndef Hadrons_MGradientFlow_GaugeFlow_hpp_
#define Hadrons_MGradientFlow_GaugeFlow_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Serialization.hpp>
#include <Hadrons/Modules/MGradientFlow/Utils.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                        Gauge Field Gradient Flow                           *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MGradientFlow)

class GaugeFlowPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(GaugeFlowPar,
                                    std::string, output,
                                    std::string, gauge,
                                    int, steps,
                                    double, step_size,
                                    int, meas_interval,
                                    bool, save_history,
                                    bool, save_rk_stages,
                                    std::string, maxTau,
                                    std::string, c_plaq,
                                    std::string, c_rect); 

    GaugeFlowPar(void)
    : save_history(false)
    , save_rk_stages(false)
    {}
};

template <typename GImpl,typename FlowAction>
class TGaugeFlow: public Module<GaugeFlowPar>
{
public:
    GAUGE_TYPE_ALIASES(GImpl,);
    typedef Evolution<FlowAction, GImpl, FIMPL> EvolutionType;
public:
    // constructor
    TGaugeFlow(const std::string name);
    // destructor
    virtual ~TGaugeFlow(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    // names of the saved gauge snapshots / RK stages (single source for
    // getOutput, setup and execute, built from the integer step)
    std::string snapshotName(const int step);
    std::string stageName(const int stage, const int step);
    std::vector<std::string> savedNames(void);
};

MODULE_REGISTER_TMP(WilsonFlow, ARG(TGaugeFlow<GIMPL,WilsonAction<GIMPL>>), MGradientFlow);
MODULE_REGISTER_TMP(SymanzikFlow, ARG(TGaugeFlow<GIMPL,SymanzikGaugeAction<GIMPL>>), MGradientFlow);
MODULE_REGISTER_TMP(ZeuthenFlow, ARG(TGaugeFlow<GIMPL,ZeuthenGaugeAction<GIMPL>>), MGradientFlow);
MODULE_REGISTER_TMP(CustomFlow, ARG(TGaugeFlow<GIMPL,PlaqPlusRectangleAction<GIMPL>>), MGradientFlow);

/******************************************************************************
 *                     TGaugeFlow implementation                          *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename GImpl,typename FlowAction>
TGaugeFlow<GImpl,FlowAction>::TGaugeFlow(const std::string name)
: Module<GaugeFlowPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename GImpl,typename FlowAction>
std::vector<std::string> TGaugeFlow<GImpl,FlowAction>::getInput(void)
{
    std::vector<std::string> in = {par().gauge};
    
    return in;
}

template <typename GImpl,typename FlowAction>
std::vector<std::string> TGaugeFlow<GImpl,FlowAction>::getOutput(void)
{
    std::vector<std::string> out = {getName(),getName()+"_U"};
    std::vector<std::string> saved = savedNames();

    out.insert(out.end(), saved.begin(), saved.end());
    
    return out;
}

// output names ////////////////////////////////////////////////////////////////
template <typename GImpl,typename FlowAction>
std::string TGaugeFlow<GImpl,FlowAction>::snapshotName(const int step)
{
    return getName()+"_U_t"+flowTimeString(step, par().step_size);
}

// RK stage W1/W2 of the step starting at flow time step*step_size
template <typename GImpl,typename FlowAction>
std::string TGaugeFlow<GImpl,FlowAction>::stageName(const int stage, const int step)
{
    return getName()+"_W"+std::to_string(stage)+"_t"+flowTimeString(step, par().step_size);
}

template <typename GImpl,typename FlowAction>
std::vector<std::string> TGaugeFlow<GImpl,FlowAction>::savedNames(void)
{
    std::vector<std::string> names;

    checkFlowSteps(par().steps, par().step_size, par().meas_interval);
    // Retain fixed-step gauge snapshots only when requested. Adaptive flow has
    // a runtime-dependent schedule and therefore only publishes its final
    // field through getName()+"_U".
    if (par().save_history && par().maxTau.empty()) {
        names.push_back(snapshotName(0));
        for (int step: flowSaveSteps(par().steps, par().step_size, par().meas_interval)) {
            names.push_back(snapshotName(step));
        }
    }
    if (par().save_rk_stages && par().save_history && par().maxTau.empty()) {
        for (int step = 0; step < par().steps; step++) {
            names.push_back(stageName(1, step));
            names.push_back(stageName(2, step));
        }
    }
    checkUniqueNames(names);

    return names;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename GImpl,typename FlowAction>
void TGaugeFlow<GImpl,FlowAction>::setup(void)
{
    if (par().save_rk_stages &&
        (!par().save_history || !par().maxTau.empty() || par().meas_interval != 1)) {
        HADRONS_ERROR(Argument, "save_rk_stages requires fixed-step save_history with meas_interval = 1");
    }
    // checked before any std::stod(c_plaq/c_rect) below
    if constexpr (std::is_same_v<FlowAction, PlaqPlusRectangleAction<GImpl>>) {
        if (par().c_plaq.empty() || par().c_rect.empty()) {
            HADRONS_ERROR(Argument, "to use PlaqPlusRectangleAction (CustomFlow), pass some value to both c_plaq and c_rect");
        }
    }

    double mTau = -1.0;
    if(!par().maxTau.empty()) {
        mTau = std::stod(par().maxTau);
    }
    envCreateLat(GaugeField, getName()+"_U");
    for (auto &n: savedNames()) {
        envCreateLat(GaugeField, n);
    }
    envCreate(HadronsSerializable, getName(), 1, 0);
    if constexpr (std::is_same_v<FlowAction, PlaqPlusRectangleAction<GImpl>>) {
        envTmp(EvolutionType, "evolve", 1, envGetGrid(GaugeField), 
            std::stod(par().c_plaq), std::stod(par().c_rect), par().step_size, 
            mTau, par().step_size, this);
    } else {
        envTmp(EvolutionType, "evolve", 1, envGetGrid(GaugeField), 3.0, 
            par().step_size, mTau, par().step_size, this);
    }
    // create action -> if c_plaq and c_rect are used, called PlaqPlusRectangleAction
    if constexpr (std::is_same_v<FlowAction, PlaqPlusRectangleAction<GImpl>>) {
        envTmp(FlowAction, "action", 1, std::stod(par().c_plaq),std::stod(par().c_rect));  
    } else {
        envTmp(FlowAction, "action", 1, 3.0);
    }
}

// execution ///////////////////////////////////////////////////////////////////
template <typename GImpl,typename FlowAction>
void TGaugeFlow<GImpl,FlowAction>::execute(void)
{
    envGetTmp(EvolutionType, evolve);
    envGetTmp(FlowAction, action);
    auto &out    = envGet(HadronsSerializable, getName());
    auto &result = out.template hold<GaugeResult>();
    auto &U   = envGet(GaugeField, par().gauge);
    auto &Uwf = envGet(GaugeField, getName()+"_U");

    double mTau = -1.0;
    LOG(Message) << "Setting up gauge gradient flow on '" << par().gauge << "' with " << par().steps
                 << " step" << ((par().steps != 1) ? "s." : ".") << std::endl;
    if(!par().maxTau.empty()) {
        LOG(Message) << "Using adaptive algorithm with maxTau = " << par().maxTau << std::endl;
        mTau = std::stod(par().maxTau);
    }
    std::cout << action.LogParameters();
    
    Uwf = U;
    if (par().save_history && par().maxTau.empty()) {
        auto &U0 = envGet(GaugeField, snapshotName(0));
        U0 = Uwf;
    }
    double flowt = 0.0;
    LOG(Message) << "Step 0 (tau = "<< flowt << ")" << std::endl;
    LOG(Message) << "Compute observables" << std::endl;
    startTimer("observables");
    evolve.gauge_status(Uwf,result,flowt);
    stopTimer("observables"); 
    // if steps = 0, give the status of gauge field without flowing
    if (par().steps != 0) { 
        if (mTau > 0) {
            unsigned int step = 0;
            do {
                step++;
                flowt += evolve.epsilon;
                LOG(Message) << "Step " << step << " (tau = "<< flowt << ")" << std::endl;
                startTimer("evolution");
                evolve.evolve_gauge_adaptive(Uwf);
                stopTimer("evolution");
                if (step % par().meas_interval == 0) {
                    LOG(Message) << "Compute observables" << std::endl;
                    startTimer("observables");
                    evolve.gauge_status(Uwf,result,flowt);
                    stopTimer("observables"); 
                }
            } while (evolve.taus < mTau);
        } else {
            for (int step = 1; step <= par().steps; step++) {
                // flow time from the integer step (same as the object names)
                flowt = par().step_size*step;
                LOG(Message) << "Step " << step << " (tau = "<< flowt << ")" << std::endl;
                startTimer("evolution");
                std::vector<GaugeField> &Wi = evolve.evolve_gauge(Uwf);
                stopTimer("evolution");
                if (step % par().meas_interval == 0) {
                    LOG(Message) << "Compute observables" << std::endl;
                    startTimer("observables");
                    evolve.gauge_status(Uwf,result,flowt);
                    stopTimer("observables");
                }
                if (par().save_history &&
                    isFlowSaveStep(step, par().steps, par().meas_interval)) {
                    auto &Us = envGet(GaugeField, snapshotName(step));
                    Us = Uwf;
                }
                if (par().save_rk_stages) {
                    auto &W1 = envGet(GaugeField, stageName(1, step - 1));
                    auto &W2 = envGet(GaugeField, stageName(2, step - 1));
                    W1 = Wi[1];
                    W2 = Wi[2];
                }
            }
        }
    }
    saveResult(par().output,"gauge_obs",result);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MGradientFlow_GaugeFlow_hpp_
