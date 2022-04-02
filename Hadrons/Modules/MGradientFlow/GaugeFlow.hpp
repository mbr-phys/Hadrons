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
                                    std::string, maxTau); 
};

template <typename GImpl,typename FlowAction>
class TGaugeFlow: public Module<GaugeFlowPar>
{
public:
    INHERIT_GIMPL_TYPES(GImpl);
    class Result : Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(Result,
                                        std::vector<RealD>, flowtime,
                                        std::vector<RealD>, plaquette,
                                        std::vector<RealD>, rectangle,
                                        std::vector<RealD>, clover,
                                        std::vector<RealD>, topcharge,
                                        std::vector<RealD>, action);
    };
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
    // action
    FlowAction SG = FlowAction(3.0);
    // gauge observable measurements at flow time
    void status(double time, GaugeField &Umu, Result &result, int step);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(WilsonFlow, ARG(TGaugeFlow<GIMPL,WilsonGaugeAction<GIMPL>>), MGradientFlow);
MODULE_REGISTER_TMP(SymanzikFlow, ARG(TGaugeFlow<GIMPL,SymanzikGaugeAction<GIMPL>>), MGradientFlow);
MODULE_REGISTER_TMP(ZeuthenFlow, ARG(TGaugeFlow<GIMPL,ZeuthenGaugeAction<GIMPL>>), MGradientFlow);

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
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename GImpl,typename FlowAction>
void TGaugeFlow<GImpl,FlowAction>::setup(void)
{
    envCreateLat(GaugeField, getName());
}

template <typename GImpl,typename FlowAction>
void TGaugeFlow<GImpl,FlowAction>::status(double time, GaugeField &Umu, Result &result, int step)
{
    RealD Q = WilsonLoops<GImpl>::TopologicalCharge(Umu);
    RealD plaq = WilsonLoops<GImpl>::avgPlaquette(Umu);
    RealD rect = WilsonLoops<GImpl>::avgRectangle(Umu);
    RealD clov = avgClover<GImpl,ComplexField,GaugeField,GaugeLinkField>(Umu);
    RealD act = SG.S(Umu);

    if (par().output.length()) {
        result.flowtime[step]  = time;
        result.plaquette[step] = plaq;
        result.rectangle[step] = rect;
        result.clover[step]    = clov;
        result.topcharge[step] =    Q;
        result.action[step]    =  act;
    } else {
        LOG(Message) << "flow time = " << std::setprecision(3) << std::fixed << time 
                     << " top. charge: " << std::setprecision(16) << std::scientific << Q
                     << " plaquette: " << plaq
                     << " rectangle: " << rect
                     << " clover: " << clov
                     << " action: " << act << std::endl;
    }
}

// execution ///////////////////////////////////////////////////////////////////
template <typename GImpl,typename FlowAction>
void TGaugeFlow<GImpl,FlowAction>::execute(void)
{
    std::string type = SG.action_name();
    std::string ga = "GaugeAction";
    std::string::size_type i = type.find(ga);
    if (i != std::string::npos) {
        type.erase(i, ga.length());
    }

    LOG(Message) << "Setting up " << type << " Flow on '" << par().gauge << "' with " << par().steps
                 << " step" << ((par().steps > 1) ? "s." : ".") << std::endl;

    RealD mTau = -1.0;
    if(!par().maxTau.empty()) {
        LOG(Message) << "Using adaptive algorithm with maxTau = " << par().maxTau << std::endl;
        mTau = (RealD)std::stoi(par().maxTau);
    }

    Result result;
    result.flowtime.resize(1+par().steps);
    result.plaquette.resize(1+par().steps);
    result.rectangle.resize(1+par().steps);
    result.clover.resize(1+par().steps);
    result.topcharge.resize(1+par().steps);
    result.action.resize(1+par().steps);

    auto &U   = envGet(GaugeField, par().gauge);
    auto &Uwf = envGet(GaugeField, getName());

    Uwf = U;
    double time = 0;
    status(time,U,result,0);
    Evolution<GImpl,FlowAction> evolve(3.0, par().step_size, mTau, par().step_size);
    if (mTau > 0) {
        unsigned int step = 0;
        do {
            step++;
            evolve.evolve_gauge_adaptive(Uwf);
            if (step % par().meas_interval == 0) {
                status(evolve.taus,Uwf,result,step);
            }
        } while (evolve.taus < mTau);
    } else {
        for (unsigned int step = 1; step <= par().steps; step++) {
            evolve.evolve_gauge(Uwf);
            if (step % par().meas_interval == 0) {
                status(step*par().step_size,Uwf,result,step);
            }
        }
    }
    if (par().output.length()) {
        saveResult(par().output,"flow",result);
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MGradientFlow_GaugeFlow_hpp_
