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
#include <Hadrons/Modules/MGradientFlow/Utils.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                      Propagator Field Gradient Flow                        *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MGradientFlow)

class FermionFlowPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(FermionFlowPar,
                                    std::string, output,
                                    std::string, q1,
                                    std::string, q2,
                                    std::string, gauge,
                                    int, steps,
                                    double, step_size,
                                    int, meas_interval,
                                    std::string, maxTau);
};

template <typename FImpl1,typename FImpl2,typename GImpl,typename FlowAction>
class TFermionFlow: public Module<FermionFlowPar>
{
public:
    FERM_TYPE_ALIASES(FImpl1, 1);
    FERM_TYPE_ALIASES(FImpl2, 2);
    INHERIT_GIMPL_TYPES(GImpl);
    class GaugeResult : Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(GaugeResult,
                                        std::vector<double>, plaquette,
                                        std::vector<double>, rectangle,
                                        std::vector<double>, clover,
                                        std::vector<double>, topcharge,
                                        std::vector<double>, action);
    };
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
    // action
    FlowAction SG = FlowAction(3.0);
    // gauge measurements at each flow time
    // void status(double time, GaugeField &Umu, auto &result, int step);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(WilsonFermionFlow,ARG(TFermionFlow<FIMPL,FIMPL,GIMPL,WilsonGaugeAction<GIMPL>>),MGradientFlow);

/******************************************************************************
 *                     TFermionFlow implementation                          *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl1,typename FImpl2,typename GImpl,typename FlowAction>
TFermionFlow<FImpl1,FImpl2,GImpl,FlowAction>::TFermionFlow(const std::string name)
: Module<FermionFlowPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl1,typename FImpl2,typename GImpl,typename FlowAction>
std::vector<std::string> TFermionFlow<FImpl1,FImpl2,GImpl,FlowAction>::getInput(void)
{
    std::vector<std::string> in = {par().gauge,par().q1};
    if (!((par().q1.compare(par().q2) == 0) || par().q2.empty())) in.push_back(par().q2);
    
    return in;
}

template <typename FImpl1,typename FImpl2,typename GImpl,typename FlowAction>
std::vector<std::string> TFermionFlow<FImpl1,FImpl2,GImpl,FlowAction>::getOutput(void)
{
    std::vector<std::string> out = {getName(),getName()+"_U"};
    for (int i = 1; i <= par().steps; i++) 
    {
        if ((i % par().meas_interval == 0) || (i == par().steps)) {
            std::stringstream st;
            st << i;
            out.push_back(getName()+"_q1_"+st.str());
            if (!((par().q1.compare(par().q2) == 0) || par().q2.empty())) {
                out.push_back(getName()+"_q2_"+st.str());}
        }
    }
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl1,typename FImpl2,typename GImpl,typename FlowAction>
void TFermionFlow<FImpl1,FImpl2,GImpl,FlowAction>::setup(void)
{
    envCreateLat(GaugeField, getName()+"_U");

    envTmpLat(PropagatorField1, "q1wf");
    envTmpLat(PropagatorField2, "q2wf");

    for (int i = 1; i <= par().steps; i++) 
    {
        if ((i % par().meas_interval == 0) || (i == par().steps)) {
            std::stringstream st;
            st << i;
            envCreateLat(PropagatorField1, getName()+"_q1_"+st.str());
            if (!((par().q1.compare(par().q2) == 0) || par().q2.empty())) {
                envCreateLat(PropagatorField2, getName()+"_q2_"+st.str());}
        }
    }
    envCreate(HadronsGenericSerializable, getName(), 1, 0);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl1,typename FImpl2,typename GImpl,typename FlowAction>
void TFermionFlow<FImpl1,FImpl2,GImpl,FlowAction>::execute(void)
{
    std::string type = SG.action_name();
    std::string ga = "GaugeAction";
    std::string::size_type i = type.find(ga);
    if (i != std::string::npos) {
        type.erase(i, ga.length());
    }

    LOG(Message) << "Setting up " << type << " Fermion Flow on '" << par().gauge << "' Gauge Field and "  << par().q1 
                 << (!(par().q2.empty()) ? ", "+par().q2+" Fermion Propagators " : " Fermion Propagator ")
                 << "with " << par().steps << " step" << ((par().steps > 1) ? "s." : ".") << std::endl;

    double mTau = -1.0;
    if(!par().maxTau.empty()) {
        LOG(Message) << "Using adaptive algorithm with maxTau = " << par().maxTau << std::endl;
        mTau = (double)std::stoi(par().maxTau);
    }

    auto &out     = envGet(HadronsGenericSerializable, getName());
    auto &Uresult = out.template hold<GaugeResult>();

    Uresult.plaquette.resize(par().steps);
    Uresult.rectangle.resize(par().steps);
    Uresult.clover.resize(par().steps);
    Uresult.topcharge.resize(par().steps);
    Uresult.action.resize(par().steps);

    auto &U   = envGet(GaugeField, par().gauge);
    auto &Uwf = envGet(GaugeField, getName()+"_U");
    Uwf = U;

    auto &q1   = envGet(PropagatorField1, par().q1);
    envGetTmp(PropagatorField1, q1wf);
    q1wf = q1;

    std::string q2ref = par().q2;
    if (((par().q1.compare(par().q2) == 0) || par().q2.empty())) q2ref = par().q1;
    auto &q2   = envGet(PropagatorField2, q2ref);
    envGetTmp(PropagatorField2, q2wf);
    q2wf = q2;

    double time = 0;
    Evolution<GImpl,FlowAction> evolve(3.0, par().step_size, mTau, par().step_size);
    if (mTau > 0) {
        unsigned int step = 0;
        do {
            step++;
            evolve.template evolve_prop_adaptive<PropagatorField1,PropagatorField2>(Uwf,q1wf,q2wf); 
            evolve.gauge_status(Uwf,Uresult,step-1);
            if (step % par().meas_interval == 0) {
                // need adaptive way to create right number of output PropagatorFields:
                // current algorithm should output less than old method so we shouldn't lose
                // any steps, but will output unneccessary empty propagators right now
                std::stringstream st;
                st << step;
                auto &q1i = envGet(PropagatorField1, getName()+"_q1_"+st.str());
                q1i = q1wf;
                if (!((par().q1.compare(par().q2) == 0) || par().q2.empty())) {
                    auto &q2i = envGet(PropagatorField2, getName()+"_q2_"+st.str());
                    q2i = q2wf;
                }
            }
        } while (evolve.taus < mTau);
    } else {
        for (unsigned int step = 1; step <= par().steps; step++) {
            evolve.template evolve_prop<PropagatorField1,PropagatorField2>(Uwf,q1wf,q2wf); 
            evolve.gauge_status(Uwf,Uresult,step-1);
            if ((step % par().meas_interval == 0) || (step == par().steps)) {
                std::stringstream st;
                st << step;
                auto &q1i = envGet(PropagatorField1, getName()+"_q1_"+st.str());
                q1i = q1wf; 
                if (!((par().q1.compare(par().q2) == 0) || par().q2.empty())) {
                    auto &q2i = envGet(PropagatorField2, getName()+"_q2_"+st.str());
                    q2i = q2wf;
                }
            }
        }
    }
    saveResult(par().output,"gauge_obs",Uresult);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MGradientFlow_FermionFlow_hpp_
