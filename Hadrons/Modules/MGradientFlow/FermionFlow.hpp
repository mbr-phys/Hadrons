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
                                    std::string, q3,
                                    std::string, gauge,
                                    int, steps,
                                    double, step_size,
                                    int, meas_interval,
                                    std::string, maxTau); 
};

template <typename FImpl1,typename FImpl2,typename FImpl3,typename GImpl,typename FlowAction>
class TFermionFlow: public Module<FermionFlowPar>
{
public:
    FERM_TYPE_ALIASES(FImpl1, 1);
    FERM_TYPE_ALIASES(FImpl2, 2);
    FERM_TYPE_ALIASES(FImpl2, 3);
    INHERIT_GIMPL_TYPES(GImpl);
    class GaugeResult : Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(GaugeResult,
                                        std::vector<RealD>, flowtime,
                                        std::vector<RealD>, plaquette,
                                        std::vector<RealD>, rectangle,
                                        std::vector<RealD>, clover,
                                        std::vector<RealD>, topcharge,
                                        std::vector<RealD>, action);
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
    void status(double time, GaugeField &Umu, auto &result, int step);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(WilsonFermionFlow,ARG(TFermionFlow<FIMPL,FIMPL,FIMPL,GIMPL,WilsonGaugeAction<GIMPL>>),MGradientFlow);

/******************************************************************************
 *                     TFermionFlow implementation                          *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl1,typename FImpl2,typename FImpl3,typename GImpl,typename FlowAction>
TFermionFlow<FImpl1,FImpl2,FImpl3,GImpl,FlowAction>::TFermionFlow(const std::string name)
: Module<FermionFlowPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl1,typename FImpl2,typename FImpl3,typename GImpl,typename FlowAction>
std::vector<std::string> TFermionFlow<FImpl1,FImpl2,FImpl3,GImpl,FlowAction>::getInput(void)
{
    std::vector<std::string> in = {par().gauge,par().q1};
    if (!((par().q1.compare(par().q2) == 0) || par().q2.empty())) in.push_back(par().q2);
    if (!((par().q1.compare(par().q3) == 0) || (par().q2.compare(par().q3) == 0) || par().q3.empty())) in.push_back(par().q3);
    
    return in;
}

template <typename FImpl1,typename FImpl2,typename FImpl3,typename GImpl,typename FlowAction>
std::vector<std::string> TFermionFlow<FImpl1,FImpl2,FImpl3,GImpl,FlowAction>::getOutput(void)
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
            if (!((par().q1.compare(par().q3) == 0) || (par().q2.compare(par().q3) == 0) || par().q3.empty())) {
                out.push_back(getName()+"_q3_"+st.str());}
        }
    }
    
    std::cout << "TESTING: q1 - " << par().q1 << " q2 - " << par().q2 << " q3 - " << par().q3 << " "; 
    for (int j = 0; j < out.size(); j++) std::cout << out[j] << " ";
    std::cout << "DONE." << std::endl;
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl1,typename FImpl2,typename FImpl3,typename GImpl,typename FlowAction>
void TFermionFlow<FImpl1,FImpl2,FImpl3,GImpl,FlowAction>::setup(void)
{
    envCreateLat(GaugeField, getName()+"_U");

    envTmpLat(PropagatorField1, "q1wf");
    envTmpLat(PropagatorField2, "q2wf");
    envTmpLat(PropagatorField3, "q3wf");

    for (int i = 1; i <= par().steps; i++) 
    {
        if ((i % par().meas_interval == 0) || (i == par().steps)) {
            std::stringstream st;
            st << i;
            envCreateLat(PropagatorField1, getName()+"_q1_"+st.str());
            if (!((par().q1.compare(par().q2) == 0) || par().q2.empty())) {
                envCreateLat(PropagatorField2, getName()+"_q2_"+st.str());}
            if (!((par().q1.compare(par().q3) == 0) || (par().q2.compare(par().q3) == 0) || par().q3.empty())) {
                envCreateLat(PropagatorField3, getName()+"_q3_"+st.str());}
        }
    }
    envCreate(HadronsGenericSerializable, getName(), 1, 0);
}

template <typename FImpl1,typename FImpl2,typename FImpl3,typename GImpl,typename FlowAction>
void TFermionFlow<FImpl1,FImpl2,FImpl3,GImpl,FlowAction>::status(double time, GaugeField &Umu, auto &Uresult, int step)
{
    RealD Q = WilsonLoops<GImpl>::TopologicalCharge(Umu);
    RealD plaq = WilsonLoops<GImpl>::avgPlaquette(Umu);
    RealD rect = WilsonLoops<GImpl>::avgRectangle(Umu);
    RealD clov = avgClover<GImpl,ComplexField,GaugeField,GaugeLinkField>(Umu);
    RealD act = SG.S(Umu);

//    if (par().output.length()) {
    Uresult.flowtime[step-1]  = time;
    Uresult.plaquette[step-1] = plaq;
    Uresult.rectangle[step-1] = rect;
    Uresult.clover[step-1]    = clov;
    Uresult.topcharge[step-1] =    Q;
    Uresult.action[step-1]    =  act;
 //   } else {
 //       LOG(Message) << "flow time = " << std::setprecision(3) << std::fixed << time
 //                    << " top. charge: " << std::setprecision(16) << std::scientific << Q                      
 //                    << " plaquette: " << plaq
 //                    << " rectangle: " << rect
 //                    << " clover: " << clov
 //                    << " action: " << act << std::endl;
 //   }
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl1,typename FImpl2,typename FImpl3,typename GImpl,typename FlowAction>
void TFermionFlow<FImpl1,FImpl2,FImpl3,GImpl,FlowAction>::execute(void)
{
    std::string type = SG.action_name();
    std::string ga = "GaugeAction";
    std::string::size_type i = type.find(ga);
    if (i != std::string::npos) {
        type.erase(i, ga.length());
    }

    LOG(Message) << "Setting up " << type << " Fermion Flow on '" << par().gauge << "' Gauge Field and " 
                 << par().q1 << ", " << par().q2 << " Fermion Propagators with " 
                 << par().steps << " step" << ((par().steps > 1) ? "s." : ".") << std::endl;

    RealD mTau = -1.0;
    if(!par().maxTau.empty()) {
        // is this working for fermion flow?
        LOG(Message) << "Using adaptive algorithm with maxTau = " << par().maxTau << std::endl;
        mTau = (RealD)std::stoi(par().maxTau);
    }

    auto &out     = envGet(HadronsGenericSerializable, getName());
    auto &Uresult = out.template hold<GaugeResult>();

    Uresult.flowtime.resize(par().steps);
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

    std::string q3ref = par().q3;
    if (((par().q1.compare(par().q3) == 0) || (par().q2.compare(par().q3) == 0) || par().q3.empty())) q3ref = par().q1;
    auto &q3   = envGet(PropagatorField3, q3ref);
    envGetTmp(PropagatorField3, q3wf);
    q3wf = q3;

    double time = 0;
//    status(time,U,Uresult,0);
    Evolution<GImpl,FlowAction> evolve(3.0, par().step_size, mTau, par().step_size);
    if (mTau > 0) {
        unsigned int step = 0;
        do {
            step++;
            status(evolve.taus,Uwf,Uresult,step);
            evolve.template evolve_prop_adaptive<PropagatorField1,PropagatorField2,PropagatorField3>(Uwf,q1wf,q2wf,q3wf);
            if (step % par().meas_interval == 0) {
                // need adaptive way to create right number of output PropagatorFields:
                // current algorithm should output less than standard method so we shouldn't lose
                // any steps, but will output unneccessary empty propagators to be contracted right now
                std::stringstream st;
                st << step;
                auto &q1i = envGet(PropagatorField1, getName()+"_q1_"+st.str());
                q1i = q1wf;
                if (!((par().q1.compare(par().q2) == 0) || par().q2.empty())) {
                    auto &q2i = envGet(PropagatorField2, getName()+"_q2_"+st.str());
                    q2i = q2wf;
                }
                if (!((par().q1.compare(par().q3) == 0) || (par().q2.compare(par().q3) == 0) || par().q3.empty())) {
                    auto &q3i = envGet(PropagatorField3, getName()+"_q3_"+st.str());
                    q3i = q3wf;
                }
            }
        } while (evolve.taus < mTau);
    } else {
        for (unsigned int step = 1; step <= par().steps; step++) {
            status(step*par().step_size,Uwf,Uresult,step);
            evolve.template evolve_prop<PropagatorField1,PropagatorField2,PropagatorField3>(Uwf,q1wf,q2wf,q3wf);
            if ((step % par().meas_interval == 0) || (step == par().steps)) {
                std::stringstream st;
                st << step;
                auto &q1i = envGet(PropagatorField1, getName()+"_q1_"+st.str());
                q1i = q1wf;
                if (!((par().q1.compare(par().q2) == 0) || par().q2.empty())) {
                    auto &q2i = envGet(PropagatorField2, getName()+"_q2_"+st.str());
                    q2i = q2wf;
                }
                if (!((par().q1.compare(par().q3) == 0) || (par().q2.compare(par().q3) == 0) || par().q3.empty())) {
                    auto &q3i = envGet(PropagatorField3, getName()+"_q3_"+st.str());
                    q3i = q3wf;
                }
            }
        }
    }
    saveResult(par().output,"gauge_obs",Uresult);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MGradientFlow_FermionFlow_hpp_
