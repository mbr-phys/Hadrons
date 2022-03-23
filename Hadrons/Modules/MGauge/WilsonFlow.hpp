/*
 * WilsonFlow.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
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
#ifndef Hadrons_MGauge_WilsonFlow_hpp_
#define Hadrons_MGauge_WilsonFlow_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                               Wilson Flow                                  *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MGauge)

class WilsonFlowPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(WilsonFlowPar,
                                    std::string, gauge,
                                    int, steps,
                                    double, step_size,
                                    int, meas_interval,
                                    std::string, maxTau); 
};

template <typename GImpl>
class TWilsonFlow: public Module<WilsonFlowPar>
{
public:
    GAUGE_TYPE_ALIASES(GImpl,);
public:
    // constructor
    TWilsonFlow(const std::string name);
    // destructor
    virtual ~TWilsonFlow(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(WilsonFlow, TWilsonFlow<GIMPL>, MGauge);

/******************************************************************************
 *                     TWilsonFlow implementation                          *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename GImpl>
TWilsonFlow<GImpl>::TWilsonFlow(const std::string name)
: Module<WilsonFlowPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename GImpl>
std::vector<std::string> TWilsonFlow<GImpl>::getInput(void)
{
    std::vector<std::string> in = {par().gauge};
    
    return in;
}

template <typename GImpl>
std::vector<std::string> TWilsonFlow<GImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename GImpl>
void TWilsonFlow<GImpl>::setup(void)
{
    envCreateLat(GaugeField, getName());
    envTmpLat(GaugeField, "Umu");
}

// execution ///////////////////////////////////////////////////////////////////
template <typename GImpl>
void TWilsonFlow<GImpl>::execute(void)
{
    LOG(Message) << "Setting up Wilson Flow on '" << par().gauge << "' with " << par().steps
                 << " step" << ((par().steps > 1) ? "s." : ".") << std::endl;

    int mTau = -1;
    if(!par().maxTau.empty()) {
        LOG(Message) << "Using adaptive algorithm with " << std::endl;
        mTau = std::stoi(par().maxTau);
    }
    WilsonFlow<GImpl>  Wflow(par().steps, par().step_size, par().meas_interval);
    auto               &U   = envGet(GaugeField, par().gauge);
    auto               &Uwf = envGet(GaugeField, getName());

    envGetTmp(GaugeField, Umu);
    Umu = U;
    LOG(Message) << "flow time = 0," 
                 << "         plaquette = " << WilsonLoops<GImpl>::avgPlaquette(U) << std::endl
                 << "    energy density = " << Wflow.energyDensityPlaquette(0,U)   << std::endl;

    Uwf = U;
    if (mTau > 0) {
        WF.smear_adaptive(Uwf, Umu, mTau);
    } else {
        WF.smear(Uwf, Umu);
    }
    LOG(Message) << "flow time = " << par().steps * par().step_size << " ,"
                 << "         plaquette = " << WilsonLoops<GImpl>::avgPlaquette(Uwf)         << std::endl
                 << "    energy density = " << Wflow.energyDensityPlaquette(par().steps,Uwf) << std::endl;
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MGauge_WilsonFlow_hpp_
