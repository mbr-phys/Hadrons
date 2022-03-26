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
private:
    typedef typename GImpl::GaugeLinkField GaugeMat;
    typedef typename GImpl::GaugeField GaugeLorentz;
    typedef typename GImpl::ComplexField ComplexField;
    // clover
    void siteClover(ComplexField &Clov, const GaugeLorentz &U);
    RealD avgClover(const GaugeLorentz &Umu);
    void status(double time, GaugeField &Umu, WilsonGaugeAction<GImpl> &SG);
    void evolve_step(GaugeField &U, WilsonGaugeAction<GImpl> &SG);
    void evolve_step_adaptive(GaugeField &U, RealD maxTau, RealD &epsilon, RealD &taus, WilsonGaugeAction<GImpl> &SG);
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

// clover //////////////////////////////////////////////////////////////////////
template <typename GImpl>
void TWilsonFlow<GImpl>::siteClover(ComplexField &Clov, const GaugeLorentz &U)
{
    GaugeMat Fmn(U.Grid()), Cmn(U.Grid()), scaledUnit(U.Grid()), Umu(U.Grid());
    Clov = Zero();
    for (int mu = 1; mu < Nd; mu++) {
        for (int nu = 0; nu < mu; nu++) {
            Umu = PeekIndex<LorentzIndex>(U, mu);
            scaledUnit = (1.0/Nc) * (adj(Umu) * Umu);
            WilsonLoops<GImpl>::FieldStrength(Fmn, U, mu, nu);
            Cmn = Fmn - trace(Fmn) * scaledUnit;
            Clov = Clov - trace(Cmn * Cmn);
        }
    }
}

template <typename GImpl>
RealD TWilsonFlow<GImpl>::avgClover(const GaugeLorentz &Umu) 
{
    ComplexField Clov(Umu.Grid());

    siteClover(Clov, Umu);
    auto Tc = sum(Clov);
    auto c = TensorRemove(Tc);

    double vol = Umu.Grid()->gSites();

    return c.real() / vol;
}

template <typename GImpl>
void TWilsonFlow<GImpl>::evolve_step(GaugeField &U, WilsonGaugeAction<GImpl> &SG) 
{
    GaugeField Z(U.Grid());
    GaugeField tmp(U.Grid());
    SG.deriv(U, Z);                                
    Z *= 0.25;                                          // Z0 = 1/4 * F(U)
    GImpl::update_field(Z, U, -2.0*par().step_size);    // U = W1 = exp(ep*Z0)*W0

    Z *= -17.0/8.0;
    SG.deriv(U, tmp); Z += tmp;                         // -17/32*Z0 + Z1
    Z *= 8.0/9.0;                                       // Z = -17/36*Z0 +8/9*Z1
    GImpl::update_field(Z, U, -2.0*par().step_size);    // U_= W2 = exp(ep*Z)*W1

    Z *= -4.0/3.0;
    SG.deriv(U, tmp); Z += tmp;                         // 4/3*(17/36*Z0 -8/9*Z1) + Z2
    Z *= 3.0/4.0;                                       // Z = 17/36*Z0 -8/9*Z1 +3/4*Z2
    GImpl::update_field(Z, U, -2.0*par().step_size);    // V(t+e) = exp(ep*Z)*W2
}

template <typename GImpl>
void TWilsonFlow<GImpl>::evolve_step_adaptive(GaugeField &U, RealD maxTau, RealD &epsilon, RealD &taus, WilsonGaugeAction<GImpl> &SG) 
{
    if (maxTau - taus < epsilon){
        epsilon = maxTau-taus;
    }
    //std::cout << GridLogMessage << "Integration epsilon : " << epsilon << std::endl;
    GaugeField Z(U.Grid());
    GaugeField Zprime(U.Grid());
    GaugeField tmp(U.Grid()), Uprime(U.Grid());
    Uprime = U;
    SG.deriv(U, Z);
    Zprime = -Z;
    Z *= 0.25;                                  // Z0 = 1/4 * F(U)
    GImpl::update_field(Z, U, -2.0*epsilon);    // U = W1 = exp(ep*Z0)*W0

    Z *= -17.0/8.0;
    SG.deriv(U, tmp); Z += tmp;                 // -17/32*Z0 +Z1
    Zprime += 2.0*tmp;
    Z *= 8.0/9.0;                               // Z = -17/36*Z0 +8/9*Z1
    GImpl::update_field(Z, U, -2.0*epsilon);    // U_= W2 = exp(ep*Z)*W1

    Z *= -4.0/3.0;
    SG.deriv(U, tmp); Z += tmp;                 // 4/3*(17/36*Z0 -8/9*Z1) +Z2
    Z *= 3.0/4.0;                               // Z = 17/36*Z0 -8/9*Z1 +3/4*Z2
    GImpl::update_field(Z, U, -2.0*epsilon);    // V(t+e) = exp(ep*Z)*W2      
    
    // Ramos
    GImpl::update_field(Zprime, Uprime, -2.0*epsilon); // V'(t+e) = exp(ep*Z')*W0
    // Compute distance as norm^2 of the difference
    GaugeField diffU = U - Uprime;
    RealD diff = norm2(diffU);   
    // adjust integration step  

    taus += epsilon;
    //std::cout << GridLogMessage << "Adjusting integration step with distance: " << diff << std::endl;

    epsilon = epsilon*0.95*std::pow(1e-4/diff,1./3.);
    //std::cout << GridLogMessage << "New epsilon : " << epsilon << std::endl;
}

template <typename GImpl>
void TWilsonFlow<GImpl>::status(double time, GaugeField &Umu, WilsonGaugeAction<GImpl> &SG)
{
    LOG(Message) << "flow time = " << std::setprecision(3) << std::fixed << time 
                 << " top. charge: " << std::setprecision(16) << std::scientific << WilsonLoops<GImpl>::TopologicalCharge(Umu)
                 << " plaquette: " << std::setprecision(16) << WilsonLoops<GImpl>::avgPlaquette(Umu) 
                 << " rectangle: " << std::setprecision(16) << WilsonLoops<GImpl>::avgRectangle(Umu) 
                 << " clover: " << std::setprecision(16) << avgClover(Umu)
                 << " action: " << std::setprecision(16) << SG.S(Umu) << std::endl;
}

// execution ///////////////////////////////////////////////////////////////////
template <typename GImpl>
void TWilsonFlow<GImpl>::execute(void)
{
    LOG(Message) << "Setting up Wilson Flow on '" << par().gauge << "' with " << par().steps
                 << " step" << ((par().steps > 1) ? "s." : ".") << std::endl;

    RealD mTau = -1.0;
    if(!par().maxTau.empty()) {
        LOG(Message) << "Using adaptive algorithm with maxTau = " << par().maxTau << std::endl;
        mTau = (RealD)std::stoi(par().maxTau);
    }
//    Grid::WilsonFlow<GImpl>  Wflow(par().steps, par().step_size, par().meas_interval);
    auto               &U   = envGet(GaugeField, par().gauge);
    auto               &Uwf = envGet(GaugeField, getName());

    envGetTmp(GaugeField, Umu);
    Umu = U;

    Uwf = U;
    double time = 0;
    WilsonGaugeAction<GImpl> SG(3.0);
    status(time,U,SG);
    if (mTau > 0) {
        RealD epsilon = par().step_size;
        RealD taus = par().step_size;
        unsigned int step = 0;
        do {
            step++;
            evolve_step_adaptive(Uwf, mTau, epsilon, taus, SG);
            if (step % par().meas_interval == 0) {
                status(taus,Uwf,SG);
            }
        } while (taus < mTau);
    } else {
        for (unsigned int step = 1; step <= par().steps; step++) {
            evolve_step(Uwf,SG);
            if (step % par().meas_interval == 0) {
                status(step*par().step_size,Uwf,SG);
            }
        }
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MGauge_WilsonFlow_hpp_
