/*
 * GradientFlow.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
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
#ifndef Hadrons_MGauge_GradientFlow_hpp_
#define Hadrons_MGauge_GradientFlow_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                               Gradient Flow                                  *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MGauge)

template <class GImpl>
class ZeuthenGaugeAction {
public:
    INHERIT_GIMPL_TYPES(GImpl);
    
    RealD beta;
    SymanzikGaugeAction<GImpl> SG;

    ZeuthenGaugeAction(RealD b): beta(b),SG(SymanzikGaugeAction<GImpl>(b)) {};

    virtual std::string action_name(){return "ZeuthenGaugeAction";}

    virtual RealD S(const GaugeField &U) {
        return SG.S(U);
    };

    virtual void deriv(const GaugeField &Umu, GaugeField &dSdU) {
                                              //  beta = 3.0, cl = -1.0/12.0 -> Symanzik
        RealD factor_p = 5.0/RealD(Nc)*0.5;   //   5.0 = beta*(1.0-8.0*cl)
        RealD factor_r = -0.25/RealD(Nc)*0.5; // -0.25 = beta*cl

        GridBase *grid = Umu.Grid();

        std::vector<GaugeLinkField> U (Nd,grid);
        std::vector<GaugeLinkField> U2(Nd,grid);

        for(int mu=0;mu<Nd;mu++){
            U[mu] = PeekIndex<LorentzIndex>(Umu,mu);
            WilsonLoops<GImpl>::RectStapleDouble(U2[mu],U[mu],mu);
        }

        GaugeLinkField dSdU_mu(grid);
        GaugeLinkField staple(grid);
        GaugeLinkField tmp(grid),tmq(grid),tmr(grid);

        for (int mu=0; mu < Nd; mu++){
            // Staple in direction mu
            WilsonLoops<GImpl>::Staple(staple,Umu,mu);
            tmp = Ta(U[mu]*staple)*factor_p;

            WilsonLoops<GImpl>::RectStaple(Umu,staple,U2,U,mu);
            tmp = tmp + Ta(U[mu]*staple)*factor_r;

            tmq = (adj(Cshift(U[mu],mu,-1)) * Cshift(tmp,mu,-1) * Cshift(U[mu],mu,-1));
            tmr = (U[mu] * Cshift(tmp,mu,1) * adj(U[mu]));

            dSdU_mu = 5.0/6.0*tmp + 1.0/12.0*tmq + 1.0/12.0*tmr;
            PokeIndex<LorentzIndex>(dSdU, dSdU_mu, mu);
        }
    };
};

class GradientFlowPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(GradientFlowPar,
                                    std::string, gauge,
                                    int, steps,
                                    double, step_size,
                                    int, meas_interval,
                                    std::string, maxTau); 
};

template <typename GImpl,typename FlowAction>
class TGradientFlow: public Module<GradientFlowPar>
{
public:
    GAUGE_TYPE_ALIASES(GImpl,);
public:
    // constructor
    TGradientFlow(const std::string name);
    // destructor
    virtual ~TGradientFlow(void) {};
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
    // action
    FlowAction SG = FlowAction(3.0);
    // clover
    void siteClover(ComplexField &Clov, const GaugeLorentz &U);
    RealD avgClover(const GaugeLorentz &Umu);
    // gauge observable measurements at flow time
    void status(double time, GaugeField &Umu);
    // gauge field evolutions
    void evolve_step(GaugeField &U);
    void evolve_step_adaptive(GaugeField &U, RealD maxTau, RealD &epsilon, RealD &taus);
};

MODULE_REGISTER_TMP(WilsonFlow, ARG(TGradientFlow<GIMPL,WilsonGaugeAction<GIMPL>>), MGauge);
MODULE_REGISTER_TMP(SymanzikFlow, ARG(TGradientFlow<GIMPL,SymanzikGaugeAction<GIMPL>>), MGauge);
MODULE_REGISTER_TMP(ZeuthenFlow, ARG(TGradientFlow<GIMPL,ZeuthenGaugeAction<GIMPL>>), MGauge);

/******************************************************************************
 *                     TGradientFlow implementation                          *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename GImpl,typename FlowAction>
TGradientFlow<GImpl,FlowAction>::TGradientFlow(const std::string name)
: Module<GradientFlowPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename GImpl,typename FlowAction>
std::vector<std::string> TGradientFlow<GImpl,FlowAction>::getInput(void)
{
    std::vector<std::string> in = {par().gauge};
    
    return in;
}

template <typename GImpl,typename FlowAction>
std::vector<std::string> TGradientFlow<GImpl,FlowAction>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename GImpl,typename FlowAction>
void TGradientFlow<GImpl,FlowAction>::setup(void)
{
    envCreateLat(GaugeField, getName());
    envTmpLat(GaugeField, "Umu");
}

// clover //////////////////////////////////////////////////////////////////////
template <typename GImpl,typename FlowAction>
void TGradientFlow<GImpl,FlowAction>::siteClover(ComplexField &Clov, const GaugeLorentz &U)
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

template <typename GImpl,typename FlowAction>
RealD TGradientFlow<GImpl,FlowAction>::avgClover(const GaugeLorentz &Umu) 
{
    ComplexField Clov(Umu.Grid());

    siteClover(Clov, Umu);
    auto Tc = sum(Clov);
    auto c = TensorRemove(Tc);

    double vol = Umu.Grid()->gSites();

    return c.real() / vol;
}

template <typename GImpl,typename FlowAction>
void TGradientFlow<GImpl,FlowAction>::evolve_step(GaugeField &U) 
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

template <typename GImpl,typename FlowAction>
void TGradientFlow<GImpl,FlowAction>::evolve_step_adaptive(GaugeField &U, RealD maxTau, RealD &epsilon, RealD &taus) 
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

template <typename GImpl,typename FlowAction>
void TGradientFlow<GImpl,FlowAction>::status(double time, GaugeField &Umu)
{
    LOG(Message) << "flow time = " << std::setprecision(3) << std::fixed << time 
                 << " top. charge: " << std::setprecision(16) << std::scientific << WilsonLoops<GImpl>::TopologicalCharge(Umu)
                 << " plaquette: " << WilsonLoops<GImpl>::avgPlaquette(Umu) 
                 << " rectangle: " << WilsonLoops<GImpl>::avgRectangle(Umu) 
                 << " clover: " << avgClover(Umu)
                 << " action: " << SG.S(Umu) << std::endl;
}

// execution ///////////////////////////////////////////////////////////////////
template <typename GImpl,typename FlowAction>
void TGradientFlow<GImpl,FlowAction>::execute(void)
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
    auto               &U   = envGet(GaugeField, par().gauge);
    auto               &Uwf = envGet(GaugeField, getName());

    envGetTmp(GaugeField, Umu);
    Umu = U;

    Uwf = U;
    double time = 0;
    status(time,U);
    if (mTau > 0) {
        RealD epsilon = par().step_size;
        RealD taus = par().step_size;
        unsigned int step = 0;
        do {
            step++;
            evolve_step_adaptive(Uwf, mTau, epsilon, taus);
            if (step % par().meas_interval == 0) {
                status(taus,Uwf);
            }
        } while (taus < mTau);
    } else {
        for (unsigned int step = 1; step <= par().steps; step++) {
            evolve_step(Uwf);
            if (step % par().meas_interval == 0) {
                status(step*par().step_size,Uwf);
            }
        }
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MGauge_GradientFlow_hpp_
