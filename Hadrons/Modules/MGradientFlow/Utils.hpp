/*
 * Utils.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
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
#ifndef Hadrons_MGradientFlow_Utils_hpp_
#define Hadrons_MGradientFlow_Utils_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>

BEGIN_HADRONS_NAMESPACE

BEGIN_MODULE_NAMESPACE(MGradientFlow)

// additional action(s) /////////////////////////////////////////////////////
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

// clover //////////////////////////////////////////////////////////////////////
template <typename GImpl, typename ComplexField, typename GaugeLorentz, typename GaugeMat>
void siteClover(ComplexField &Clov, const GaugeLorentz &U)
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

template <typename GImpl, typename ComplexField, typename GaugeLorentz, typename GaugeMat>
RealD avgClover(const GaugeLorentz &Umu) 
{
    ComplexField Clov(Umu.Grid());

    siteClover<GImpl,ComplexField,GaugeLorentz,GaugeMat>(Clov, Umu);
    auto Tc = sum(Clov);
    auto c = TensorRemove(Tc);

    double vol = Umu.Grid()->gSites();

    return c.real() / vol;
}

// field evolution /////////////////////////////////////////////////////////
template <typename GImpl, typename FlowAction>
class Evolution {
    public:
        INHERIT_GIMPL_TYPES(GImpl);

        RealD epsilon, maxTau, taus;
        FlowAction SG;
        Evolution(RealD beta, RealD step, RealD mTau, RealD ts) : 
            SG(FlowAction(beta)), epsilon(step), maxTau(mTau), taus(ts) {};

        virtual void evolve_gauge(GaugeField &U) {
            GaugeField Z(U.Grid());
            GaugeField tmp(U.Grid());
            SG.deriv(U, Z);                                
            Z *= 0.25;                                  // Z0 = 1/4 * F(U)
            GImpl::update_field(Z, U, -2.0*epsilon);    // U = W1 = exp(ep*Z0)*W0

            Z *= -17.0/8.0;
            SG.deriv(U, tmp); Z += tmp;                 // -17/32*Z0 + Z1
            Z *= 8.0/9.0;                               // Z = -17/36*Z0 +8/9*Z1
            GImpl::update_field(Z, U, -2.0*epsilon);    // U_= W2 = exp(ep*Z)*W1

            Z *= -4.0/3.0;
            SG.deriv(U, tmp); Z += tmp;                 // 4/3*(17/36*Z0 -8/9*Z1) + Z2
            Z *= 3.0/4.0;                               // Z = 17/36*Z0 -8/9*Z1 +3/4*Z2
            GImpl::update_field(Z, U, -2.0*epsilon);    // V(t+e) = exp(ep*Z)*W2
        };

        virtual void evolve_gauge_adaptive(GaugeField &U) {
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
        };

        template <typename FImpl>
        FImpl generic_laplace(RealD a, RealD b, GaugeField &Umu, const FImpl& x_in, int skip_axis) {
            std::vector<GaugeLinkField> U(Nd,Umu.Grid());
            for (int mu = 0; mu < Nd; mu++) {
                U[mu] = PeekIndex<LorentzIndex>(Umu,mu);
            }

            FImpl tmp = x_in;
            FImpl x_out = x_in;
            x_out *= a;
            for (int mu = 0; mu < Nd; mu++) {
                if (mu != skip_axis) {
                    x_out += b * (-2.0*tmp + U[mu]*Cshift(tmp,mu,1) + Cshift(adj(U[mu])*tmp,mu,-1));
                }
            }
            return x_out;
        };

        template <typename FImpl>
        void laplace_flow(GaugeField &W0, GaugeField &W1, GaugeField &W2, FImpl &prop) {
            // is it -2*epsilon or just epsilon here?
            FImpl psi1 = prop + (-2.0*epsilon/4.0)*generic_laplace<FImpl>(0.0, 1.0, W0, prop, -1);
            FImpl psi2 = prop + 8.0*(-2.0*epsilon/9.0)*generic_laplace<FImpl>(0.0, 1.0, W1, psi1, -1) - 2.0*(-2.0*epsilon/9.0)*generic_laplace(0.0, 1.0, W0, prop, -1);
            FImpl psi3 = psi1 + 3.0*(-2.0*epsilon/4.0)*generic_laplace<FImpl>(0.0, 1.0, W2, psi2, -1);

            prop = psi3;
        };

        template <typename FImpl1, typename FImpl2>
        void evolve_prop(GaugeField &U, FImpl1 &q1, FImpl2 &q2) {
            GaugeField Z(U.Grid());
            GaugeField tmp(U.Grid());
            GaugeField W0(U.Grid()),W1(U.Grid()),W2(U.Grid());//W3(U.Grid());
            W0 = U;
            SG.deriv(U, Z);                                
            Z *= 0.25;                                  // Z0 = 1/4 * F(U)
            GImpl::update_field(Z, U, -2.0*epsilon);    // U = W1 = exp(ep*Z0)*W0
            W1 = U;

            Z *= -17.0/8.0;
            SG.deriv(U, tmp); Z += tmp;                 // -17/32*Z0 + Z1
            Z *= 8.0/9.0;                               // Z = -17/36*Z0 +8/9*Z1
            GImpl::update_field(Z, U, -2.0*epsilon);    // U = W2 = exp(ep*Z)*W1
            W2 = U;

            Z *= -4.0/3.0;
            SG.deriv(U, tmp); Z += tmp;                 // 4/3*(17/36*Z0 -8/9*Z1) + Z2
            Z *= 3.0/4.0;                               // Z = 17/36*Z0 -8/9*Z1 +3/4*Z2
            GImpl::update_field(Z, U, -2.0*epsilon);    // V(t+e) = exp(ep*Z)*W2
            //W3 = U;

            // gauge boundary conditions in here?

            laplace_flow<FImpl1>(W0,W1,W2,q1);
            laplace_flow<FImpl2>(W0,W1,W2,q2);
        };

// is this something to add, and if so, does prop evolution change?
//        virtual void evolve_prop_adaptive(GaugeField &U) {
//        };
};

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MGradientFlow_Utils_hpp_
