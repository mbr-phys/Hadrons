/*
 * Utils_OBC.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
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
#ifndef Hadrons_MGradientFlow_Utils_OBC_hpp_
#define Hadrons_MGradientFlow_Utils_OBC_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/GaugeBoundary.hpp>
#include <cmath>
#include <limits>

BEGIN_HADRONS_NAMESPACE

BEGIN_MODULE_NAMESPACE(MGradientFlow)

namespace OBC
{

// Keep this schema in lockstep with MGradientFlow::GaugeResult in Utils.hpp.
// Utils_OBC is intentionally self-contained so it can replace Utils.hpp.
class GaugeResult : Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(GaugeResult,
                                    std::vector<double>,    flowtime,
                                    std::vector<double>,    plaquette,
                                    std::vector<double>,    rectangle,
                                    std::vector<double>,    clover,
                                    std::vector<double>,    topocharge,
                                    std::vector<double>,    action,
                                    std::vector<ComplexD>,  polyakovX,
                                    std::vector<ComplexD>,  polyakovY,
                                    std::vector<ComplexD>,  polyakovZ,
                                    std::vector<ComplexD>,  polyakovT,
                                    std::string,            boundaryType,
                                    int,                    timeDirection,
                                    double,                 cG,
                                    std::string,            observableScheme,
                                    bool,                   temporalPolyakovAvailable);
};

// Wilson flow implemented through the padded plaquette/rectangle force path.
// With c1 = 0 this is PlaqPlusRectangleAction(beta, 0).
template <class GImpl>
class WilsonAction : public RBCGaugeAction<GImpl> {
public:
    INHERIT_GIMPL_TYPES(GImpl);

    WilsonAction(RealD beta,
                 const OpenGaugeBoundaryParameters &boundaryParameters =
                     OpenGaugeBoundaryParameters())
    : RBCGaugeAction<GImpl>(beta, 0.0, boundaryParameters)
    {}

    virtual std::string action_name()
    {
        return "WilsonAction";
    }
};

// additional action(s) /////////////////////////////////////////////////////
template <class GImpl>
class ZeuthenGaugeAction {
private:
    typename WilsonLoops<GImpl>::StapleAndRectStapleAllWorkspace workspace;
public:
    INHERIT_GIMPL_TYPES(GImpl);
    
    double beta;
    OpenGaugeBoundaryParameters boundaryParameters;
    SymanzikGaugeAction<GImpl> SG;

    ZeuthenGaugeAction(
        double b,
        const OpenGaugeBoundaryParameters &boundaryParameters_ =
            OpenGaugeBoundaryParameters())
    : beta(b)
    , boundaryParameters(boundaryParameters_)
    , SG(SymanzikGaugeAction<GImpl>(b, boundaryParameters_))
    {
        if (boundaryParameters.open) {
            HADRONS_ERROR(Implementation,
                          "Zeuthen gauge flow with open boundaries is deferred: "
                          "the covariant improvement stencil is still periodic");
        }
    };

    virtual std::string action_name(){return "ZeuthenGaugeAction";}

    virtual double S(const GaugeField &U) {
        return SG.S(U);
    };

    virtual std::string LogParameters(){
        std::stringstream sstream;
        sstream << GridLogMessage << "["<<action_name() <<"] beta: " << beta << std::endl;
        sstream << SG.LogParameters();
        return sstream.str();
    }


    virtual void deriv(const GaugeField &Umu, GaugeField &dSdU) {
                                                //  beta = 3.0, cl = -1.0/12.0 -> Symanzik
        double factor_p = 5.0/double(Nc)*0.5;   //   5.0 = beta*(1.0-8.0*cl)
        double factor_r = -0.25/double(Nc)*0.5; // -0.25 = beta*cl

        GridBase *grid = Umu.Grid();

        std::vector<GaugeLinkField> U (Nd,grid);

        for(int mu=0;mu<Nd;mu++){
            U[mu] = PeekIndex<LorentzIndex>(Umu,mu);
        }

        std::vector<GaugeLinkField> RectStaple(Nd,grid), Staple(Nd,grid);
        WilsonLoops<GImpl>::StapleAndRectStapleAll(Staple, RectStaple, U, workspace);


        GaugeLinkField dSdU_mu(grid);
        GaugeLinkField tmq(grid),tmr(grid);

        for (int mu=0; mu < Nd; mu++){
            dSdU_mu = Ta(U[mu]*Staple[mu])*factor_p;
            dSdU_mu = dSdU_mu + Ta(U[mu]*RectStaple[mu])*factor_r;

            tmq = (adj(Cshift(U[mu],mu,-1)) * Cshift(dSdU_mu,mu,-1) * Cshift(U[mu],mu,-1));
            tmr = (U[mu] * Cshift(dSdU_mu,mu,1) * adj(U[mu]));

            dSdU_mu = 5.0/6.0*dSdU_mu + 1.0/12.0*tmq + 1.0/12.0*tmr;
            PokeIndex<LorentzIndex>(dSdU, dSdU_mu, mu);
        }
    };
};



// field evolution /////////////////////////////////////////////////////////////
template <typename FlowAction, typename GImpl, typename FImpl>
class Evolution {
    FERM_TYPE_ALIASES(FImpl,);
    typedef typename GImpl::GaugeLinkField GaugeLinkField;

    private:
        OpenGaugeBoundaryParameters boundaryParameters_;
        std::string gaugeBoundaryType_;
        double cG_;
        GridBase *grid_;
        GaugeLinkField linkBuf_;
        GaugeField zBuf1_, zBuf2_, uBuf1_, uBuf2_;
        LatticeComplex openTemporalMask_;
        TimerArray &timer_;
        std::vector<GaugeField> Wi_;

        static OpenGaugeBoundaryParameters makeBoundaryParameters(
            const std::string &gaugeBoundaryType, const double cG)
        {
            if (!GaugeBoundary::isValidType(gaugeBoundaryType) ||
                (gaugeBoundaryType == GaugeBoundary::unknown)) {
                HADRONS_ERROR(Argument,
                              "gauge boundary type must be 'periodic' or 'open'");
            }
            if (GaugeBoundary::isOpen(gaugeBoundaryType)) {
                if (!std::isfinite(cG) || cG < 0.0) {
                    HADRONS_ERROR(Argument,
                                  "cG must be finite and non-negative for open gauge boundaries");
                }
                return OpenGaugeBoundaryParameters::Open(Nd - 1, cG);
            }

            return OpenGaugeBoundaryParameters();
        }

        static FlowAction makeBetaAction(
            const double beta,
            const OpenGaugeBoundaryParameters &boundaryParameters)
        {
            return FlowAction(beta, boundaryParameters);
        }

        static FlowAction makeCoefficientAction(
            const double cPlaq, const double cRect,
            const OpenGaugeBoundaryParameters &boundaryParameters)
        {
            return FlowAction(cPlaq, cRect, boundaryParameters);
        }

        RealD loopWeightSum(
            const GaugeField &U,
            const std::vector<typename OpenGaugeBoundary<GImpl>::Loop> &loops) const
        {
            OpenGaugeBoundary<GImpl> boundary(boundaryParameters_);
            RealD total = 0.0;
            for (const auto &loop : loops) {
                if (loop.coefficient == 0.0) continue;
                auto weight = boundary.loopWeight(U.Grid(), loop.path);
                total += real(TensorRemove(sum(weight)));
            }
            return total;
        }

        double actionCompatibleLoopAverage(
            const GaugeField &U,
            const std::vector<typename OpenGaugeBoundary<GImpl>::Loop> &loops) const
        {
            OpenGaugeBoundary<GImpl> boundary(boundaryParameters_);
            const RealD weight = loopWeightSum(U, loops);
            GRID_ASSERT(weight > 0.0);
            return 1.0 - boundary.action(U, loops)/weight;
        }

        void maskTemporalEndpoints(GaugeLinkField &field) const
        {
            if (!GaugeBoundary::isOpen(gaugeBoundaryType_)) return;
            field = openTemporalMask_ * field;
        }

        void fieldStrengthOBC(GaugeLinkField &field, const GaugeField &U,
                              const int mu, const int nu) const
        {
            WilsonLoops<GImpl>::FieldStrength(field, U, mu, nu);
            maskTemporalEndpoints(field);
        }

        double topologicalChargeOBC(const GaugeField &U) const
        {
            GRID_ASSERT(Nd == 4);
            const int t = boundaryParameters_.timeDir;
            std::vector<int> spatial;
            for (int mu = 0; mu < Nd; ++mu) if (mu != t) spatial.push_back(mu);
            GRID_ASSERT(spatial.size() == 3);

            GaugeLinkField Bx(U.Grid()), By(U.Grid()), Bz(U.Grid());
            GaugeLinkField Ex(U.Grid()), Ey(U.Grid()), Ez(U.Grid());
            fieldStrengthOBC(Bx, U, spatial[1], spatial[2]);
            fieldStrengthOBC(By, U, spatial[2], spatial[0]);
            fieldStrengthOBC(Bz, U, spatial[0], spatial[1]);
            fieldStrengthOBC(Ex, U, t, spatial[0]);
            fieldStrengthOBC(Ey, U, t, spatial[1]);
            fieldStrengthOBC(Ez, U, t, spatial[2]);

            const double pi = std::acos(-1.0);
            ComplexField q(U.Grid());
            q = (8.0/(32.0*pi*pi))*trace(Bx*Ex + By*Ey + Bz*Ez);
            return TensorRemove(sum(q)).real();
        }

        double avgCloverOBC(const GaugeField &U)
        {
            ComplexField clov(U.Grid());
            GaugeLinkField field(U.Grid()), scaledUnit(U.Grid());
            clov = Zero();
            scaledUnit = 1.0/Nc;
            for (int mu = 1; mu < Nd; ++mu) {
                for (int nu = 0; nu < mu; ++nu) {
                    fieldStrengthOBC(field, U, mu, nu);
                    field -= trace(field)*scaledUnit;
                    clov -= trace(field*field);
                }
            }
            const int nt = U.Grid()->GlobalDimensions()[boundaryParameters_.timeDir];
            GRID_ASSERT(nt > 2);
            return TensorRemove(sum(clov)).real() /
                   (U.Grid()->gSites()*(double(nt - 2)/double(nt)));
        }

    public:
        double epsilon, maxTau, taus;
        FlowAction SG;

        // The open-boundary flow is defined on fields which vanish on the
        // two temporal endpoint slices.  Keep a scalar mask on the common
        // four-dimensional grid; it is used only when bc == 0.
        void setupOpenTemporalMask()
        {
            Lattice<iScalar<vInteger>> t(grid_);
            LatticeCoordinate(t, Nd - 1);

            LatticeComplex zero(grid_);
            zero = Zero();
            openTemporalMask_ = 1.0;
            openTemporalMask_ = where(t == 0, zero, openTemporalMask_);
            openTemporalMask_ = where(t == grid_->GlobalDimensions()[Nd - 1] - 1,
                                      zero, openTemporalMask_);
        }

        template <typename Field>
        void applyOpenTemporalBoundary(Field &field) const
        {
            field = openTemporalMask_ * field;
        }

        // constructor with beta
        Evolution(GridBase *grid, double beta, double step, double mTau,
                  double ts, const std::string &gaugeBoundaryType, double cG,
                  TimerArray *timer = nullptr)
        : boundaryParameters_(makeBoundaryParameters(gaugeBoundaryType, cG))
        , gaugeBoundaryType_(gaugeBoundaryType)
        , cG_(GaugeBoundary::isOpen(gaugeBoundaryType) ? cG : 1.0)
        , grid_(grid), linkBuf_(grid), zBuf1_(grid), zBuf2_(grid)
        , uBuf1_(grid), uBuf2_(grid), openTemporalMask_(grid), timer_(*timer)
        , epsilon(step), maxTau(mTau), taus(ts)
        , SG(makeBetaAction(beta, boundaryParameters_))
        {
            Wi_.resize(5, grid);
            setupOpenTemporalMask();
        }

        // constructor with c_plaq, c_rect
        Evolution(GridBase *grid, double c_plaq, double c_rect, double step,
                  double mTau, double ts,
                  const std::string &gaugeBoundaryType, double cG,
                  TimerArray *timer = nullptr)
        : boundaryParameters_(makeBoundaryParameters(gaugeBoundaryType, cG))
        , gaugeBoundaryType_(gaugeBoundaryType)
        , cG_(GaugeBoundary::isOpen(gaugeBoundaryType) ? cG : 1.0)
        , grid_(grid), linkBuf_(grid), zBuf1_(grid), zBuf2_(grid)
        , uBuf1_(grid), uBuf2_(grid), openTemporalMask_(grid), timer_(*timer)
        , epsilon(step), maxTau(mTau), taus(ts)
        , SG(makeCoefficientAction(c_plaq, c_rect, boundaryParameters_))
        {
            Wi_.resize(5, grid);
            setupOpenTemporalMask();
        }

        const OpenGaugeBoundaryParameters &gaugeBoundaryParameters(void) const
        {
            return boundaryParameters_;
        }

        const std::string &gaugeBoundaryType(void) const
        {
            return gaugeBoundaryType_;
        }

        double cG(void) const
        {
            return cG_;
        }

        // clover //////////////////////////////////////////////////////////////////////
        void siteClover(ComplexField &Clov, const GaugeField &U)
        {
            GaugeLinkField scaledUnit(U.Grid());
            Clov = Zero();
            scaledUnit = 1.0/Nc;
            for (int mu = 1; mu < Nd; mu++) {
                for (int nu = 0; nu < mu; nu++) {
                    linkBuf_ = PeekIndex<LorentzIndex>(U, mu);
                    WilsonLoops<GImpl>::FieldStrength(linkBuf_, U, mu, nu);
                    linkBuf_ -= trace(linkBuf_) * scaledUnit;
                    Clov -= trace(linkBuf_ * linkBuf_);
                }
            }
        }

        double avgClover(const GaugeField &Umu) 
        {
            ComplexField Clov(Umu.Grid());

            siteClover(Clov, Umu);
            auto Tc = sum(Clov);
            auto c = TensorRemove(Tc);

            double vol = Umu.Grid()->gSites();

            return c.real() / vol;
        }

        // polyakov loop in mu direction  //////////////////////////////////////////////
        ComplexD avgPolyakovLoopMu(const GaugeField &Umu, int mu) { // assuming Nd=4
            GaugeLinkField P(Umu.Grid());
            ComplexD out;

            double vol = Umu.Grid()->gSites();

            linkBuf_ = peekLorentz(Umu,mu);
            P = linkBuf_;
            for (int t=1; t < Umu.Grid()->GlobalDimensions()[mu]; t++) {
                P = GImpl::CovShiftForward(linkBuf_,mu,P);
            }
            RealD norm = 1.0/(Nc*vol);
            out = sum(trace(P))*norm;
            return out;
        }


        void gauge_RK(GaugeField &U) {
            Wi_[0] = U;                                     // W0
            timer_.startTimer("evolution_deriv");
            SG.deriv(U, zBuf1_);       
            timer_.stopTimer("evolution_deriv");                         
            zBuf1_ *= 0.25;                                 // Z0 = 1/4 * F(U)
            GImpl::update_field(zBuf1_, U, -2.0*epsilon);   // U = W1 = exp(ep*Z0)*W0
            Wi_[1] = U;                                     // W1

            zBuf1_ *= -17.0/8.0;
            timer_.startTimer("evolution_deriv");
            SG.deriv(U,uBuf1_); 
            timer_.stopTimer("evolution_deriv");
            zBuf1_ += uBuf1_;                                // -17/32*Z0 + Z1
            zBuf1_ *= 8.0/9.0;                               // Z = -17/36*Z0 +8/9*Z1
            GImpl::update_field(zBuf1_, U, -2.0*epsilon);    // U = W2 = exp(ep*Z)*W1
            Wi_[2] = U;                                      // W2

            zBuf1_ *= -4.0/3.0;
            timer_.startTimer("evolution_deriv");
            SG.deriv(U,uBuf1_); 
            timer_.stopTimer("evolution_deriv");
            zBuf1_ += uBuf1_;                                // 4/3*(17/36*Z0 -8/9*Z1) + Z2
            zBuf1_ *= 3.0/4.0;                               // Z = 17/36*Z0 -8/9*Z1 +3/4*Z2
            GImpl::update_field(zBuf1_, U, -2.0*epsilon);    // V(t+e) = exp(ep*Z)*W2
            Wi_[3] = U;                                      // W3
        };

        void gauge_RK_adaptive(GaugeField &U) {


            if (maxTau - taus < epsilon){
                epsilon = maxTau-taus;
            }
            uBuf2_ = U;
            Wi_[0] = U;                                        // W0
            SG.deriv(U, zBuf1_);
            zBuf2_ = -zBuf1_;
            zBuf1_ *= 0.25;                                    // Z0 = 1/4 * F(U)
            GImpl::update_field(zBuf1_, U, -2.0*epsilon);      // U = W1 = exp(ep*Z0)*W0
            Wi_[1] = U;                                        // W1

            zBuf1_ *= -17.0/8.0;
            SG.deriv(U, uBuf1_); zBuf1_ += uBuf1_;             // -17/32*Z0 +Z1
            zBuf2_ += 2.0*uBuf1_;
            zBuf1_ *= 8.0/9.0;                                 // Z = -17/36*Z0 +8/9*Z1
            GImpl::update_field(zBuf1_, U, -2.0*epsilon);      // U = W2 = exp(ep*Z)*W1
            Wi_[2] = U;                                        // W2

            zBuf1_ *= -4.0/3.0;
            SG.deriv(U, uBuf1_); zBuf1_ += uBuf1_;             // 4/3*(17/36*Z0 -8/9*Z1) +Z2
            zBuf1_ *= 3.0/4.0;                                 // Z = 17/36*Z0 -8/9*Z1 +3/4*Z2
            GImpl::update_field(zBuf1_, U, -2.0*epsilon);      // V(t+e) = exp(ep*Z)*W2      
            Wi_[3] = U;                                        // W3
            
            GImpl::update_field(zBuf2_, uBuf2_, -2.0*epsilon); // V'(t+e) = exp(ep*Z')*W0
            Wi_[4] = uBuf2_;                                   // Uprime
        };

        void adaptive_eps(const GaugeField& U, const GaugeField& Uprime) {
            // Compute distance as norm^2 of the difference
            uBuf1_ = U - Uprime;
            double diff = norm2(uBuf1_);   
            // adjust integration step  

            taus += epsilon;
            epsilon = epsilon*0.95*std::pow(1e-4/diff,1./3.);
        };

        std::vector<GaugeField> & evolve_gauge(GaugeField &U) {
            gauge_RK(U);
            U = Wi_[3];
            return Wi_;
        };
        
        void evolve_gauge_adaptive(GaugeField &U) {
            gauge_RK_adaptive(U);
            adaptive_eps(Wi_[3],Wi_[4]);
            U = Wi_[3];
        };

        
        void gauge_apply_boundary(GaugeField &Umu, std::vector<int> bc) {
            GaugeLinkField tmp1(Umu.Grid());
            GaugeLinkField tmp2(Umu.Grid());
            GaugeLinkField tmp3(Umu.Grid());
            Lattice<iScalar<vInteger>> coord(Umu.Grid());

            for (int mu = 0; mu < Nd; mu++) {
                LatticeCoordinate(coord,mu);
                
                tmp1 = PeekIndex<LorentzIndex>(Umu,mu);
                tmp2 = (double)bc[mu]*tmp1;
                int dimSize = Umu.Grid()->GlobalDimensions()[mu] - 1;
                tmp3 = where((coord == dimSize), tmp2, tmp1);
                PokeIndex<LorentzIndex>(Umu, tmp3, mu);
            }
        };

        template <typename Field>
        Field generic_laplace(double a, double b, GaugeField &Umu, const Field& x_in,
                               int skip_axis, int bc) {
            double Nx = Nd;
            if (skip_axis != -1) Nx--;

            Field x = x_in;
            if (bc == 0) {
                // Mask before applying periodic Grid shifts.  The shifted
                // values which would wrap t = 0 <-> T-1 are then zero.
                applyOpenTemporalBoundary(x);
            }

            Field x_out = (a + -2.0*Nx*b) * x;
            for (int mu = 0; mu < Nd; mu++) {
                if (mu != skip_axis) {
                    GaugeLinkField U = PeekIndex<LorentzIndex>(Umu, mu);
                    x_out += b*(GImpl::CovShiftForward(U,mu,x) + GImpl::CovShiftBackward(U,mu,x));
                }
            }
            if (bc == 0) {
                applyOpenTemporalBoundary(x_out);
            }
            return x_out;
        };

        template <typename Field>
        void laplace_flow(GaugeField &W0, GaugeField &W1, GaugeField &W2, Field &field, int bc) {
            if (bc == 0) {
                applyOpenTemporalBoundary(field);
            }

            Field psi1 = field + (epsilon/4.0)*generic_laplace(0.0, 1.0, W0, field, -1, bc);
            if (bc == 0) {
                applyOpenTemporalBoundary(psi1);
            }
            Field psi2 = field + (8.0*epsilon/9.0)*generic_laplace(0.0, 1.0, W1, psi1, -1, bc) - (2.0*epsilon/9.0)*generic_laplace(0.0, 1.0, W0, field, -1, bc);
            if (bc == 0) {
                applyOpenTemporalBoundary(psi2);
            }
            Field psi3 = psi1 + (3.0*epsilon/4.0)*generic_laplace(0.0, 1.0, W2, psi2, -1, bc);
            if (bc == 0) {
                applyOpenTemporalBoundary(psi3);
            }

            field = psi3;
        };

        // Adjoint (backward) fermion flow using reverse RK stepping
        // Implements: (delta_s + Delta_s) xi(tau;s) = 0, flowing from s+eps to s
        // Formula (reverse Hermitian adjoint of forward RK):
        //   lambda_3 = xi_{s+eps}
        //   lambda_2 = 3/4 Delta_2 lambda_3
        //   lambda_1 = lambda_3 + 8/9 Delta_1 lambda_2
        //   lambda_0 = lambda_1 + lambda_2 + 1/4 Delta_0(lambda_1 - 8/9 lambda_2)
        //   xi_s = lambda_0
        template <typename Field>
        void adjoint_laplace_flow(GaugeField &W0, GaugeField &W1, GaugeField &W2,
                                  Field &xi, int bc) {
            // lambda_3 = xi_{s+eps} (input)
            Field lambda3 = xi;
            if (bc == 0) {
                applyOpenTemporalBoundary(lambda3);
            }

            // lambda_2 = 3/4 Delta_2 lambda_3
            Field lambda2 = (3.0*epsilon/4.0)*generic_laplace(0.0, 1.0, W2, lambda3, -1, bc);
            if (bc == 0) {
                applyOpenTemporalBoundary(lambda2);
            }

            // lambda_1 = lambda_3 + 8/9 Delta_1 lambda_2
            Field lambda1 = lambda3 + (8.0*epsilon/9.0)*generic_laplace(0.0, 1.0, W1, lambda2, -1, bc);
            if (bc == 0) {
                applyOpenTemporalBoundary(lambda1);
            }

            // lambda_0 = lambda_1 + lambda_2+ 1/4 Delta_0(lambda_1 - 8/9 lambda_2)
            Field tmp = lambda1 - 8.0/9.0*lambda2;
            Field lambda0 = lambda1 + lambda2 + (epsilon/4.0)*generic_laplace(0.0, 1.0, W0, tmp, -1, bc);
            if (bc == 0) {
                applyOpenTemporalBoundary(lambda0);
            }

            // xi_s = lambda_0 (output)
            xi = lambda0;
        };

        std::vector<GaugeField> & evolve_gaugeFF(GaugeField &U, std::vector<int> &bc) {
            gauge_RK(U);
            U = 1.0*Wi_[3];

            gauge_apply_boundary(Wi_[0],bc);
            gauge_apply_boundary(Wi_[1],bc);
            gauge_apply_boundary(Wi_[2],bc);

            return Wi_;
        };

        // gauge field status //////////////////////////////////////////////////////////
        void gauge_status(GaugeField &Umu, GaugeResult &result, double flowt)
        {
            const bool open = GaugeBoundary::isOpen(gaugeBoundaryType_);
            result.boundaryType = gaugeBoundaryType_;
            result.timeDirection = boundaryParameters_.timeDir;
            result.cG = cG_;
            result.observableScheme = open ? "obc-action-compatible-v1" : "periodic-v1";
            result.temporalPolyakovAvailable = !open;

            double Q = open ? topologicalChargeOBC(Umu)
                            : WilsonLoops<GImpl>::TopologicalCharge(Umu);
            double plaq;
            double rect;
            double clov;
            if (open) {
                using Boundary = OpenGaugeBoundary<GImpl>;
                plaq = actionCompatibleLoopAverage(Umu, Boundary::plaquetteLoops(1.0));
                rect = actionCompatibleLoopAverage(
                    Umu, Boundary::plaquetteRectangleLoops(0.0, 1.0));
                clov = avgCloverOBC(Umu);
            } else {
                plaq = WilsonLoops<GImpl>::avgPlaquette(Umu);
                rect = WilsonLoops<GImpl>::avgRectangle(Umu);
                clov = avgClover(Umu);
            }
            double act = SG.S(Umu);
            ComplexD polyX = avgPolyakovLoopMu(Umu,0);
            ComplexD polyY = avgPolyakovLoopMu(Umu,1);
            ComplexD polyZ = avgPolyakovLoopMu(Umu,2);
            ComplexD polyT = open
                ? ComplexD(std::numeric_limits<RealD>::quiet_NaN(),
                           std::numeric_limits<RealD>::quiet_NaN())
                : avgPolyakovLoopMu(Umu,3);

            result.flowtime.push_back(flowt);
            result.plaquette.push_back(plaq);
            result.rectangle.push_back(rect);
            result.clover.push_back(clov);
            result.topocharge.push_back(Q);
            result.action.push_back(act);
            result.polyakovX.push_back(polyX);
            result.polyakovY.push_back(polyY);
            result.polyakovZ.push_back(polyZ);
            result.polyakovT.push_back(polyT);
        };
};

} // namespace OBC

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MGradientFlow_Utils_OBC_hpp_
