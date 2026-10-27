/*
 * HQEOperators.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2026
 *
 * Author: Matthew Black <matthewkblack@protonmail.com>
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
#ifndef Hadrons_MContraction_HQEOperators_hpp_
#define Hadrons_MContraction_HQEOperators_hpp_

#include <Hadrons/Global.hpp>

BEGIN_HADRONS_NAMESPACE

/*
 * Lattice building blocks for HQE / heavy-meson operator insertions.
 * ------------------------------------------------------------------
 * Shared by MUtilities::CovariantDerivative, MGauge::FieldStrength,
 * MUtilities::GaugeInsertion and MContraction::HeavyMeson3pt, so that all of
 * them use the same conventions:
 *
 *  - covariant derivatives on fields in the fundamental (Euclidean, lattice
 *    units, acting to the right):
 *      symmetric: D f(x) = 1/2 [U_mu(x) f(x+mu) - U_mu^dag(x-mu) f(x-mu)]
 *      forward:   D f(x) = U_mu(x) f(x+mu) - f(x)
 *      backward:  D f(x) = f(x) - U_mu^dag(x-mu) f(x-mu)
 *    For mu = Tp the link on the last time slice is multiplied by the fermion
 *    temporal boundary phase bc (+1 periodic, -1 antiperiodic).
 *  - Laplacian: sum_{mu in dirs} [U_mu(x) f(x+mu) + U_mu^dag(x-mu) f(x-mu) - 2 f(x)]
 *  - field strength: Grid clover WilsonLoops::FieldStrength, projected to be
 *    traceless. It is anti-hermitian, F_{mu nu} ~ i a^2 g G_{mu nu} up to the
 *    sign convention for the gauge potential. Stored for mu < nu with
 *    F_{nu mu} = -F_{mu nu}.
 *  - adjoint covariant derivative of a colour-matrix field (symmetric):
 *      D M(x) = 1/2 [U_mu(x) M(x+mu) U_mu^dag(x) - U_mu^dag(x-mu) M(x-mu) U_mu(x-mu)]
 *  - sigma_{mu nu} = (i/2)[gamma_mu, gamma_nu] = i gamma_mu gamma_nu (mu != nu),
 *    with Grid's Euclidean gamma matrices.
 *  - chromomagnetic: sigma.B f = sum_{i<j<Tp} sigma_ij F_ij f
 *    chromoelectric: sigma.E f = sum_{j<Tp}   sigma_jt F_jt f
 */

namespace HQE
{
    enum class DerivType {symmetric, forward, backward};

    inline DerivType derivType(const std::string &s)
    {
        if (s.empty() or (s == "symmetric"))
        {
            return DerivType::symmetric;
        }
        else if (s == "forward")
        {
            return DerivType::forward;
        }
        else if (s == "backward")
        {
            return DerivType::backward;
        }
        HADRONS_ERROR(Argument, "unknown derivative type '" + s
                      + "' (expected symmetric, forward or backward)");

        return DerivType::symmetric;
    }

    // storage index of F_{mu nu}, symmetric in (mu, nu), -1 on the diagonal
    inline int fsIndex(const int mu, const int nu)
    {
        static const int idx[4][4] = {{-1,  0,  1,  2},
                                      { 0, -1,  3,  4},
                                      { 1,  3, -1,  5},
                                      { 2,  4,  5, -1}};

        return idx[mu][nu];
    }

    // F_{mu nu} = fsSign(mu, nu)*F[fsIndex(mu, nu)]
    inline double fsSign(const int mu, const int nu)
    {
        return (mu < nu) ? 1. : -1.;
    }

    // gauge link in direction mu, with the fermion temporal boundary phase
    template <typename GImpl>
    void boundaryLink(typename GImpl::GaugeLinkField &link,
                      const typename GImpl::GaugeField &U, const int mu,
                      const int bc)
    {
        link = PeekIndex<LorentzIndex>(U, mu);
        if ((mu == Tp) and (bc != 1))
        {
            typename GImpl::GaugeLinkField flipped(U.Grid());
            Lattice<iScalar<vInteger>>     t(U.Grid());
            const unsigned int             tLast = U.Grid()->GlobalDimensions()[Tp] - 1;

            LatticeCoordinate(t, Tp);
            flipped = static_cast<RealD>(bc)*link;
            link    = where(t >= tLast, flipped, link);
        }
    }

    // covariant derivative of a fundamental field (link must include any bc)
    template <typename GImpl, typename Field>
    void covDeriv(Field &out, const Field &in,
                  const typename GImpl::GaugeLinkField &link, const int mu,
                  const DerivType type)
    {
        switch (type)
        {
            case DerivType::symmetric:
                out = 0.5*(GImpl::CovShiftForward(link, mu, in)
                           - GImpl::CovShiftBackward(link, mu, in));
                break;
            case DerivType::forward:
                out = GImpl::CovShiftForward(link, mu, in) - in;
                break;
            case DerivType::backward:
                out = in - GImpl::CovShiftBackward(link, mu, in);
                break;
        }
    }

    // nearest-neighbour covariant Laplacian over the directions in dirs
    template <typename GImpl, typename Field>
    void laplacian(Field &out, const Field &in,
                   const typename GImpl::GaugeField &U,
                   const std::vector<int> &dirs, const int bc)
    {
        typename GImpl::GaugeLinkField link(U.Grid());

        out = Zero();
        for (int mu: dirs)
        {
            boundaryLink<GImpl>(link, U, mu, bc);
            out += GImpl::CovShiftForward(link, mu, in)
                   + GImpl::CovShiftBackward(link, mu, in) - 2.*in;
        }
    }

    // traceless clover field strength F_{mu nu}
    template <typename GImpl>
    void fieldStrength(typename GImpl::GaugeLinkField &F,
                       const typename GImpl::GaugeField &U,
                       const int mu, const int nu)
    {
        typename GImpl::GaugeLinkField unit(U.Grid());

        WilsonLoops<GImpl>::FieldStrength(F, U, mu, nu);
        unit = 1./Nc;
        F    = F - trace(F)*unit;
    }

    // all six components, F[fsIndex(mu, nu)] = F_{mu nu} for mu < nu
    template <typename GImpl>
    void allFieldStrengths(std::vector<typename GImpl::GaugeLinkField> &F,
                           const typename GImpl::GaugeField &U)
    {
        F.clear();
        F.resize(6, typename GImpl::GaugeLinkField(U.Grid()));
        for (int mu = 0; mu < Nd; ++mu)
        for (int nu = mu + 1; nu < Nd; ++nu)
        {
            fieldStrength<GImpl>(F[fsIndex(mu, nu)], U, mu, nu);
        }
    }

    // symmetric covariant derivative of an adjoint (colour-matrix) field
    template <typename GImpl>
    void adjCovDeriv(typename GImpl::GaugeLinkField &out,
                     const typename GImpl::GaugeLinkField &M,
                     const typename GImpl::GaugeField &U, const int mu)
    {
        typename GImpl::GaugeLinkField link(U.Grid()), tmp(U.Grid());

        link = PeekIndex<LorentzIndex>(U, mu);
        tmp  = adj(link)*M*link;
        out  = 0.5*(link*Cshift(M, mu, 1)*adj(link) - Cshift(tmp, mu, -1));
    }

    // sum_{mu in dirs, mu != nu} [D_mu, F_{mu nu}]
    template <typename GImpl>
    void fieldStrengthDivergence(typename GImpl::GaugeLinkField &out,
                                 const std::vector<typename GImpl::GaugeLinkField> &F,
                                 const typename GImpl::GaugeField &U,
                                 const int nu, const std::vector<int> &dirs)
    {
        typename GImpl::GaugeLinkField tmp(U.Grid());

        out = Zero();
        for (int mu: dirs)
        {
            if (mu == nu)
            {
                continue;
            }
            adjCovDeriv<GImpl>(tmp, F[fsIndex(mu, nu)], U, mu);
            out += fsSign(mu, nu)*tmp;
        }
    }

    // out = sigma_{mu nu} in, mu != nu
    template <typename Field>
    void applySigma(Field &out, const Field &in, const int mu, const int nu)
    {
        out = Gamma::gmu[nu]*in;
        out = Gamma::gmu[mu]*out;
        out = timesI(out);
    }

    // out = sum_{i<j<Tp} sigma_ij F_ij in
    template <typename GImpl, typename Field>
    void applySigmaB(Field &out, const Field &in,
                     const std::vector<typename GImpl::GaugeLinkField> &F)
    {
        Field tmp(in.Grid());

        out = Zero();
        for (int i = 0; i < Tp; ++i)
        for (int j = i + 1; j < Tp; ++j)
        {
            applySigma(tmp, in, i, j);
            out += F[fsIndex(i, j)]*tmp;
        }
    }

    // out = sum_{j<Tp} sigma_jt F_jt in
    template <typename GImpl, typename Field>
    void applySigmaE(Field &out, const Field &in,
                     const std::vector<typename GImpl::GaugeLinkField> &F)
    {
        Field tmp(in.Grid());

        out = Zero();
        for (int j = 0; j < Tp; ++j)
        {
            applySigma(tmp, in, j, Tp);
            out += F[fsIndex(j, Tp)]*tmp;
        }
    }
}

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_HQEOperators_hpp_
