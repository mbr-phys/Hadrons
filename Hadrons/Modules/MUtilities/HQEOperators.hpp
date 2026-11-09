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
#ifndef Hadrons_MUtilities_HQEOperators_hpp_
#define Hadrons_MUtilities_HQEOperators_hpp_

#include <Hadrons/Global.hpp>

BEGIN_HADRONS_NAMESPACE

/*
 * Lattice building blocks for HQE / heavy-meson operator insertions.
 * ------------------------------------------------------------------
 * Shared by MGauge::FieldStrength, MUtilities::GaugeInsertion and MContraction::HeavyMeson3pt, so that all of them use the same conventions:
 *
 *  - covariant derivatives on fundamental fields: see MUtilities/DerivHelpers.hpp
 *  - field strength: Grid clover WilsonLoops::FieldStrength, projected to be traceless. It is anti-hermitian, F_{mu nu} ~ i a^2 g G_{mu nu} up to the 
 *    sign convention for the gauge potential. The clover is exactly antisymmetric, F_{nu mu} = -F_{mu nu}.
 *  - chromoelectric and chromomagnetic components, each stored as 3 colour-matrix fields:
 *          E[j] = F_{j t}                                  (j = 0, 1, 2)
 *          B[k] = F_{ij}, (i, j, k) cyclic permutation of (0, 1, 2), i.e. B = (F_12, F_20, F_01)
 *  - adjoint covariant derivative of a colour-matrix field (symmetric):
 *          D M(x) = 1/2 [U_mu(x) M(x+mu) U_mu^dag(x) - U_mu^dag(x-mu) M(x-mu) U_mu(x-mu)]
 *  - sigma_{mu nu} = (i/2)[gamma_mu, gamma_nu] = i gamma_mu gamma_nu (mu != nu), with Grid's Euclidean gamma matrices.
 *  - chromomagnetic: sigma.B f = sum_{i<j<Tp} sigma_ij F_ij f = sum_k sigma_ij B[k] f  (cyclic i, j, k)
 *    chromoelectric: sigma.E f = sum_{j<Tp}   sigma_jt F_jt f = sum_j sigma_jt E[j] f
 *  - divergence of E: div E = sum_{i<Tp} [D_i, E[i]]
 *  - curl of E, stored like B: curlE[k] = [D_i, E[j]] - [D_j, E[i]], (i, j, k) cyclic, so that
 *    sum_{i != j < Tp} sigma_ij [D_i, E[j]] f = applySigmaB(f, curlE)
 */

namespace HQE
{
    // the pair (i, j) with (i, j, k) a cyclic permutation of (0, 1, 2)
    inline int cyclicI(const int k) { return (k + 1) % 3; }
    inline int cyclicJ(const int k) { return (k + 2) % 3; }

    // traceless clover field strength F_{mu nu}
    template <typename GImpl>
    void fieldStrength(typename GImpl::GaugeLinkField &F, const typename GImpl::GaugeField &U, const int mu, const int nu)
    {
        typename GImpl::GaugeLinkField unit(U.Grid());

        WilsonLoops<GImpl>::FieldStrength(F, U, mu, nu);
        unit = 1./Nc;
        F = F - trace(F)*unit;
    }

    // E[j] = F_{j t}
    template <typename GImpl>
    void chromoE(std::vector<typename GImpl::GaugeLinkField> &E, const typename GImpl::GaugeField &U)
    {
        E.resize(3, typename GImpl::GaugeLinkField(U.Grid()));
        for (int j = 0; j < Tp; ++j)
        {
            fieldStrength<GImpl>(E[j], U, j, Tp);
        }
    }

    // B[k] = F_{ij}, (i, j, k) cyclic
    template <typename GImpl>
    void chromoB(std::vector<typename GImpl::GaugeLinkField> &B, const typename GImpl::GaugeField &U)
    {
        B.resize(3, typename GImpl::GaugeLinkField(U.Grid()));
        for (int k = 0; k < Tp; ++k)
        {
            fieldStrength<GImpl>(B[k], U, cyclicI(k), cyclicJ(k));
        }
    }

    // symmetric covariant derivative of an adjoint (colour-matrix) field
    template <typename GImpl>
    void adjCovDeriv(typename GImpl::GaugeLinkField &out, const typename GImpl::GaugeLinkField &M, const typename GImpl::GaugeField &U, const int mu)
    {
        typename GImpl::GaugeLinkField link(U.Grid()), tmp(U.Grid());

        link = PeekIndex<LorentzIndex>(U, mu);
        tmp  = adj(link)*M*link;
        out  = 0.5*(link*Cshift(M, mu, 1)*adj(link) - Cshift(tmp, mu, -1));
    }

    // div E = sum_i [D_i, E[i]]
    template <typename GImpl>
    void divE(typename GImpl::GaugeLinkField &out, const std::vector<typename GImpl::GaugeLinkField> &E, const typename GImpl::GaugeField &U)
    {
        typename GImpl::GaugeLinkField tmp(U.Grid());

        out = Zero();
        for (int i = 0; i < Tp; ++i)
        {
            adjCovDeriv<GImpl>(tmp, E[i], U, i);
            out += tmp;
        }
    }

    // curlE[k] = [D_i, E[j]] - [D_j, E[i]], (i, j, k) cyclic
    template <typename GImpl>
    void curlE(std::vector<typename GImpl::GaugeLinkField> &out, const std::vector<typename GImpl::GaugeLinkField> &E, const typename GImpl::GaugeField &U)
    {
        typename GImpl::GaugeLinkField tmp(U.Grid());

        out.resize(3, typename GImpl::GaugeLinkField(U.Grid()));
        for (int k = 0; k < Tp; ++k)
        {
            adjCovDeriv<GImpl>(out[k], E[cyclicJ(k)], U, cyclicI(k));
            adjCovDeriv<GImpl>(tmp, E[cyclicI(k)], U, cyclicJ(k));
            out[k] -= tmp;
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

    // out = sum_k sigma_ij B[k] in, (i, j, k) cyclic
    template <typename GImpl, typename Field>
    void applySigmaB(Field &out, const Field &in, const std::vector<typename GImpl::GaugeLinkField> &B)
    {
        Field tmp(in.Grid());

        out = Zero();
        for (int k = 0; k < Tp; ++k)
        {
            applySigma(tmp, in, cyclicI(k), cyclicJ(k));
            out += B[k]*tmp;
        }
    }

    // out = sum_j sigma_jt E[j] in
    template <typename GImpl, typename Field>
    void applySigmaE(Field &out, const Field &in, const std::vector<typename GImpl::GaugeLinkField> &E)
    {
        Field tmp(in.Grid());

        out = Zero();
        for (int j = 0; j < Tp; ++j)
        {
            applySigma(tmp, in, j, Tp);
            out += E[j]*tmp;
        }
    }
}

END_HADRONS_NAMESPACE

#endif // Hadrons_MUtilities_HQEOperators_hpp_
