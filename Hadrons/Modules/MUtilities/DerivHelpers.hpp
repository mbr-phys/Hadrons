/*
 * DerivHelpers.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
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
#ifndef Hadrons_MUtilities_DerivHelpers_hpp_
#define Hadrons_MUtilities_DerivHelpers_hpp_

#include <Hadrons/Global.hpp>

BEGIN_HADRONS_NAMESPACE

/*
 * Covariant derivatives of fields in the fundamental
 * --------------------------------------------------
 * Shared by MUtilities::CovariantDerivative and MContraction::HeavyMeson3pt.
 * Euclidean, lattice units, acting to the right:
 *          symmetric: D f(x) = 1/2 [U_mu(x) f(x+mu) - U_mu^dag(x-mu) f(x-mu)]
 *          forward:   D f(x) = U_mu(x) f(x+mu) - f(x)
 *          backward:  D f(x) = f(x) - U_mu^dag(x-mu) f(x-mu)
 * For mu = Tp the link on the last time slice is multiplied by the fermion temporal boundary phase bc (+1 periodic, -1 antiperiodic).
 * Laplacian: sum_{mu in dirs} [U_mu(x) f(x+mu) + U_mu^dag(x-mu) f(x-mu) - 2 f(x)]
 */

enum class DerivType {symmetric, forward, backward};

// gauge link in direction mu, with the fermion temporal boundary phase
template <typename GImpl>
void boundaryLink(typename GImpl::GaugeLinkField &link, const typename GImpl::GaugeField &U, const int mu, const int bc)
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
void covDeriv(Field &out, const Field &in, const typename GImpl::GaugeLinkField &link, const int mu, const DerivType type)
{
    switch (type)
    {
        case DerivType::symmetric:
            out = 0.5*(GImpl::CovShiftForward(link, mu, in) - GImpl::CovShiftBackward(link, mu, in));
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
void laplacian(Field &out, const Field &in, const typename GImpl::GaugeField &U, const std::vector<int> &dirs, const int bc)
{
    typename GImpl::GaugeLinkField link(U.Grid());

    out = Zero();
    for (int mu: dirs)
    {
        boundaryLink<GImpl>(link, U, mu, bc);
        out += GImpl::CovShiftForward(link, mu, in) + GImpl::CovShiftBackward(link, mu, in) - 2.*in;
    }
}

END_HADRONS_NAMESPACE

#endif // Hadrons_MUtilities_DerivHelpers_hpp_
