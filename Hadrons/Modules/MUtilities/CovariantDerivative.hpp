/*
 * CovariantDerivative.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
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
#ifndef Hadrons_MUtilities_CovariantDerivative_hpp_
#define Hadrons_MUtilities_CovariantDerivative_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Modules/MUtilities/DerivHelpers.hpp>

BEGIN_HADRONS_NAMESPACE

/*
 * Covariant derivative of a fermion or propagator field
 * -----------------------------------------------------
 * type = symmetric | forward | backward:
 *     out = D_mu q, with mu the single entry of dirs
 * type = laplacian:
 *     out = sum_{mu in dirs} [U_mu(x) q(x+mu) + U_mu^dag(x-mu) q(x-mu) - 2 q(x)]
 *     (dirs = "0 1 2" gives the spatial Laplacian)
 *
 * Conventions are those of MUtilities/DerivHelpers.hpp. 
 * bc is the fermion temporal boundary phase applied to the last time slice links 
 * when the time direction is involved.
 *
 * Parameters:
 * - q:     input field
 * - gauge: gauge field (use the flowed one at the flow time of q)
 * - type:  symmetric, forward, backward or laplacian
 * - dirs:  direction(s), e.g. "3" or "0 1 2"
 * - bc:    temporal boundary phase (+1 or -1)
 *
 * There are no default values: all parameters must be set.
 */

/******************************************************************************
 *                       TCovariantDerivative                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MUtilities)

class CovariantDerivativePar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(CovariantDerivativePar,
                                    std::string, q,
                                    std::string, gauge,
                                    std::string, type,
                                    std::string, dirs,
                                    int,         bc);
};

inline DerivType covariantDerivativeType(const std::string &s)
{
    if (s == "symmetric")
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
                  + "' (expected symmetric, forward, backward or laplacian)");

    return DerivType::symmetric;
}

template <typename FImpl, typename GImpl, typename Field>
class TCovariantDerivative: public Module<CovariantDerivativePar>
{
public:
    GAUGE_TYPE_ALIASES(GImpl,);
public:
    // constructor
    TCovariantDerivative(const std::string name);
    // destructor
    virtual ~TCovariantDerivative(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(CovariantDerivativePropagator,
                    ARG(TCovariantDerivative<FIMPL, GIMPL, FIMPL::PropagatorField>),
                    MUtilities);
MODULE_REGISTER_TMP(CovariantDerivativeFermion,
                    ARG(TCovariantDerivative<FIMPL, GIMPL, FIMPL::FermionField>),
                    MUtilities);

/******************************************************************************
 *                  TCovariantDerivative implementation                       *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl, typename GImpl, typename Field>
TCovariantDerivative<FImpl, GImpl, Field>::TCovariantDerivative(const std::string name)
: Module<CovariantDerivativePar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl, typename GImpl, typename Field>
std::vector<std::string> TCovariantDerivative<FImpl, GImpl, Field>::getInput(void)
{
    std::vector<std::string> in = {par().q, par().gauge};

    return in;
}

template <typename FImpl, typename GImpl, typename Field>
std::vector<std::string> TCovariantDerivative<FImpl, GImpl, Field>::getOutput(void)
{
    std::vector<std::string> out = {getName()};

    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl, typename GImpl, typename Field>
void TCovariantDerivative<FImpl, GImpl, Field>::setup(void)
{
    if (par().q.empty() or par().gauge.empty())
    {
        HADRONS_ERROR(Argument, "q and gauge must be set");
    }
    if ((par().bc != 1) and (par().bc != -1))
    {
        HADRONS_ERROR(Argument, "bc must be +1 or -1");
    }
    std::vector<int> dirs = strToVec<int>(par().dirs);

    if (dirs.empty())
    {
        HADRONS_ERROR(Argument, "dirs must be set");
    }
    for (int mu: dirs)
    {
        if ((mu < 0) or (mu >= Nd))
        {
            HADRONS_ERROR(Argument, "direction out of range: " + par().dirs);
        }
    }
    if (par().type != "laplacian")
    {
        covariantDerivativeType(par().type);
        if (dirs.size() != 1)
        {
            HADRONS_ERROR(Argument, "type '" + par().type + "' needs exactly one direction in dirs");
        }
    }
    envCreateLat(Field, getName());
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl, typename GImpl, typename Field>
void TCovariantDerivative<FImpl, GImpl, Field>::execute(void)
{
    std::vector<int> dirs = strToVec<int>(par().dirs);

    LOG(Message) << "Applying covariant derivative (" << par().type
                 << ", dirs '" << par().dirs << "', bc " << par().bc
                 << ") to '" << par().q << "' with gauge field '"
                 << par().gauge << "'" << std::endl;

    auto &q  = envGet(Field, par().q);
    auto &U   = envGet(GaugeField, par().gauge);
    auto &out = envGet(Field, getName());

    if (par().type == "laplacian")
    {
        laplacian<GImpl>(out, q, U, dirs, par().bc);
    }
    else
    {
        GaugeLinkField link(U.Grid());

        boundaryLink<GImpl>(link, U, dirs[0], par().bc);
        covDeriv<GImpl>(out, q, link, dirs[0], covariantDerivativeType(par().type));
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MUtilities_CovariantDerivative_hpp_
