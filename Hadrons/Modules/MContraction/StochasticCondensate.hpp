/*
 * StochasticCondensate.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
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
#ifndef Hadrons_MContraction_StochasticCondensate_hpp_
#define Hadrons_MContraction_StochasticCondensate_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Serialization.hpp>

BEGIN_HADRONS_NAMESPACE

/*
 * Stochastic condensate from noise-source contractions
 * ----------------------------------------------------
 * Computes: -⟨eta^dagger Gamma phi⟩ + c_fl ⟨eta^dagger eta⟩ (c_fl only for scalar channel at positive GF time)
 * 
 * Parameters:
 * - eta: noise field (FermionField or PropagatorField)
 * - phi: solution field (same type as eta)
 * - gamma: gamma matrix insertion (default: "Identity")
 * - c_fl: flow-time O(a) improvement coefficient (default: 0; only used for gamma="Identity")
 * 
 * - Scalar condensate: gamma="Identity", c_fl=0.5 (tree-level Wilson) or 0 (DWF)
 * - Pseudoscalar: gamma="Gamma5", c_fl=0
 * - Derivative condensate: use DslashField first, then gamma="Identity"
 */

/******************************************************************************
 *                         TStochasticCondensate                              *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MContraction)

class StochasticCondensatePar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(StochasticCondensatePar,
                                    std::string, eta,
                                    std::string, phi,
                                    Gamma::Algebra, gamma,
                                    double, c_fl,
                                    std::string, output);
};

template <typename FImpl, typename Field>
class TStochasticCondensate: public Module<StochasticCondensatePar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    class Result: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(Result,
                                        Gamma::Algebra, gamma,
                                        double, c_fl,
                                        Complex, condensate);
    };
public:
    // constructor
    TStochasticCondensate(const std::string name);
    // destructor
    virtual ~TStochasticCondensate(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    virtual std::vector<std::string> getOutputFiles(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(StochasticCondensateFermion, 
                    ARG(TStochasticCondensate<FIMPL, FIMPL::FermionField>), 
                    MContraction);
MODULE_REGISTER_TMP(StochasticCondensatePropagator, 
                    ARG(TStochasticCondensate<FIMPL, FIMPL::PropagatorField>), 
                    MContraction);

/******************************************************************************
 *                     TStochasticCondensate implementation                   *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl, typename Field>
TStochasticCondensate<FImpl, Field>::TStochasticCondensate(const std::string name)
: Module<StochasticCondensatePar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl, typename Field>
std::vector<std::string> TStochasticCondensate<FImpl, Field>::getInput(void)
{
    std::vector<std::string> in = {par().eta, par().phi};
    return in;
}

template <typename FImpl, typename Field>
std::vector<std::string> TStochasticCondensate<FImpl, Field>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    return out;
}

template <typename FImpl, typename Field>
std::vector<std::string> TStochasticCondensate<FImpl, Field>::getOutputFiles(void)
{
    std::vector<std::string> output;
    if (!par().output.empty())
        output.push_back(resultFilename(par().output));
    return output;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl, typename Field>
void TStochasticCondensate<FImpl, Field>::setup(void)
{
    envCreate(HadronsSerializable, getName(), 1, 0);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl, typename Field>
void TStochasticCondensate<FImpl, Field>::execute(void)
{
    LOG(Message) << "Computing stochastic condensate '" << getName() 
                 << "' using eta='" << par().eta << "' and phi='" << par().phi 
                 << "' with gamma=" << par().gamma
                 << " and c_fl=" << par().c_fl << "." << std::endl;

    auto &eta = envGet(Field, par().eta);
    auto &phi = envGet(Field, par().phi);
    
    Gamma G(par().gamma);
    
    LatticeComplex integrand = -localInnerProduct(eta, closure(G * phi));
    
    // Add c_fl term only for scalar channel (gamma = Identity)
    if (par().c_fl != 0.0 && par().gamma == Gamma::Algebra::Identity) {
        LatticeComplex eta_norm2 = localInnerProduct(eta, eta);
        integrand += par().c_fl * eta_norm2;
    } else if (par().c_fl != 0.0 && par().gamma != Gamma::Algebra::Identity) {
        LOG(Warning) << "c_fl term ignored for non-scalar gamma structure" << std::endl;
    }
    
    Complex condensate = TensorRemove(sum(integrand));
    
    Result result;
    result.gamma = par().gamma;
    result.c_fl = par().c_fl;
    result.condensate = condensate;
    
    auto &out = envGet(HadronsSerializable, getName());
    out = result;
    
    saveResult(par().output, "condensate", result);
    
    LOG(Message) << "Condensate = " << condensate << std::endl;
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_StochasticCondensate_hpp_
