/*
 * DslashField.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
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
#ifndef Hadrons_MUtilities_DslashField_hpp_
#define Hadrons_MUtilities_DslashField_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Serialization.hpp>

BEGIN_HADRONS_NAMESPACE

/*
 * Covariant derivative operator Dslash = gamma^mu D_mu
 * ----------------------------------------------------
 * Applies the symmetric covariant derivative to a fermion field:
 *   Dslash psi(x) = 1/2 sum_mu gamma_mu [U_mu(x) psi(x+mu) - U^dagger_mu(x-mu) psi(x-mu)]
 * 
 * Parameters:
 * - input: input field (FermionField or PropagatorField)
 * - gauge: gauge field at same flow time as input
 * - output: name for output field (default: <name>_out)
 */

/******************************************************************************
 *                            TDslashField                                    *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MUtilities)

class DslashFieldPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(DslashFieldPar,
                                    std::string, input,
                                    std::string, gauge,
                                    std::string, output);
};

template <typename FImpl, typename Field>
class TDslashField: public Module<DslashFieldPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TDslashField(const std::string name);
    // destructor
    virtual ~TDslashField(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
    
private:
    // Helper function to compute Dslash on a field
    void computeDslash(Field &out, const Field &in, const GaugeField &U);
};

MODULE_REGISTER_TMP(DslashFieldFermion, 
                    ARG(TDslashField<FIMPL, FIMPL::FermionField>), 
                    MUtilities);
MODULE_REGISTER_TMP(DslashFieldPropagator, 
                    ARG(TDslashField<FIMPL, FIMPL::PropagatorField>), 
                    MUtilities);

/******************************************************************************
 *                       TDslashField implementation                          *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl, typename Field>
TDslashField<FImpl, Field>::TDslashField(const std::string name)
: Module<DslashFieldPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl, typename Field>
std::vector<std::string> TDslashField<FImpl, Field>::getInput(void)
{
    std::vector<std::string> in = {par().input, par().gauge};
    return in;
}

template <typename FImpl, typename Field>
std::vector<std::string> TDslashField<FImpl, Field>::getOutput(void)
{
    std::string outName = par().output;
    if (outName.empty()) {
        outName = getName() + "_out";
    }
    std::vector<std::string> out = {outName};
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl, typename Field>
void TDslashField<FImpl, Field>::setup(void)
{
    std::string outName = par().output;
    if (outName.empty()) {
        outName = getName() + "_out";
    }
    envCreateLat(Field, outName);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl, typename Field>
void TDslashField<FImpl, Field>::execute(void)
{
    LOG(Message) << "Computing Dslash on '" << par().input 
                 << "' using gauge field '" << par().gauge << "'." << std::endl;

    auto &in = envGet(Field, par().input);
    auto &U = envGet(GaugeField, par().gauge);
    
    std::string outName = par().output;
    if (outName.empty()) {
        outName = getName() + "_out";
    }
    auto &out = envGet(Field, outName);
    
    // Compute Dslash = gamma^mu D_mu (symmetric covariant derivative)
    computeDslash(out, in, U);
    
    LOG(Message) << "Dslash computation complete. Output: '" << outName << "'" << std::endl;
}

// helper: compute Dslash //////////////////////////////////////////////////////
template <typename FImpl, typename Field>
void TDslashField<FImpl, Field>::computeDslash(Field &out, const Field &in, const GaugeField &U)
{
    // Dslash psi(x) = 1/2 sum_mu gamma_mu [U_mu(x) psi(x+mu) - U^dagger_mu(x-mu) psi(x-mu)]
    // This follows the NPRUtils::dslash convention
    
    out = Zero();
    Field tmp(in.Grid());
    
    for (int mu = 0; mu < Nd; mu++) {
        auto U_mu = peekLorentz(U, mu);
        
        tmp = FImpl::CovShiftForward(U_mu, mu, in);
        tmp = tmp - FImpl::CovShiftBackward(U_mu, mu, in);
        out += Gamma::gmu[mu] * tmp;
    }
    
    out = 0.5 * out;
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MUtilities_DslashField_hpp_
