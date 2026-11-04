/*
 * GaugeInsertion.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
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
#ifndef Hadrons_MUtilities_GaugeInsertion_hpp_
#define Hadrons_MUtilities_GaugeInsertion_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Modules/MUtilities/HQEOperators.hpp>

BEGIN_HADRONS_NAMESPACE

/*
 * Gluonic insertion on a fermion or propagator field
 * --------------------------------------------------
 * mode = matrix: out = gamma M(x) q(x), M a colour-matrix field (e.g. from MGauge::FieldStrength)
 * mode = sigmaB: out = gamma sum_{i<j<3} sigma_ij F_ij q   (chromomagnetic)
 * mode = sigmaE: out = gamma sum_{j<3}   sigma_j3 F_j3 q   (chromoelectric)
 *
 * F is the traceless clover field strength of 'gauge' and sigma_{mu nu} = (i/2)[gamma_mu, gamma_nu]
 * (see MUtilities/HQEOperators.hpp)
 *
 * Parameters:
 * - q:      input field
 * - gauge:  gauge field (sigmaB, sigmaE)
 * - matrix: colour-matrix field (matrix mode)
 * - mode:   matrix, sigmaB or sigmaE
 * - gamma:  gamma matrix applied after the insertion (Identity for none)
 *
 * There are no default values: q, mode, gamma and the field required by the
 * mode (matrix or gauge) must be set.
 */

/******************************************************************************
 *                          TGaugeInsertion                                   *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MUtilities)

class GaugeInsertionPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(GaugeInsertionPar,
                                    std::string,    q,
                                    std::string,    gauge,
                                    std::string,    matrix,
                                    std::string,    mode,
                                    Gamma::Algebra, gamma);
};

template <typename FImpl, typename GImpl, typename Field>
class TGaugeInsertion: public Module<GaugeInsertionPar>
{
public:
    GAUGE_TYPE_ALIASES(GImpl,);
public:
    // constructor
    TGaugeInsertion(const std::string name);
    // destructor
    virtual ~TGaugeInsertion(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(GaugeInsertionPropagator, ARG(TGaugeInsertion<FIMPL, GIMPL, FIMPL::PropagatorField>), MUtilities);
MODULE_REGISTER_TMP(GaugeInsertionFermion, ARG(TGaugeInsertion<FIMPL, GIMPL, FIMPL::FermionField>), MUtilities);

/******************************************************************************
 *                     TGaugeInsertion implementation                         *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl, typename GImpl, typename Field>
TGaugeInsertion<FImpl, GImpl, Field>::TGaugeInsertion(const std::string name)
: Module<GaugeInsertionPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl, typename GImpl, typename Field>
std::vector<std::string> TGaugeInsertion<FImpl, GImpl, Field>::getInput(void)
{
    std::vector<std::string> in = {par().q};

    if (par().mode == "matrix")
    {
        in.push_back(par().matrix);
    }
    else
    {
        in.push_back(par().gauge);
    }

    return in;
}

template <typename FImpl, typename GImpl, typename Field>
std::vector<std::string> TGaugeInsertion<FImpl, GImpl, Field>::getOutput(void)
{
    std::vector<std::string> out = {getName()};

    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl, typename GImpl, typename Field>
void TGaugeInsertion<FImpl, GImpl, Field>::setup(void)
{
    if ((par().mode != "matrix") and (par().mode != "sigmaB")
        and (par().mode != "sigmaE"))
    {
        HADRONS_ERROR(Argument, "mode must be 'matrix', 'sigmaB' or 'sigmaE'");
    }

    if (par().q.empty())
    {
        HADRONS_ERROR(Argument, "q must be set");
    }

    if ((par().mode == "matrix") and par().matrix.empty())
    {
        HADRONS_ERROR(Argument, "matrix must be set in matrix mode");
    }

    if ((par().mode != "matrix") and par().gauge.empty())
    {
        HADRONS_ERROR(Argument, "gauge must be set in sigmaB/sigmaE mode");
    }

    if (par().gamma == Gamma::Algebra::undef)
    {
        HADRONS_ERROR(Argument, "gamma must be set (Identity for none)");
    }

    envCreateLat(Field, getName());
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl, typename GImpl, typename Field>
void TGaugeInsertion<FImpl, GImpl, Field>::execute(void)
{
    LOG(Message) << "Applying gluonic insertion '" << par().mode << "' with gamma " << par().gamma << " to '" << par().q << "'" << std::endl;

    auto  &q   = envGet(Field, par().q);
    auto  &out = envGet(Field, getName());
    Gamma g(par().gamma);

    if (par().mode == "matrix")
    {
        auto &M = envGet(GaugeLinkField, par().matrix);

        out = M*q;
    }
    else
    {
        auto &U = envGet(GaugeField, par().gauge);
        std::vector<GaugeLinkField> F;

        if (par().mode == "sigmaB")
        {
            HQE::chromoB<GImpl>(F, U);
            HQE::applySigmaB<GImpl>(out, q, F);
        }
        else
        {
            HQE::chromoE<GImpl>(F, U);
            HQE::applySigmaE<GImpl>(out, q, F);
        }
    }
    out = g*out;
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MUtilities_GaugeInsertion_hpp_
