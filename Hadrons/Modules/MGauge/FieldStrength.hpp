/*
 * FieldStrength.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
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
#ifndef Hadrons_MGauge_FieldStrength_hpp_
#define Hadrons_MGauge_FieldStrength_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Modules/MUtilities/HQEOperators.hpp>

BEGIN_HADRONS_NAMESPACE

/*
 * Clover field strength as a colour-matrix field
 * ----------------------------------------------
 * mode = clover: out = F_{mu nu}  (traceless, anti-hermitian clover)
 * mode = divE:   out = sum_i [D_i, F_{i t}]  (adjoint covariant divergence of the chromoelectric field)
 * mode = DE:     out = [D_mu, F_{nu t}]       (adjoint covariant derivative of one chromoelectric component)
 *
 * Conventions are those of MUtilities/HQEOperators.hpp. 
 * The output is a GaugeLinkField that can be inserted on a propagator with MUtilities::GaugeInsertion (mode = matrix).
 *
 * Parameters:
 * - gauge: gauge field 
 * - mode:  clover, divE or DE
 * - mu:    first index (clover, DE)
 * - nu:    second index (clover; DE: spatial component of E)
 *
 * There are no default values: gauge, mode and (for clover, DE) mu, nu must be set.
 */

/******************************************************************************
 *                           TFieldStrength                                   *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MGauge)

class FieldStrengthPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(FieldStrengthPar,
                                    std::string,  gauge,
                                    std::string,  mode,
                                    unsigned int, mu,
                                    unsigned int, nu);
};

template <typename GImpl>
class TFieldStrength: public Module<FieldStrengthPar>
{
public:
    GAUGE_TYPE_ALIASES(GImpl,);
public:
    // constructor
    TFieldStrength(const std::string name);
    // destructor
    virtual ~TFieldStrength(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(FieldStrength, TFieldStrength<GIMPL>, MGauge);

/******************************************************************************
 *                     TFieldStrength implementation                          *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename GImpl>
TFieldStrength<GImpl>::TFieldStrength(const std::string name)
: Module<FieldStrengthPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename GImpl>
std::vector<std::string> TFieldStrength<GImpl>::getInput(void)
{
    std::vector<std::string> in = {par().gauge};

    return in;
}

template <typename GImpl>
std::vector<std::string> TFieldStrength<GImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};

    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename GImpl>
void TFieldStrength<GImpl>::setup(void)
{
    if (par().gauge.empty())
    {
        HADRONS_ERROR(Argument, "gauge must be set");
    }
    if ((par().mode != "clover") and (par().mode != "divE") and (par().mode != "DE"))
    {
        HADRONS_ERROR(Argument, "mode must be 'clover', 'divE' or 'DE'");
    }
    if ((par().mode == "DE") and ((par().mu >= Nd) or (par().nu >= Tp)))
    {
        HADRONS_ERROR(Argument, "DE: mu must be in {0, 1, 2, 3} and nu in {0, 1, 2}");
    }

    if (par().mode == "clover")
    {
        if ((par().mu >= Nd) or (par().nu >= Nd))
        {
            HADRONS_ERROR(Argument, "mu and nu must be in {0, 1, 2, 3}");
        }
        if (par().mu == par().nu)
        {
            HADRONS_ERROR(Argument, "mu and nu must differ");
        }
    }

    envCreateLat(GaugeLinkField, getName());
}

// execution ///////////////////////////////////////////////////////////////////
template <typename GImpl>
void TFieldStrength<GImpl>::execute(void)
{
    auto &U   = envGet(GaugeField, par().gauge);
    auto &out = envGet(GaugeLinkField, getName());

    if (par().mode == "clover")
    {
        LOG(Message) << "Computing clover field strength F_{" << par().mu << par().nu << "} of '" << par().gauge << "'" << std::endl;
        HQE::fieldStrength<GImpl>(out, U, par().mu, par().nu);
    }
    else if (par().mode == "DE")
    {
        std::vector<GaugeLinkField> E;

        LOG(Message) << "Computing [D_" << par().mu << ", F_{" << par().nu << " t}] of '" << par().gauge << "'" << std::endl;
        HQE::chromoE<GImpl>(E, U);
        HQE::adjCovDeriv<GImpl>(out, E[par().nu], U, par().mu);
    }
    else
    {
        std::vector<GaugeLinkField> E;

        LOG(Message) << "Computing chromoelectric divergence sum_i [D_i, F_{i t}] of '" << par().gauge << "'" << std::endl;
        HQE::chromoE<GImpl>(E, U);
        HQE::divE<GImpl>(out, E, U);
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MGauge_FieldStrength_hpp_
