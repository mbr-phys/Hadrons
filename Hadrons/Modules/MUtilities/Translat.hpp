/*
 * Translat.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2022
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
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
#ifndef Hadrons_MUtilities_Translat_hpp_
#define Hadrons_MUtilities_Translat_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                            Field Translat                                  *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MUtilities)

class TranslatPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(TranslatPar,
                                    std::string, field,
                                    std::vector<unsigned int>, xvec);
};

template <typename Field>
class TTranslat: public Module<TranslatPar>
{
public:
    // constructor
    TTranslat(const std::string name);
    // destructor
    virtual ~TTranslat(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(TranslatGaugeField, TTranslat<GIMPL::GaugeField>, MUtilities);
MODULE_REGISTER_TMP(TranslatPropagatorField, TTranslat<FIMPL::PropagatorField>, MUtilities);

/******************************************************************************
 *                     TTranslat implementation                          *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename Field>
TTranslat<Field>::TTranslat(const std::string name)
: Module<TranslatPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename Field>
std::vector<std::string> TTranslat<Field>::getInput(void)
{
    std::vector<std::string> in = {par().field};
    
    return in;
}

template <typename Field>
std::vector<std::string> TTranslat<Field>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename Field>
void TTranslat<Field>::setup(void)
{
    envCreateLat(Field, getName());
    envTmpLat(Field, "tmp");
}

// execution ///////////////////////////////////////////////////////////////////
template <typename Field>
void TTranslat<Field>::execute(void)
{
    LOG(Message) << "Shifting '" << par().field << "' by vector [" 
                 << par().xvec[0] << " " << par().xvec[1] << " "
                 << par().xvec[2] << " " << par().xvec[3] << "]." << std::endl;

    auto &U   = envGet(Field, par().field);
    auto &Utr = envGet(Field, getName());
    envGetTmp(Field, tmp);
    tmp = U;
    
    for (int dir = 0; dir < Nd; dir++) {
        int length = par().xvec[dir];
        if (length > 0) {
            for (int j = 0; j < length; j++) {
                tmp = Cshift(tmp,dir,1);
            }
        }
    }
    Utr = tmp;
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MUtilities_Translat_hpp_
