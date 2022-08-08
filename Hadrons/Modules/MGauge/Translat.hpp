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
#ifndef Hadrons_MGauge_Translat_hpp_
#define Hadrons_MGauge_Translat_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                            Gauge translat                                  *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MGauge)

class TranslatPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(TranslatPar,
                                    std::string, gauge,
                                    std::vector<unsigned int>, xvec);
};

template <typename GImpl>
class TTranslat: public Module<TranslatPar>
{
public:
    INHERIT_GIMPL_TYPES(GImpl);
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

MODULE_REGISTER_TMP(Translat, TTranslat<GIMPL>, MGauge);

/******************************************************************************
 *                     TTranslat implementation                          *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename GImpl>
TTranslat<GImpl>::TTranslat(const std::string name)
: Module<TranslatPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename GImpl>
std::vector<std::string> TTranslat<GImpl>::getInput(void)
{
    std::vector<std::string> in = {par().gauge};
    
    return in;
}

template <typename GImpl>
std::vector<std::string> TTranslat<GImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename GImpl>
void TTranslat<GImpl>::setup(void)
{
    envCreateLat(GaugeField, getName());
}

// execution ///////////////////////////////////////////////////////////////////
template <typename GImpl>
void TTranslat<GImpl>::execute(void)
{
    LOG(Message) << "Shifting '" << par().gauge << "' by vector [" 
                 << par().xvec[0] << " " << par().xvec[1] << " "
                 << par().xvec[2] << " " << par().xvec[3] << "]." << std::endl;

    auto &U   = envGet(GaugeField, par().gauge);
    auto &Utr = envGet(GaugeField, getName());
    
    std::vector<GaugeLinkField> Umu(Nd,U.Grid()), tmp(Nd,U.Grid());
    for (int mu = 0; mu < Nd; mu++) {
        Umu[mu] = PeekIndex<LorentzIndex>(U,mu);
    }
    for (int dir = 0; dir < Nd; dir++) {
        int length = par().xvec[dir];
        if (length > 0) {
            for (int nu = 0; nu < Nd; nu++) {
                for (int j = 0; j < length; j++) {
                    tmp[nu] = Cshift(Umu[nu],dir,1);
                    Umu[nu] = tmp[nu];
                }
            }
        }
    }
    for (int mu = 0; mu < Nd; mu++) {
        PokeIndex<LorentzIndex>(Utr, Umu[mu], mu);
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MGauge_Translat_hpp_
