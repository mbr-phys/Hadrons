/*
* MomentumProjectorList.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
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

#ifndef Hadrons_MUtilities_MomentumProjectorList_hpp_
#define Hadrons_MUtilities_MomentumProjectorList_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Serialization.hpp>
#include <Grid/algorithms/blas/MomentumProject.h>

BEGIN_HADRONS_NAMESPACE

/*
 * Momentum projector from an explicit list of momenta
 * ---------------------------------------------------
 * Same products as MomentumProjectorSphere, so any consumer with a
 * momProjector parameter (Gamma3pt, Meson, StochasticMeson,
 * ContractConservedSink) can use either:
 *   - <name>:         the MomentumProject object,
 *   - <name>_momList: the momenta, as integer vectors.
 *
 * Parameters:
 * - momenta: spatial momenta in units of 2 pi / L, one "px py pz" string each.
 * - sign:    +1 or -1 (required); the phase is exp(sign * 2 pi i p.x / L).
 *            The two choices are complex conjugates of each other: sign picks
 *            whether x is treated as a source-type or a sink-type position.
 *            -1 matches MomentumProjectorSphere (and so Gamma3pt's labels),
 *            +1 matches MSink::Point.
 *
 * The phase fields are temporaries; only the BLAS copy inside the
 * MomentumProject object stays resident.
 */

/******************************************************************************
 *                           MomentumProjectorList                            *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MUtilities)

class MomentumProjectorListPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(MomentumProjectorListPar,
                                    std::vector<std::string>, momenta,
                                    int,                      sign);
};

template <typename ComplexField>
class TMomentumProjectorList: public Module<MomentumProjectorListPar>
{
public:
    using MPType = MomentumProject<ComplexField, ComplexField>;
public:
    // constructor
    TMomentumProjectorList(const std::string name);
    // destructor
    virtual ~TMomentumProjectorList(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
protected:
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(MomentumProjectorList, TMomentumProjectorList<FIMPL::ComplexField>, MUtilities);

/******************************************************************************
 *                    TMomentumProjectorList implementation                   *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename ComplexField>
TMomentumProjectorList<ComplexField>::TMomentumProjectorList(const std::string name)
: Module<MomentumProjectorListPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename ComplexField>
std::vector<std::string> TMomentumProjectorList<ComplexField>::getInput(void)
{
    std::vector<std::string> in;

    return in;
}

template <typename ComplexField>
std::vector<std::string> TMomentumProjectorList<ComplexField>::getOutput(void)
{
    std::vector<std::string> out = {getName(), getName() + "_momList"};

    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename ComplexField>
void TMomentumProjectorList<ComplexField>::setup(void)
{
    const unsigned int nd = env().getNd() - 1;

    if (par().momenta.empty())
    {
        HADRONS_ERROR(Argument, "momenta cannot be empty");
    }
    if ((par().sign != 1) && (par().sign != -1))
    {
        HADRONS_ERROR(Argument, "sign must be +1 or -1");
    }
    env().template createObject<MPType>(getName(), Environment::Storage::standard, 0);
    envCreate(std::vector<std::vector<int>>, getName() + "_momList", 1, 0);

    auto &momList = envGet(std::vector<std::vector<int>>, getName() + "_momList");

    for (auto &s: par().momenta)
    {
        std::vector<int> p = strToVec<int>(s);

        if (p.size() != nd)
        {
            HADRONS_ERROR(Size, "momentum '" + s + "' does not have "
                          + std::to_string(nd) + " components");
        }
        momList.push_back(p);
    }
    envTmpLat(ComplexField, "coor");
    envTmp(std::vector<ComplexField>, "phases", 1, momList.size(),
           envGetGrid(ComplexField));
}

// execution ///////////////////////////////////////////////////////////////////
template <typename ComplexField>
void TMomentumProjectorList<ComplexField>::execute(void)
{
    const Complex i(0.0, 1.0);
    auto          &MP      = envGet(MPType, getName());
    auto          &momList = envGet(std::vector<std::vector<int>>, getName() + "_momList");

    envGetTmp(ComplexField, coor);
    envGetTmp(std::vector<ComplexField>, phases);
    LOG(Message) << "Building momentum projector for " << momList.size()
                 << " momenta, phase exp(" << ((par().sign > 0) ? "+" : "-")
                 << "2 pi i p.x/L)" << std::endl;
    startTimer("Phases");
    for (unsigned int m = 0; m < momList.size(); ++m)
    {
        ComplexField &ph = phases[m];
        const auto   &p  = momList[m];

        ph = Zero();
        for (unsigned int mu = 0; mu < p.size(); ++mu)
        {
            LatticeCoordinate(coor, mu);
            ph = ph + (static_cast<Real>(p[mu]) / static_cast<Real>(env().getDim(mu))) * coor;
        }
        ph = exp(static_cast<Real>(2*M_PI*par().sign)*i*ph);
    }
    stopTimer("Phases");

    startTimer("ImportPhases");
    MP.Allocate(momList.size(), env().getGrid());
    MP.ImportMomenta(phases);
    stopTimer("ImportPhases");
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MUtilities_MomentumProjectorList_hpp_
