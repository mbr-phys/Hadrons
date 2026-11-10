/*
 * TimeBlockSpinColorDiagonal.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
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
#ifndef Hadrons_MNoise_TimeBlockSpinColorDiagonal_hpp_
#define Hadrons_MNoise_TimeBlockSpinColorDiagonal_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/DilutedNoise.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                       TimeBlockSpinColorDiagonal                           *
 ******************************************************************************
 * Spin-colour diagonal Z4 noise diluted in blocks of time slices.
 * - nsrc: number of base noises
 * - blockStarts: first time slice of each block, starting at 0 and strictly
 *   increasing; block b covers [blockStarts[b], blockStarts[b + 1]), the last
 *   block ends at Nt.
 * Component i is base noise i / nBlock on block i % nBlock (see
 * TimeBlockNoise), giving nsrc*nBlock components.
 */
BEGIN_MODULE_NAMESPACE(MNoise)

class TimeBlockSpinColorDiagonalPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(TimeBlockSpinColorDiagonalPar,
                                    unsigned int,              nsrc,
                                    std::vector<unsigned int>, blockStarts);
};

template <typename FImpl>
class TTimeBlockSpinColorDiagonal: public Module<TimeBlockSpinColorDiagonalPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TTimeBlockSpinColorDiagonal(const std::string name);
    // destructor
    virtual ~TTimeBlockSpinColorDiagonal(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(TimeBlockSpinColorDiagonal, TTimeBlockSpinColorDiagonal<FIMPL>, MNoise);
MODULE_REGISTER_TMP(ZTimeBlockSpinColorDiagonal, TTimeBlockSpinColorDiagonal<ZFIMPL>, MNoise);

/******************************************************************************
 *                 TTimeBlockSpinColorDiagonal implementation                 *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TTimeBlockSpinColorDiagonal<FImpl>::TTimeBlockSpinColorDiagonal(const std::string name)
: Module<TimeBlockSpinColorDiagonalPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TTimeBlockSpinColorDiagonal<FImpl>::getInput(void)
{
    std::vector<std::string> in;
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TTimeBlockSpinColorDiagonal<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TTimeBlockSpinColorDiagonal<FImpl>::setup(void)
{
    envCreateDerived(SpinColorDiagonalNoise<FImpl>, 
                     TimeBlockNoise<FImpl>,
                     getName(), 1, envGetGrid(FermionField), par().nsrc, par().blockStarts);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TTimeBlockSpinColorDiagonal<FImpl>::execute(void)
{
    auto &noise = envGetDerived(SpinColorDiagonalNoise<FImpl>, TimeBlockNoise<FImpl>, getName());

    LOG(Message) << "Generating time-block diluted spin-color diagonal noise with"
                 << " nsrc = " << par().nsrc << " and " << noise.nBlock()
                 << " blocks:";
    for (unsigned int b = 0; b < noise.nBlock(); ++b)
    {
        LOG(Message) << " [" << noise.blockBegin(b) << ", " << noise.blockEnd(b) << ")";
    }
    LOG(Message) << std::endl;
    noise.generateNoise(rng4d());
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MNoise_TimeBlockSpinColorDiagonal_hpp_
