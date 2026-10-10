/*
 * LoadOpenQcd.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2023
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Fabian Joswig <fabian.joswig@ed.ac.uk>
 * Author: Fabian Joswig <fabian.joswig@wwu.de>
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
#ifndef Hadrons_MIO_LoadOpenQcd_hpp_
#define Hadrons_MIO_LoadOpenQcd_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/GaugeBoundary.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Serialization.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 Load an OpenQcd configuration

 file    Namestem of the files to read in. In contrast to the LoadNersc
         module, the separator between the namestem and the configuration
         number is 'n' instead of '.' compliant with the openQCD convention.         
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MIO)

class LoadOpenQcdPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(LoadOpenQcdPar,
                                    std::string, file,
                                    std::string, gaugeBC,
                                    double, cG,
                                    int, timeDir);
};

template <typename GImpl>
class TLoadOpenQcd: public Module<LoadOpenQcdPar>
{
public:
    GAUGE_TYPE_ALIASES(GImpl,);
public:
    // constructor
    TLoadOpenQcd(const std::string name);
    // destructor
    virtual ~TLoadOpenQcd(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(LoadOpenQcd,  TLoadOpenQcd<GIMPL>,  MIO);

/******************************************************************************
*                       TLoadOpenQcd implementation                             *
******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename GImpl>
TLoadOpenQcd<GImpl>::TLoadOpenQcd(const std::string name)
: Module<LoadOpenQcdPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename GImpl>
std::vector<std::string> TLoadOpenQcd<GImpl>::getInput(void)
{
    std::vector<std::string> in;
    
    return in;
}

template <typename GImpl>
std::vector<std::string> TLoadOpenQcd<GImpl>::getOutput(void)
{
    return {getName(), GaugeBoundary::metadataName(getName())};
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename GImpl>
void TLoadOpenQcd<GImpl>::setup(void)
{
    if (par().gaugeBC.empty()) {
        HADRONS_ERROR(Argument,
                      "gaugeBC must be explicitly set to 'unknown', "
                      "'periodic', or 'open'");
    }
    if (!GaugeBoundary::isValidType(par().gaugeBC)) {
        HADRONS_ERROR(Argument,
                      "gaugeBC must be 'unknown', 'periodic', or 'open'");
    }
    if (par().timeDir < 0 || par().timeDir >= Nd) {
        HADRONS_ERROR(Argument, "timeDir is outside the lattice dimensions");
    }
    if (!std::isfinite(par().cG) || par().cG < 0.0) {
        HADRONS_ERROR(Argument,
                      "cG must be explicitly set to a finite non-negative value");
    }

    const bool userSpecified = par().gaugeBC != GaugeBoundary::unknown;
    const double effectiveCG = (par().gaugeBC == GaugeBoundary::open)
                             ? par().cG : 1.0;
    envCreateLat(GaugeField, getName());
    envCreate(GaugeBoundaryInfo, GaugeBoundary::metadataName(getName()), 1,
              par().gaugeBC, par().timeDir, effectiveCG, userSpecified);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename GImpl>
void TLoadOpenQcd<GImpl>::execute(void)
{
    FieldMetaData header;
    std::string   fileName = par().file + "n"
                             + std::to_string(vm().getTrajectory());
    LOG(Message) << "Loading OpenQcd configuration from file '" << fileName
                 << "'" << std::endl;

    LOG(Message) << "Declared temporal gauge boundary: " << par().gaugeBC;
    if (par().gaugeBC == GaugeBoundary::open) {
        LOG(Message) << " (timeDir=" << par().timeDir
                     << ", cG=" << par().cG << ")";
    }
    LOG(Message) << std::endl;

    auto &U = envGet(GaugeField, getName());
    OpenQcdIO::readConfiguration(U, header, fileName);

    const RealD naivePlaquette = WilsonLoops<GImpl>::avgPlaquette(U);
    LOG(Message) << "Loaded configuration naive plaquette (Grid-normalised): "
                 << naivePlaquette << std::endl;
    LOG(Message) << "OpenQCD header plaquette (Grid-normalised): "
                 << header.plaquette << std::endl;

    if (par().gaugeBC == GaugeBoundary::unknown) {
        LOG(Warning) << "The openQCD configuration header does not encode the "
                     << "gauge boundary type; downstream modules must use an "
                     << "explicit boundary declaration." << std::endl;
    }

    if (par().gaugeBC == GaugeBoundary::open) {
        using Boundary = OpenGaugeBoundary<GImpl>;
        const Coordinate dimensions = U.Grid()->GlobalDimensions();
        const int nt = dimensions[par().timeDir];
        if (nt < 2) {
            HADRONS_ERROR(Size,
                          "open temporal boundary requires extent at least 2");
        }

        const RealD openQcdWeight =
            GaugeBoundary::openQcdPlaquetteWeight(dimensions, par().timeDir);
        const RealD actionWeight =
            GaugeBoundary::actionPlaquetteWeight(dimensions, par().timeDir,
                                                  par().cG);

        const auto loops = Boundary::plaquetteLoops(1.0);
        Boundary openQcdBoundary(
            OpenGaugeBoundaryParameters::Open(par().timeDir, 1.0));
        const RealD openQcdAction = openQcdBoundary.action(U, loops);
        const RealD weightedGrid = 1.0 - openQcdAction / openQcdWeight;
        const RealD weightedRaw = RealD(GImpl::num_colours) * weightedGrid;

        Boundary actionBoundary(
            OpenGaugeBoundaryParameters::Open(par().timeDir, par().cG));
        const RealD cGAction = actionBoundary.action(U, loops);
        const RealD actionWeighted = 1.0 - cGAction / actionWeight;

        LOG(Message) << "openQCD OBC-weighted plaquette (Grid-normalised): "
                     << weightedGrid << std::endl;
        LOG(Message) << "openQCD OBC-weighted plaquette (raw colour trace): "
                     << weightedRaw << std::endl;
        LOG(Message) << "cG action-weighted plaquette (Grid-normalised, cG="
                     << par().cG << "): " << actionWeighted << std::endl;
        if (std::abs(par().cG - 1.0) <= 1.0e-14) {
            LOG(Message) << "For cG=1 the openQCD-weighted and action-weighted "
                         << "plaquettes agree by construction." << std::endl;
        } else {
            LOG(Message) << "The cG action-weighted value differs from the "
                         << "standard openQCD plaq_wsum diagnostic because "
                         << "plaq_wsum always uses boundary weight 1/2."
                         << std::endl;
        }
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MIO_LoadOpenQcd_hpp_
