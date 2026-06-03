/*
* PermCheck.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
*
* Copyright (C) 2015 - 2020
*
* Author: Michael Marshall <Michael.Marshall@ed.ac.uk>
* Author: Antonin Portelli <antonin.portelli@me.com>
* Author: Felix Erben <felix.erben@ed.ac.uk>
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

#ifndef Hadrons_MDistil_PermCheck_hpp_
#define Hadrons_MDistil_PermCheck_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/DilutedNoise.hpp>
#include <Hadrons/DistillationVectors.hpp>

BEGIN_HADRONS_NAMESPACE

BEGIN_MODULE_NAMESPACE(MDistil)

/******************************************************************************
 *                             PermCheck                                    *
 ******************************************************************************/


class PermCheckPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(PermCheckPar,
                                    std::string, inPath1,
                                    std::string, inPath2,
                                    std::string, outPath,
                                    std::string, noisePol,
                                    std::string, timeSources);
};

template <typename FImpl>
class TPermCheck: public Module<PermCheckPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    // constructor
    TPermCheck(const std::string name);
    // destructor
    virtual ~TPermCheck(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
protected:
    unsigned int Ls_;
};

MODULE_REGISTER_TMP(PermCheck, TPermCheck<FIMPL>, MDistil);

/******************************************************************************
 *                 TPermCheck implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TPermCheck<FImpl>::TPermCheck(const std::string name) : Module<PermCheckPar>(name) {}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TPermCheck<FImpl>::getInput(void)
{
    std::vector<std::string> in={ par().noisePol };
    return in;
}


template <typename FImpl>
std::vector<std::string> TPermCheck<FImpl>::getOutput(void)
{
    std::vector<std::string> out{ };
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TPermCheck<FImpl>::setup(void)
{
    GridCartesian * gridHD = envGetGrid(FermionField);
    GridCartesian * gridLD = envGetSliceGrid(FermionField,gridHD->Nd() -1);

    auto &dilNoise = envGet(DistillationNoise<FImpl>, par().noisePol);
    int nDL = dilNoise.dilutionSize(DistillationNoise<FImpl>::Index::l);
    int nDS = dilNoise.dilutionSize(DistillationNoise<FImpl>::Index::s);
    envTmp   (FermionField,    "fermion3dtmp1" ,1, gridLD);
    envTmp   (FermionField,    "fermion3dtmp2" ,1, gridLD);

    Grid::Coordinate coor  = gridHD->GlobalDimensions();
    coor[3] = nDL * nDS;
    Grid::GridCartesian * gridDD = Grid::SpaceTimeGrid::makeFourDimGrid( coor, GridDefaultSimd(Nd,vComplex::Nsimd()), GridDefaultMpi() );
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TPermCheck<FImpl>::execute(void)
{
    // general grid setup
    GridCartesian * gridHD = envGetGrid(FermionField);
    GridCartesian * gridLD = envGetSliceGrid(FermionField,gridHD->Nd() -1);
    const int Ntlocal{gridHD->LocalDimensions()[Tdir]};
    const int Ntfirst{gridHD->LocalStarts()[Tdir]};
    int nT=env().getDim(Tdir);

    envGetTmp(FermionField,    fermion3dtmp1);
    envGetTmp(FermionField,    fermion3dtmp2);

    int Nt{env().getDim(Tdir)};
    auto &dilNoise = envGet(DistillationNoise<FImpl>, par().noisePol);
    int nNoise = dilNoise.size();
    int nVec = dilNoise.getNl();
    int nDL = dilNoise.dilutionSize(DistillationNoise<FImpl>::Index::l);        
    int nDS = dilNoise.dilutionSize(DistillationNoise<FImpl>::Index::s);        
    int nDT = dilNoise.dilutionSize(DistillationNoise<FImpl>::Index::t);        

    std::vector<int> timeSources;
    if (!par().timeSources.empty()) 
    {
        LOG(Message) << "LOADING timeSources = " + par().timeSources << std::endl;
        timeSources = strToVec<int>(par().timeSources); 
        Nt = timeSources.size();
    } 
    else
    {
        std::iota(timeSources.begin(), timeSources.end(), 0);
    }

    LOG(Message) << "PermCheck set up with nDL = " << nDL << ", nDS = " << nDS << ", nDT = " << nDT << std::endl;

    int dk,ds,dSolve,tH,skip;
    std::array<unsigned int, 3> index;

    for (int t = 0; t < Ntlocal; t++ )
    {
        tH = t + Ntfirst;

        ScidacReader pReader;
        std::string filename1 = par().inPath1 + "." + std::to_string(vm().getTrajectory()) + "/t" + std::to_string(t) + "_pkg.bin";
        pReader.open(filename1);

        ScidacReader qReader;
        std::string filename2 = par().inPath2 + "." + std::to_string(vm().getTrajectory()) + "/t" + std::to_string(t) + "_pkg.bin";
        qReader.open(filename2);

        for (int tsrc=0; tsrc<timeSources.size(); tsrc++)
        {
            int tD = timeSources[tsrc];
            for(int id=0; id<nDL * nDS; id++)
            {
                index = dilNoise.dilutionCoordinates(id);
                dk = index[DistillationNoise<FImpl>::Index::l];
                ds = index[DistillationNoise<FImpl>::Index::s];
                dSolve = dilNoise.dilutionIndex(tD,dk,ds);
                skip = 0; 

                DistillationVectorsIo::pkgComponentReader(pReader, fermion3dtmp1, 1, nDL, nDS, nDT, dSolve, skip); 
                DistillationVectorsIo::pkgComponentReader(qReader, fermion3dtmp2, 1, nDL, nDS, nDT, dSolve, skip); 

                RealD nrm = norm2(fermion3dtmp1-fermion3dtmp2);
                LOG(Message) << "INDICES: t = " << t << ", tH = " << tH << ", tD = " << tD 
                             << ", id = " << id << ", dSolve = " << dSolve << ", skip = " << skip
                             << ", NORM = " << nrm << std::endl;

            }
        }
        pReader.close();
        qReader.close();
    }
}

END_MODULE_NAMESPACE
END_HADRONS_NAMESPACE

#endif // Hadrons_MDistil_PermCheck_hpp_
