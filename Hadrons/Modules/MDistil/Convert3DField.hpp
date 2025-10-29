/*
* Convert3DField.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
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

#ifndef Hadrons_MDistil_Convert3DField_hpp_
#define Hadrons_MDistil_Convert3DField_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/DilutedNoise.hpp>
#include <Hadrons/DistillationVectors.hpp>

BEGIN_HADRONS_NAMESPACE

BEGIN_MODULE_NAMESPACE(MDistil)

/******************************************************************************
 *                             Convert3DField                                    *
 ******************************************************************************/


class Convert3DFieldPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(Convert3DFieldPar,
                                    std::string, inPath,
                                    std::string, outPath,
                                    std::string, noisePol,
                                    std::string, timeSources);
};

template <typename FImpl>
class TConvert3DField: public Module<Convert3DFieldPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    // constructor
    TConvert3DField(const std::string name);
    // destructor
    virtual ~TConvert3DField(void) {};
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

MODULE_REGISTER_TMP(Convert3DField, TConvert3DField<FIMPL>, MDistil);

/******************************************************************************
 *                 TConvert3DField implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TConvert3DField<FImpl>::TConvert3DField(const std::string name) : Module<Convert3DFieldPar>(name) {}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TConvert3DField<FImpl>::getInput(void)
{
    std::vector<std::string> in={ par().noisePol };
    return in;
}


template <typename FImpl>
std::vector<std::string> TConvert3DField<FImpl>::getOutput(void)
{
    std::vector<std::string> out{ };
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TConvert3DField<FImpl>::setup(void)
{
    GridCartesian * gridHD = envGetGrid(FermionField);
    GridCartesian * gridLD = envGetSliceGrid(FermionField,gridHD->Nd() -1);

    auto &dilNoise = envGet(DistillationNoise<FImpl>, par().noisePol);
    int nDL = dilNoise.dilutionSize(DistillationNoise<FImpl>::Index::l);
    int nDS = dilNoise.dilutionSize(DistillationNoise<FImpl>::Index::s);
    //envTmp    (std::vector<FermionField>, "vector3d", 1, nDL * nDS, gridLD);
    envTmp   (FermionField,    "fermion3dtmp" ,1, gridLD);

    Grid::Coordinate coor  = gridHD->GlobalDimensions();
    coor[3] = nDL * nDS;
    Grid::GridCartesian * gridDD = Grid::SpaceTimeGrid::makeFourDimGrid( coor, GridDefaultSimd(Nd,vComplex::Nsimd()), GridDefaultMpi() );
    envTmp   (FermionField,    "fermionDDtmp" ,1, gridDD);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TConvert3DField<FImpl>::execute(void)
{
    // general grid setup
    GridCartesian * gridHD = envGetGrid(FermionField);
    GridCartesian * gridLD = envGetSliceGrid(FermionField,gridHD->Nd() -1);
    const int Ntlocal{gridHD->LocalDimensions()[Tdir]};
    const int Ntfirst{gridHD->LocalStarts()[Tdir]};
    int nT=env().getDim(Tdir);

    //envGetTmp(std::vector<FermionField>,    vector3d);
    envGetTmp(FermionField,    fermion3dtmp);
    envGetTmp(FermionField,    fermionDDtmp);

    int Nt{env().getDim(Tdir)};
    auto &dilNoise = envGet(DistillationNoise<FImpl>, par().noisePol);
    int nNoise = dilNoise.size();        
    int nVec = dilNoise.getNl();
    int nDL = dilNoise.dilutionSize(DistillationNoise<FImpl>::Index::l);        
    int nDS = dilNoise.dilutionSize(DistillationNoise<FImpl>::Index::s);        
    int nDT = dilNoise.dilutionSize(DistillationNoise<FImpl>::Index::t);        

    LOG(Message) << "TIMESOURCES = " + par().timeSources << std::endl;
    std::vector<int> timeSources;
    if (!par().timeSources.empty()) 
    {
        LOG(Message) << "LOADING timeSources = " + par().timeSources << std::endl;
        timeSources = strToVec<int>(par().timeSources); 
        Nt = timeSources.size();
        //nDT = timeSources.size();
    } 
    else
    {
        std::iota(timeSources.begin(), timeSources.end(), 0);
    }

    LOG(Message) << "CONVERSION set up with nDL = " << nDL << ", nDS = " << nDS << ", nDT = " << nDT << std::endl;

    int dk,ds,dSolve,tH;
    std::array<unsigned int, 3> index;



    for (int t = 0; t < Ntlocal; t++ )
    {
        tH = t + Ntfirst;
        //for (int tD = 0; tD < Nt; tD++ )
        for (int tD : timeSources)
        {
            startTimer("read I/O");
            for(int id=0; id<nDL * nDS; id++)
            {
                index = dilNoise.dilutionCoordinates(id);
                dk = index[DistillationNoise<FImpl>::Index::l];
                ds = index[DistillationNoise<FImpl>::Index::s];
                dSolve = dilNoise.dilutionIndex(tD,dk,ds);
                /* LOG(Message) << "INDICES: index = [" << index[0] << "," << index[1] << "," << index[2] 
                             << "], dk = " << dk << ", ds = " << ds << ", dSolve = " << dSolve << std::endl; */
                std::string tFileName = par().inPath;
                tFileName.append("_t");
                tFileName.append(std::to_string(tH));
                DistillationVectorsIo::readComponent(fermion3dtmp, tFileName, 1, nDL, nDS, nDT, dSolve, vm().getTrajectory());
                // this is vector 2 on timeslice tH 
                //ExtractSliceLocal(fermion3dtmp2,fermion4dtmp,0,t,Tdir);
                //vector3d[id]=fermion3dtmp;
                InsertSliceLocal(fermion3dtmp,fermionDDtmp,0,id,Tdir);
            }
            stopTimer("read I/O");

            /*startTimer("write I/O");
            std::string tFileName = par().outPath;
            tFileName.append("_tSm");
            tFileName.append(std::to_string(tD));
            tFileName.append("_tLoc");
            tFileName.append(std::to_string(tH));
            std::vector<int> tS;
            DistillationVectorsIo::write(tFileName, vector3d, "3d vector", 1, nDL, nDS, nDT, tS, false, vm().getTrajectory());
            stopTimer("write I/O");*/
            startTimer("write I/O");
            std::string tFileName = par().outPath;
            tFileName.append("_DD");
            tFileName.append("_tSm");
            tFileName.append(std::to_string(tD));
            tFileName.append("_tLoc");
            tFileName.append(std::to_string(tH));
            std::vector<int> tS;
            DistillationVectorsIo::writeComponent(tFileName, fermionDDtmp, "unsmSolve", 1, nDL, nDS, nDT, tS, 0, vm().getTrajectory());
            stopTimer("write I/O");

        }
    }

}

END_MODULE_NAMESPACE
END_HADRONS_NAMESPACE

#endif // Hadrons_MDistil_Convert3DField_hpp_
