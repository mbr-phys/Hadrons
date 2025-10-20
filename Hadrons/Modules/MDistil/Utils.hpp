#ifndef Hadrons_MDistil_Utils_hpp_
#define Hadrons_MDistil_Utils_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/A2AMatrix.hpp>
#include <Hadrons/DiskVector.hpp>
#include <Hadrons/TimerArray.hpp>
#include <Hadrons/DistilMatrix.hpp>

using namespace Grid;
using namespace Hadrons;

#define TIME_MOD(t) (((t) + nT) % nT)
#define CLOCK() std::cout << "Contractor : " << std::setw(10) << tAr.getDTimer("total")/1e6 << " s : "

#endif
