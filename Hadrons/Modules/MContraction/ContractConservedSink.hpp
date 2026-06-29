/*
 * ContractConservedSink.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2023
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
#ifndef Hadrons_MContraction_ContractConservedSink_hpp_
#define Hadrons_MContraction_ContractConservedSink_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/*
  Conserved Current insertions
 -----------------------------
 
 * options:
 - prop:       propagator. Must match the action, i.e. 5D action needs 5D propagator
 - action:     action module used for propagator solution (string)
 - source:     source module for the quark, used to remove contact terms (string)
 - mom:        momentum phase to put at sink (string)
 - current:    current to be inserted (Current::Vector or Current::Axial)
 - dir:        index for direction to take current in (unsigned int, 0 1 2 3)
*/

/******************************************************************************
 *                              ContractConservedSink                                  *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MContraction)

class ContractConservedSinkPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(ContractConservedSinkPar,
                                    std::string,    prop,   // Name of the propagator we are checking Ward identity
                                    std::string,    action,
                                    std::string,    source,
                                    std::string,    mom,
                                    Current,        current,
                                    Gamma::Algebra, gSrc,
                                    unsigned int,   dir);
};

template <typename FImpl>
class TContractConservedSink: public Module<ContractConservedSinkPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    class Result: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(Result,
                                        Current,              current, 
                                        Gamma::Algebra,       gamma_src,
                                        std::vector<Complex>, corr);
    };
public:
    // constructor
    TContractConservedSink(const std::string name);
    // destructor
    virtual ~TContractConservedSink(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
protected:
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
    // Perform Slice Sum and then save delta
    void SliceOut(std::vector<Complex> &Out, SlicedComplex &Sum, const ComplexField &f, bool bDiff) const
    {
        sliceSum(f, Sum, Tp);
        const auto nt = Sum.size();
        for (size_t t = 0; t < nt; ++t)
        {
            Out[t] = TensorRemove(bDiff ? Sum[t] - Sum[(t-1+nt)%nt] : Sum[t]);
        }
    }
private:
    unsigned int Ls_;
    std::string momphName_;
};

MODULE_REGISTER_TMP(ContractConservedSink, TContractConservedSink<FIMPL>, MContraction);
MODULE_REGISTER_TMP(ZContractConservedSink, TContractConservedSink<ZFIMPL>, MContraction);

/******************************************************************************
 *                     TContractConservedSink implementation                           *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TContractConservedSink<FImpl>::TContractConservedSink(const std::string name)
: Module<ContractConservedSinkPar>(name)
, momphName_ (name + "_momph")
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TContractConservedSink<FImpl>::getInput(void)
{
    return { par().prop, par().action, par().source };
}

template <typename FImpl>
std::vector<std::string> TContractConservedSink<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TContractConservedSink<FImpl>::setup(void)
{
    // The propagator can be 4d or 5d, but must match the action
    const unsigned int ActionLs_{ env().getObjectLs(par().action) };
    Ls_ = env().getObjectLs( par().prop );
    if (Ls_ != ActionLs_)
    {
        std::string sError{ "Ls mismatch: propagator Ls="};
        sError.append( std::to_string( Ls_ ) );
        sError.append( ", action Ls=" );
        sError.append( std::to_string( ActionLs_ ) );
        HADRONS_ERROR(Size, sError);
    }
    // These temporaries are always 4d
    envTmpLat(PropagatorField, "tmp");
    envTmpLat(ComplexField, "tmp_current");
    envCreateLat(HadronsSerializable, getName(), 1, 0);

    envTmpLat(LatticeComplex, "coor");
    envCacheLat(LatticeComplex, momphName_);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TContractConservedSink<FImpl>::execute(void)
{
    LOG(Message) << "Performing conserved current insertion for propagator " << par().prop << std::endl;
    auto &prop = envGet(PropagatorField, par().prop);

    LOG(Message) << "Action " << par().action << std::endl;
    auto &act = envGet(FMat, par().action);

    LOG(Message) << "Physical source " << par().source << std::endl;
    auto &phys_source = envGet(PropagatorField, par().source);

    if ((par().current != Current::Vector) && (par().current != Current::Axial))
    {
        HADRONS_ERROR(Argument, "par().current should either be Current::Vector or Current::Axial");
    }

    if (par().dir > 3) 
    {
        HADRONS_ERROR(Argument, "par().dir should be an unsigned int = 0, 1, 2, or 3");
    }

    Result result;
    result.gamma_src = par().gSrc;
    result.current = par().current;
    const int nt {env().getDim(Tp)};
    result.corr.resize(nt,0.);

    envGetTmp(PropagatorField, tmp);
    envGetTmp(ComplexField, tmp_current);
    SlicedComplex sumGC(nt);

    LOG(Message) << "Getting conserved current sink" << std::endl;    
    act.ContractConservedCurrent(prop, prop, tmp, phys_source, par().current, par().dir);

    // include phase
    if (!par().mom.empty()) 
    {
        LOG(Message) << "Projecting to momentum [" << par().mom << "]" << std::endl;

        auto &ph = envGet(LatticeComplex, momphName_);

        Complex           i(0.0,1.0);
        std::vector<Real> p;

        envGetTmp(LatticeComplex, coor);
        p  = strToVec<Real>(par().mom);
        ph = Zero();
        for (unsigned int mu = 0; mu < p.size(); mu++)
        {
            LatticeCoordinate(coor, mu);
            ph = ph + (p[mu]/env().getDim(mu))*coor;
        }
        ph = exp((Real)(2*M_PI)*i*ph);

        tmp = ph*tmp;
    }

    //  trace it out
    LOG(Message) << "Contracting with gSrc = " << par().gSrc << std::endl;
    tmp_current = trace(par().gSrc*tmp);
    SliceOut(result.corr, sumGC, tmp_current, false);

    saveResult(par().output, "conservedSink", result);
    auto &out = envGet(HadronsSerializable, getName());
    out = result;

}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_ContractConservedSink_hpp_
