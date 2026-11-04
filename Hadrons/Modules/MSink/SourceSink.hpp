/*
 * SourceSink.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2023
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Raoul Hodgson <raoul.hodgson@ed.ac.uk.com>
 * Author: Raoul Hodgson <raoul.hodgson@ed.ac.uk>
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

#ifndef Hadrons_MSink_SourceSink_hpp_
#define Hadrons_MSink_SourceSink_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                                 SourceSink                                 *
 ******************************************************************************/
/*
 * Sink function using a source field:
 *   sink(q)[t] = sum_{x, x_t = t} S(x) q(x)
 * with S = source       (sourceOp = "none"), or
 *      S = adj(source)  (sourceOp = "adj").
 * Use "adj" to close a propagator against a complex noise source (e.g.
 * MSource::Z2), so that the noise average gives eta(z) adj(eta)(z') -> delta.
 * sourceOp has no default value: any other value is an error.
 */
BEGIN_MODULE_NAMESPACE(MSink)

class SourceSinkPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(SourceSinkPar,
                                    std::string, source,
                                    std::string, sourceOp);
};

template <typename FImpl>
class TSourceSink: public Module<SourceSinkPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    SINK_TYPE_ALIASES();
public:
    // constructor
    TSourceSink(const std::string name);
    // destructor
    virtual ~TSourceSink(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
protected:
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(SourceSink,       TSourceSink<FIMPL>,        MSink);

/******************************************************************************
 *                          TSourceSink implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TSourceSink<FImpl>::TSourceSink(const std::string name)
: Module<SourceSinkPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TSourceSink<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().source};
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TSourceSink<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TSourceSink<FImpl>::setup(void)
{
    if ((par().sourceOp != "adj") and (par().sourceOp != "none"))
    {
        HADRONS_ERROR(Argument, "sourceOp must be 'adj' or 'none' (got '"
                      + par().sourceOp + "')");
    }
    envCreate(SinkFn, getName(), 1, nullptr);
}


// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TSourceSink<FImpl>::execute(void)
{
    LOG(Message) << "Setting up sink function with source '" << par().source
                 << "' as the sink (sourceOp '" << par().sourceOp << "')" << std::endl;

    PropagatorField &source  = envGet(PropagatorField, par().source);
    const bool      adjoint  = (par().sourceOp == "adj");
    
    auto sink = [this, source, adjoint](const PropagatorField &field)
    {
        SlicedPropagator res;
        PropagatorField  tmp(field.Grid());

        if (adjoint)
        {
            tmp = adj(source)*field;
        }
        else
        {
            tmp = source*field;
        }
        
        sliceSum(tmp, res, Tp);
        
        return res;
    };
    envGet(SinkFn, getName()) = sink;
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MSink_SourceSink_hpp_
