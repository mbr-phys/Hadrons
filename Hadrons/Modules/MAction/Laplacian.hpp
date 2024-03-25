/*
 * Laplacian.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2023
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
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
#ifndef Hadrons_MAction_Laplacian_hpp_
#define Hadrons_MAction_Laplacian_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                           Laplacian operator                               *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MAction)

class LaplacianPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(LaplacianPar,
                                    std::string, gauge,
                                    double,      m2);
};

template <typename Field>
class Laplacian: public LinearOperatorBase<Field>
{
public:
    Laplacian(const double m2)
    : m2_(m2)
    {}

    virtual void OpDiag(const Field &in, Field &out)
    {
        HADRONS_ERROR(Implementation, "Laplacian method not implemented");
    }

    virtual void OpDir(const Field &in, Field &out, int dir, int disp)
    {
        HADRONS_ERROR(Implementation, "Laplacian method not implemented");
    }

    virtual void OpDirAll(const Field &in, std::vector<Field> &out)
    {
        HADRONS_ERROR(Implementation, "Laplacian method not implemented");
    }

    virtual void Op(const Field &in, Field &out)
    {
        unsigned int nd = in.Grid()->Nd();

        out = (m2_ + 2.*nd)*in;
        for (unsigned int mu = 0; mu < nd; ++mu)
        {
            out  -= Cshift(in, mu, 1);
            out  -= Cshift(in, mu, -1);
        }
    }

    virtual void AdjOp(const Field &in, Field &out)
    {
        Op(in, out);
    }

    virtual void HermOpAndNorm(const Field &in, Field &out, RealD &n1, RealD &n2)
    {
        HermOp(in,out);
        Complex dot = innerProduct(in,out); 
        n1 = real(dot);
        n2 = norm2(out);
    }

    virtual void HermOp(const Field &in, Field &out)
    {
        Op(in, out);
    }
private:
    double m2_;
};

template <typename Field>
class TLaplacian: public Module<LaplacianPar>
{
public:
    // constructor
    TLaplacian(const std::string name);
    // destructor
    virtual ~TLaplacian(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(FermionLaplacian, TLaplacian<FIMPL::FermionField>, MAction);

/******************************************************************************
 *                      TLaplacian implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename Field>
TLaplacian<Field>::TLaplacian(const std::string name)
: Module<LaplacianPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename Field>
std::vector<std::string> TLaplacian<Field>::getInput(void)
{
    std::vector<std::string> in = {par().gauge};
    
    return in;
}

template <typename Field>
std::vector<std::string> TLaplacian<Field>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename Field>
void TLaplacian<Field>::setup(void)
{
    LOG(Message) << "Setting up Laplacian operator with m^2= " << par().m2 << std::endl;
    envCreateDerived(LinearOperatorBase<Field>, ARG(Laplacian<Field>), 
                     getName(), 1, par().m2);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename Field>
void TLaplacian<Field>::execute(void)
{}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MAction_Laplacian_hpp_
