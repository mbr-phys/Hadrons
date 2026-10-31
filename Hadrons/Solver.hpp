/*
 * Solver.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
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
#ifndef Hadrons_Solver_hpp_
#define Hadrons_Solver_hpp_

#include <Hadrons/Global.hpp>

BEGIN_HADRONS_NAMESPACE

template <typename FImpl>
class Solver
{
public:
    typedef typename FImpl::FermionField                      FermionField;
    typedef FermionOperator<FImpl>                            FMat; 
    typedef std::function<void(FermionField &, 
                               const FermionField &)>         SolverFn;
    typedef std::function<void(std::vector<FermionField> &, 
                               const std::vector<FermionField> &)> SolverVFn;
public:
    Solver(SolverFn fn, FMat &mat): mat_(mat), fn_(fn), vfn_(nullptr) {}
    Solver(SolverFn fn, SolverVFn vfn, FMat &mat): mat_(mat), fn_(fn), vfn_(vfn) {}

    // Wilson-clover with fixed temporal boundaries (time phase 0) has zero rows
    // on t = 0 and t = T-1, so source components there can never be solved
    // for: they are projected out of the source, and out of the solution.
    void operator()(FermionField &sol, const FermionField &src)
    {
        auto project = boundaryProjector();

        if (!project)
        {
            fn_(sol, src);
        }
        else
        {
            FermionField projSrc(src);

            projectSource(project, projSrc);
            fn_(sol, projSrc);
            project(sol);
        }
    }

    void operator()(std::vector<FermionField> &sol, const std::vector<FermionField> &src)
    {
        if (sol.size() != src.size())
        {
            HADRONS_ERROR(Size, "source and solution vectors size mismatch");
        }
        auto project = boundaryProjector();

        if (vfn_ == nullptr)
        {
            for (unsigned int i = 0; i < sol.size(); ++i)
            {
                (*this)(sol[i], src[i]);
            }
        }
        else if (!project)
        {
            vfn_(sol, src);
        }
        else
        {
            std::vector<FermionField> projSrc(src);

            for (auto &s: projSrc)
            {
                projectSource(project, s);
            }
            vfn_(sol, projSrc);
            for (auto &s: sol)
            {
                project(s);
            }
        }
    }

    FMat & getFMat(void)
    {
        return mat_;
    }
private:
    typedef std::function<void(FermionField &)> Projector;

    // empty unless mat_ is a compact Wilson(-exp)-clover operator with fixed
    // temporal boundaries
    Projector boundaryProjector(void)
    {
        if constexpr (std::is_same<FImpl, WilsonImplF>::value
                      or std::is_same<FImpl, WilsonImplD>::value)
        {
            if (auto *m = dynamic_cast<CompactWilsonClover<FImpl> *>(&mat_))
            {
                if (m->fixedBoundaries)
                {
                    return [m](FermionField &f){ m->ApplyBoundaryProjection(f); };
                }
            }
            if (auto *m = dynamic_cast<CompactWilsonExpClover<FImpl> *>(&mat_))
            {
                if (m->fixedBoundaries)
                {
                    return [m](FermionField &f){ m->ApplyBoundaryProjection(f); };
                }
            }
        }

        return nullptr;
    }

    static void projectSource(const Projector &project, FermionField &src)
    {
        RealD before = norm2(src), removed;

        project(src);
        removed = before - norm2(src);
        if (removed > 0.)
        {
            LOG(Message) << "Fixed temporal boundaries: removed "
                         << std::sqrt(removed/before)
                         << " (relative norm) of the source on the boundary slices"
                         << std::endl;
        }
    }

    FMat      &mat_;
    SolverFn  fn_;
    SolverVFn vfn_;
};

END_HADRONS_NAMESPACE

#endif // Hadrons_Solver_hpp_
