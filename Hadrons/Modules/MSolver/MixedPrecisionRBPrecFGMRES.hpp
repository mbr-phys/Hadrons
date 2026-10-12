/*
 * MixedPrecisionRBPrecFGMRES.hpp, part of Hadrons
 *
 * Copyright (C) 2026
 *
 * Hadrons is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 */
#ifndef Hadrons_MSolver_MixedPrecisionRBPrecFGMRES_hpp_
#define Hadrons_MSolver_MixedPrecisionRBPrecFGMRES_hpp_

#include <algorithm>
#include <cmath>
#include <cstdint>

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Solver.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *       Mixed-precision Schur red-black preconditioned flexible GMRES        *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MSolver)

class MixedPrecisionRBPrecFGMRESPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(MixedPrecisionRBPrecFGMRESPar,
                                    std::string , innerAction,
                                    std::string , outerAction,
                                    unsigned int, maxInnerIteration,
                                    unsigned int, maxOuterIteration,
                                    unsigned int, restartLength,
                                    double      , innerResidual,
                                    double      , residual,
                                    std::string , outerGuesser);
};

template <typename FImplInner, typename FImplOuter>
class TMixedPrecisionRBPrecFGMRES:
    public Module<MixedPrecisionRBPrecFGMRESPar>
{
public:
    FERM_TYPE_ALIASES(FImplInner, Inner);
    FERM_TYPE_ALIASES(FImplOuter, Outer);
    SOLVER_TYPE_ALIASES(FImplOuter,);
    typedef FNonHermitianSchurOpInner SchurFMatInner;
private:
    class SchurBiCGSTABPreconditioner:
        public LinearFunction<FermionFieldInner>
    {
    public:
        using LinearFunction<FermionFieldInner>::operator();

        SchurBiCGSTABPreconditioner(SchurFMatInner &op,
                                    const double residual,
                                    const unsigned int maxIteration)
        : op_(op)
        , residual_(residual)
        , maxIteration_(maxIteration)
        {}

        virtual void operator()(const FermionFieldInner &source,
                                FermionFieldInner &solution)
        {
            solution = Zero();
            BiCGSTAB<FermionFieldInner> solver(residual_, maxIteration_, false);
            solver(op_, source, solution);

            ++applicationCount_;
            const auto iterations = std::min<unsigned long long>(
                static_cast<unsigned long long>(solver.IterationsToComplete),
                static_cast<unsigned long long>(maxIteration_));
            iterationCount_ += iterations;

            const double solutionNorm2 = norm2(solution);
            if (!std::isfinite(solutionNorm2))
            {
                HADRONS_ERROR(Program,
                    "single-precision Schur BiCGSTAB preconditioner produced "
                    "a non-finite correction");
            }
        }

        std::uint64_t getApplicationCount(void) const
        {
            return applicationCount_;
        }

        std::uint64_t getIterationCount(void) const
        {
            return iterationCount_;
        }
    private:
        SchurFMatInner &op_;
        double          residual_;
        unsigned int    maxIteration_;
        std::uint64_t   applicationCount_{0};
        std::uint64_t   iterationCount_{0};
    };
public:
    TMixedPrecisionRBPrecFGMRES(const std::string name);
    virtual ~TMixedPrecisionRBPrecFGMRES(void) = default;
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    virtual DependencyMap getObjectDependencies(void);
protected:
    virtual void setup(void);
    virtual void execute(void);
};

MODULE_REGISTER_TMP(MixedPrecisionRBPrecFGMRES,
                    ARG(TMixedPrecisionRBPrecFGMRES<FIMPLF, FIMPLD>),
                    MSolver);
MODULE_REGISTER_TMP(ZMixedPrecisionRBPrecFGMRES,
                    ARG(TMixedPrecisionRBPrecFGMRES<ZFIMPLF, ZFIMPLD>),
                    MSolver);

/******************************************************************************
 *              TMixedPrecisionRBPrecFGMRES implementation                    *
 ******************************************************************************/
template <typename FImplInner, typename FImplOuter>
TMixedPrecisionRBPrecFGMRES<FImplInner, FImplOuter>
::TMixedPrecisionRBPrecFGMRES(const std::string name)
: Module<MixedPrecisionRBPrecFGMRESPar>(name)
{}

template <typename FImplInner, typename FImplOuter>
std::vector<std::string>
TMixedPrecisionRBPrecFGMRES<FImplInner, FImplOuter>::getInput(void)
{
    std::vector<std::string> in = {par().innerAction, par().outerAction};

    if (!par().outerGuesser.empty())
    {
        in.push_back(par().outerGuesser);
    }

    return in;
}

template <typename FImplInner, typename FImplOuter>
std::vector<std::string>
TMixedPrecisionRBPrecFGMRES<FImplInner, FImplOuter>::getOutput(void)
{
    return {getName(), getName() + "_subtract"};
}

template <typename FImplInner, typename FImplOuter>
DependencyMap
TMixedPrecisionRBPrecFGMRES<FImplInner, FImplOuter>
::getObjectDependencies(void)
{
    DependencyMap dep;

    dep.insert({par().innerAction, getName()});
    dep.insert({par().innerAction, getName() + "_subtract"});
    dep.insert({par().outerAction, getName()});
    dep.insert({par().outerAction, getName() + "_subtract"});
    if (!par().outerGuesser.empty())
    {
        dep.insert({par().outerGuesser, getName()});
        dep.insert({par().outerGuesser, getName() + "_subtract"});
    }

    return dep;
}

#define SOLVER_BODY                                                                          \
ZeroGuesser<FermionFieldOuter> defaultGuesser;                                               \
LinearFunction<FermionFieldOuter> &guesser =                                                 \
    (guesserPt == nullptr) ? defaultGuesser : *guesserPt;                                    \
SchurFMatInner innerSchur(imat);                                                             \
SchurBiCGSTABPreconditioner innerPreconditioner(                                             \
    innerSchur, par().innerResidual, par().maxInnerIteration);                               \
MixedPrecisionFlexibleGeneralisedMinimalResidual<                                            \
    FermionFieldOuter, FermionFieldInner> fgmres(                                            \
        par().residual, par().maxOuterIteration, innerGrid,                                  \
        innerPreconditioner, par().restartLength);                                           \
HADRONS_DEFAULT_NON_HERMITIAN_SCHUR_SOLVE<FermionFieldOuter> schurSolver(fgmres);            \
schurSolver.subtractGuess(subGuess);                                                         \
schurSolver(omat, source, sol, guesser);                                                     \
LOG(Message) << "Mixed-precision FGMRES completed with "                                    \
             << fgmres.IterationCount << " outer iterations and "                           \
             << innerPreconditioner.getIterationCount()                                     \
             << " accumulated inner iterations in "                                         \
             << innerPreconditioner.getApplicationCount()                                   \
             << " preconditioner applications" << std::endl;

template <typename FImplInner, typename FImplOuter>
void TMixedPrecisionRBPrecFGMRES<FImplInner, FImplOuter>::setup(void)
{
    if (par().maxInnerIteration == 0)
    {
        HADRONS_ERROR(Argument, "zero maximum inner iteration");
    }
    if (par().maxOuterIteration == 0)
    {
        HADRONS_ERROR(Argument, "zero maximum outer iteration");
    }
    if (par().restartLength == 0)
    {
        HADRONS_ERROR(Argument, "zero FGMRES restart length");
    }
    if (par().restartLength > par().maxOuterIteration)
    {
        HADRONS_ERROR(Argument,
                      "FGMRES restart length exceeds maximum outer iteration");
    }
    if (!std::isfinite(par().innerResidual) || (par().innerResidual < 0.0))
    {
        HADRONS_ERROR(Argument,
                      "inner residual must be finite and non-negative");
    }
    if (!std::isfinite(par().residual) || (par().residual <= 0.0))
    {
        HADRONS_ERROR(Argument, "outer residual must be finite and positive");
    }

    const auto innerLs = env().getObjectLs(par().innerAction);
    const auto outerLs = env().getObjectLs(par().outerAction);
    if (innerLs != outerLs)
    {
        HADRONS_ERROR(Size,
                      "inner and outer actions have different fifth dimensions");
    }

    auto &imat = envGet(FMatInner, par().innerAction);
    auto &omat = envGet(FMatOuter, par().outerAction);
    auto *innerGrid = getGrid<FermionFieldInner>(true, innerLs);
    auto *outerGrid = getGrid<FermionFieldOuter>(true, outerLs);
    if (innerGrid->GlobalDimensions() != outerGrid->GlobalDimensions())
    {
        HADRONS_ERROR(Size,
                      "inner and outer actions have different global dimensions");
    }

    const std::string innerPolicy = (par().innerResidual == 0.0)
        ? "fixed-work inner mode"
        : "inner residual " + std::to_string(par().innerResidual);
    LOG(Message) << "Setting up mixed-precision non-Hermitian Schur FGMRES for "
                 << "inner/outer action '" << par().innerAction << "'/'"
                 << par().outerAction << "', outer residual " << par().residual
                 << ", maximum outer iteration " << par().maxOuterIteration
                 << ", restart length " << par().restartLength
                 << ", maximum inner iteration " << par().maxInnerIteration
                 << ", " << innerPolicy << std::endl;

    LinearFunction<FermionFieldOuter> *guesserPt = nullptr;
    if (!par().outerGuesser.empty())
    {
        guesserPt = &envGet(LinearFunction<FermionFieldOuter>,
                            par().outerGuesser);
    }

    auto makeSolver = [&imat, &omat, innerGrid, guesserPt, this](bool subGuess)
    {
        return [&imat, &omat, innerGrid, guesserPt, subGuess, this]
        (FermionFieldOuter &sol, const FermionFieldOuter &source)
        {
            SOLVER_BODY;
        };
    };

    auto solver = makeSolver(false);
    envCreate(Solver, getName(), outerLs, solver, omat);
    auto solverSubtract = makeSolver(true);
    envCreate(Solver, getName() + "_subtract", outerLs, solverSubtract, omat);
}

#undef SOLVER_BODY

template <typename FImplInner, typename FImplOuter>
void TMixedPrecisionRBPrecFGMRES<FImplInner, FImplOuter>::execute(void)
{}

END_MODULE_NAMESPACE
END_HADRONS_NAMESPACE

#endif // Hadrons_MSolver_MixedPrecisionRBPrecFGMRES_hpp_
