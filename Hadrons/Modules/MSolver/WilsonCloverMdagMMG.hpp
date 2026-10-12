/*
 * WilsonCloverMdagMMG.hpp, part of Hadrons
 *
 * Copyright (C) 2026
 */
#ifndef Hadrons_MSolver_WilsonCloverMdagMMG_hpp_
#define Hadrons_MSolver_WilsonCloverMdagMMG_hpp_

#include <cmath>
#include <memory>

#include <Hadrons/Global.hpp>
#include <Grid/algorithms/multigrid/WilsonCloverMdagMMultiGrid.h>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Solver.hpp>

BEGIN_HADRONS_NAMESPACE
BEGIN_MODULE_NAMESPACE(MSolver)

class WilsonCloverMdagMMGSetupData
{
public:
    FERM_TYPE_ALIASES(FIMPLD,);
    static constexpr int nBasis = 40;
    typedef Grid::WilsonCloverMdagMMG::PreconditionerBase<FermionField>
        Preconditioner;

    WilsonCloverMdagMMGSetupData(
        FMat &matrix, GridCartesian *grid,
        const Grid::WilsonCloverMdagMMG::Params &params,
        const bool runChecks, const double checkTolerance)
    : matrix_(matrix)
    , params_(params)
    , levels_(grid, params_)
    , preconditioner_(Grid::WilsonCloverMdagMMG::create<
          vSpinColourVectorD, vTComplexD, nBasis>(
              params_, levels_, matrix_, matrix_))
    {
        preconditioner_->setup();
        if (runChecks)
        {
            preconditioner_->runChecks(checkTolerance);
        }
    }

    FMat &getMatrix(void) { return matrix_; }
    Preconditioner &getPreconditioner(void) { return *preconditioner_; }
    void reportTimings(void) { preconditioner_->reportTimings(); }
private:
    FMat                                      &matrix_;
    Grid::WilsonCloverMdagMMG::Params          params_;
    Grid::WilsonCloverMdagMMG::LevelInfo       levels_;
    std::unique_ptr<Preconditioner>             preconditioner_;
};

class WilsonCloverMdagMMGSetupPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(
        WilsonCloverMdagMMGSetupPar,
        std::string, action,
        Grid::WilsonCloverMdagMMG::Params, mg,
        bool, runChecks,
        double, checkTolerance);
};

class TWilsonCloverMdagMMGSetup:
    public Module<WilsonCloverMdagMMGSetupPar>
{
public:
    FERM_TYPE_ALIASES(FIMPLD,);
    TWilsonCloverMdagMMGSetup(const std::string name): Module(name) {}
    virtual std::vector<std::string> getInput(void) { return {par().action}; }
    virtual std::vector<std::string> getOutput(void) { return {getName()}; }
    virtual DependencyMap getObjectDependencies(void)
    {
        return {{par().action, getName()}};
    }
protected:
    virtual void setup(void)
    {
        if ((par().mg.nLevels < 2) || (par().mg.nLevels > 4))
        {
            HADRONS_ERROR(Argument, "multigrid supports two to four levels");
        }
        if (!std::isfinite(par().checkTolerance) ||
            (par().checkTolerance <= 0.0))
        {
            HADRONS_ERROR(Argument,
                          "multigrid check tolerance must be finite and positive");
        }
        Grid::WilsonCloverMdagMMG::validateParameters(par().mg);
        auto &matrix = envGetDerived(
            FMat, CompactWilsonClover<FIMPLD>, par().action);
        envCreate(WilsonCloverMdagMMGSetupData, getName(), 1, matrix,
                  envGetGrid(FermionField), par().mg,
                  par().runChecks, par().checkTolerance);
    }
    virtual void execute(void) {}
};

MODULE_REGISTER_TMP(WilsonCloverMdagMMGSetup,
                    TWilsonCloverMdagMMGSetup, MSolver);

class WilsonCloverMdagMMGFGMRESPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(
        WilsonCloverMdagMMGFGMRESPar,
        std::string, action,
        std::string, mgSetup,
        unsigned int, maxIteration,
        unsigned int, restartLength,
        double, residual,
        std::string, guesser);
};

class TWilsonCloverMdagMMGFGMRES:
    public Module<WilsonCloverMdagMMGFGMRESPar>
{
public:
    FERM_TYPE_ALIASES(FIMPLD,);
    SOLVER_TYPE_ALIASES(FIMPLD,);
    TWilsonCloverMdagMMGFGMRES(const std::string name): Module(name) {}
    virtual std::vector<std::string> getInput(void)
    {
        std::vector<std::string> in = {par().action, par().mgSetup};
        if (!par().guesser.empty()) in.push_back(par().guesser);
        return in;
    }
    virtual std::vector<std::string> getOutput(void)
    {
        return {getName(), getName() + "_subtract"};
    }
    virtual DependencyMap getObjectDependencies(void)
    {
        DependencyMap dep;
        for (auto output: getOutput())
        {
            dep.insert({par().action, output});
            dep.insert({par().mgSetup, output});
            if (!par().guesser.empty()) dep.insert({par().guesser, output});
        }
        return dep;
    }
protected:
    virtual void setup(void);
    virtual void execute(void) {}
};

MODULE_REGISTER_TMP(WilsonCloverMdagMMGFGMRES,
                    TWilsonCloverMdagMMGFGMRES, MSolver);

inline void TWilsonCloverMdagMMGFGMRES::setup(void)
{
    if (par().maxIteration == 0)
        HADRONS_ERROR(Argument, "zero maximum iteration");
    if ((par().restartLength == 0) ||
        (par().restartLength > par().maxIteration))
        HADRONS_ERROR(Argument, "invalid FGMRES restart length");
    if (!std::isfinite(par().residual) || (par().residual <= 0.0))
        HADRONS_ERROR(Argument, "residual must be finite and positive");

    auto &matrix = envGet(FMat, par().action);
    auto &mg = envGet(WilsonCloverMdagMMGSetupData, par().mgSetup);
    if (&matrix != &mg.getMatrix())
        HADRONS_ERROR(Definition,
                      "multigrid setup was constructed for a different action");

    LinearFunction<FermionField> *guesserPt = nullptr;
    if (!par().guesser.empty())
        guesserPt = &envGet(LinearFunction<FermionField>, par().guesser);

    auto makeSolver = [&matrix, &mg, guesserPt, this](bool subtractGuess)
    {
        return [&matrix, &mg, guesserPt, subtractGuess, this]
        (FermionField &solution, const FermionField &source)
        {
            ZeroGuesser<FermionField> zeroGuesser;
            auto &guesser = (guesserPt == nullptr) ?
                static_cast<LinearFunction<FermionField> &>(zeroGuesser) :
                *guesserPt;
            FermionField initialGuess(source.Grid());
            FermionField normalSource(source.Grid());
            FermionField residual(source.Grid());

            guesser(source, solution);
            initialGuess = solution;
            matrix.Mdag(source, normalSource);
            MdagMLinearOperator<FMat, FermionField> normalOperator(matrix);
            FlexibleGeneralisedMinimalResidual<FermionField> fgmres(
                par().residual, par().maxIteration,
                mg.getPreconditioner(), par().restartLength);
            fgmres(normalOperator, normalSource, solution);
            if (subtractGuess) solution = solution - initialGuess;

            if (!subtractGuess)
            {
                matrix.M(solution, residual);
                residual = residual - source;
                LOG(Message) << "Wilson-clover MdagM MG-FGMRES full D residual "
                             << std::sqrt(norm2(residual)/norm2(source))
                             << std::endl;
            }
            LOG(Message) << "Wilson-clover MdagM MG-FGMRES outer iterations "
                         << fgmres.IterationCount << std::endl;
            mg.reportTimings();
        };
    };

    auto solver = makeSolver(false);
    envCreate(Solver, getName(), 1, solver, matrix);
    auto solverSubtract = makeSolver(true);
    envCreate(Solver, getName() + "_subtract", 1, solverSubtract, matrix);
}

END_MODULE_NAMESPACE
END_HADRONS_NAMESPACE
#endif
