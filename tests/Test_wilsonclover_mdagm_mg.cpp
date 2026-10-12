/*
 * Demonstrate the Hadrons full-field MdagM Wilson-clover multigrid solver.
 * Example: ./Test_wilsonclover_mdagm_mg --grid 4.4.4.4 --mpi 1.1.1.1
 */
#include <Hadrons/Application.hpp>
#include <Hadrons/Modules.hpp>

using namespace Grid;
using namespace Hadrons;

int main(int argc, char *argv[])
{
    Grid_init(&argc, &argv);
    HadronsLogError.Active(GridLogError.isActive());
    HadronsLogWarning.Active(GridLogWarning.isActive());
    HadronsLogMessage.Active(GridLogMessage.isActive());
    HadronsLogIterative.Active(GridLogIterative.isActive());
    HadronsLogDebug.Active(GridLogDebug.isActive());

    Application application;
    Application::GlobalPar globalPar;
    globalPar.trajCounter.start = 0;
    globalPar.trajCounter.end   = 1;
    globalPar.trajCounter.step  = 1;
    globalPar.runId             = "wilsonclover_mdagm_mg";
    application.setPar(globalPar);

    application.createModule<MGauge::Unit>("gauge");
    MSource::Point::Par sourcePar;
    sourcePar.position = "0 0 0 0";
    application.createModule<MSource::Point>("source", sourcePar);

    MAction::WilsonClover::Par actionPar;
    actionPar.gauge = "gauge";
    actionPar.mass = 0.5;
    actionPar.csw_r = 1.0;
    actionPar.csw_t = 1.0;
    actionPar.cF = 1.0;
    actionPar.boundary = "1 1 1 1";
    actionPar.twist = "0 0 0 0";
    application.createModule<MAction::WilsonClover>("WilsonClover", actionPar);

    MSolver::WilsonCloverMdagMMGSetup::Par setupPar;
    setupPar.action = "WilsonClover";
    setupPar.mg = Grid::WilsonCloverMdagMMG::Params(
        2, {{2, 2, 2, 2}}, {1.0e-8}, {2}, {4}, false,
        {1.0e-1}, {2}, {4}, 5.0e-2, 10, 50);
    setupPar.runChecks = false;
    setupPar.checkTolerance = 1.0e-11;
    application.createModule<MSolver::WilsonCloverMdagMMGSetup>(
        "MGSetup", setupPar);

    MSolver::WilsonCloverMdagMMGFGMRES::Par mgSolverPar;
    mgSolverPar.action = "WilsonClover";
    mgSolverPar.mgSetup = "MGSetup";
    mgSolverPar.maxIteration = 500;
    mgSolverPar.restartLength = 25;
    mgSolverPar.residual = 1.0e-10;
    mgSolverPar.guesser = "";
    application.createModule<MSolver::WilsonCloverMdagMMGFGMRES>(
        "MGFGMRES", mgSolverPar);

    MSolver::RBPrecBiCGSTAB::Par referencePar;
    referencePar.action = "WilsonClover";
    referencePar.maxIteration = 2000;
    referencePar.residual = 1.0e-10;
    referencePar.guesser = "";
    application.createModule<MSolver::RBPrecBiCGSTAB>(
        "BiCGSTAB", referencePar);

    MFermion::GaugeProp::Par propPar;
    propPar.source = "source";
    propPar.solver = "MGFGMRES";
    application.createModule<MFermion::GaugeProp>("prop_mg", propPar);
    propPar.solver = "BiCGSTAB";
    application.createModule<MFermion::GaugeProp>("prop_reference", propPar);

    MUtilities::NormCheckPropagator::Par checkPar;
    checkPar.field = "prop_mg";
    checkPar.reference = "prop_reference";
    checkPar.output = "";
    application.createModule<MUtilities::NormCheckPropagator>(
        "mg_comparison", checkPar);

    application.saveParameterFile("TestWilsonCloverMdagMMG.xml");
    application.run();
    Grid_finalize();
    return EXIT_SUCCESS;
}
