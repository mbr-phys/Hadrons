/*
 * Test_wilsonclover_solvers.cpp, part of Hadrons
 *
 * Compare Hadrons Wilson-clover fixed- and mixed-precision solver modules.
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
    LOG(Message) << "Grid initialized" << std::endl;

    Application application;

    Application::GlobalPar globalPar;
    globalPar.trajCounter.start = 0;
    globalPar.trajCounter.end   = 1;
    globalPar.trajCounter.step  = 1;
    globalPar.runId             = "wilsonclover_solvers";
    application.setPar(globalPar);

    application.createModule<MGauge::Unit>("gauge");

    MSource::Point::Par sourcePar;
    sourcePar.position = "0 0 0 0";
    application.createModule<MSource::Point>("source", sourcePar);

    MAction::WilsonClover::Par actionPar;
    actionPar.gauge    = "gauge";
    actionPar.mass     = 0.5;
    actionPar.csw_r    = 1.0;
    actionPar.csw_t    = 1.0;
    actionPar.cF       = 0.0;
    actionPar.boundary = "1 1 1 1";
    actionPar.twist    = "0 0 0 0";
    application.createModule<MAction::WilsonClover>("WilsonClover", actionPar);

#ifdef GRID_DEFAULT_PRECISION_DOUBLE
    MUtilities::GaugeSinglePrecisionCast::Par gaugeCastPar;
    gaugeCastPar.field = "gauge";
    application.createModule<MUtilities::GaugeSinglePrecisionCast>(
        "gauge_f", gaugeCastPar);

    actionPar.gauge = "gauge_f";
    application.createModule<MAction::WilsonCloverF>(
        "WilsonCloverF", actionPar);
#endif

    MSolver::RBPrecBiCGSTAB::Par bicgstabPar;
    bicgstabPar.action       = "WilsonClover";
    bicgstabPar.maxIteration = 2000;
    bicgstabPar.residual     = 1.0e-10;
    bicgstabPar.guesser      = "";
    application.createModule<MSolver::RBPrecBiCGSTAB>("BiCGSTAB", bicgstabPar);

    MSolver::RBPrecFGMRES::Par fgmresPar;
    fgmresPar.action         = "WilsonClover";
    fgmresPar.maxIteration   = 2000;
    fgmresPar.restartLength  = 25;
    fgmresPar.residual       = 1.0e-10;
    fgmresPar.guesser        = "";
    application.createModule<MSolver::RBPrecFGMRES>("FGMRES", fgmresPar);

#ifdef GRID_DEFAULT_PRECISION_DOUBLE
    MSolver::MixedPrecisionRBPrecFGMRES::Par mixedFgmresPar;
    mixedFgmresPar.innerAction         = "WilsonCloverF";
    mixedFgmresPar.outerAction         = "WilsonClover";
    mixedFgmresPar.maxInnerIteration   = 50;
    mixedFgmresPar.maxOuterIteration   = 2000;
    mixedFgmresPar.restartLength       = 25;
    mixedFgmresPar.innerResidual       = 1.0e-2;
    mixedFgmresPar.residual            = 1.0e-10;
    mixedFgmresPar.outerGuesser        = "";
    application.createModule<MSolver::MixedPrecisionRBPrecFGMRES>(
        "MixedFGMRES", mixedFgmresPar);
#endif

    MFermion::GaugeProp::Par propagatorPar;
    propagatorPar.source = "source";
    propagatorPar.solver = "BiCGSTAB";
    application.createModule<MFermion::GaugeProp>("prop_bicgstab", propagatorPar);

    propagatorPar.solver = "FGMRES";
    application.createModule<MFermion::GaugeProp>("prop_fgmres", propagatorPar);

#ifdef GRID_DEFAULT_PRECISION_DOUBLE
    propagatorPar.solver = "MixedFGMRES";
    application.createModule<MFermion::GaugeProp>(
        "prop_mixed_fgmres", propagatorPar);
#endif

    MUtilities::NormCheckPropagator::Par normCheckPar;
    normCheckPar.field     = "prop_fgmres";
    normCheckPar.reference = "prop_bicgstab";
    normCheckPar.output    = "";
    application.createModule<MUtilities::NormCheckPropagator>(
        "solver_comparison", normCheckPar);

#ifdef GRID_DEFAULT_PRECISION_DOUBLE
    normCheckPar.field = "prop_mixed_fgmres";
    application.createModule<MUtilities::NormCheckPropagator>(
        "mixed_solver_comparison", normCheckPar);
#endif

    application.saveParameterFile("TestWilsonCloverSolvers.xml");
    application.run();

    LOG(Message) << "Grid is finalizing now" << std::endl;
    Grid_finalize();

    return EXIT_SUCCESS;
}
