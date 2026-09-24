#include <Hadrons/Application.hpp>
#include <Hadrons/Modules.hpp>

using namespace Grid;
using namespace Hadrons;

int main(int argc, char *argv[])
{
    // initialization
    Grid_init(&argc, &argv);
    HadronsLogError.Active(GridLogError.isActive());
    HadronsLogWarning.Active(GridLogWarning.isActive());
    HadronsLogMessage.Active(GridLogMessage.isActive());
    HadronsLogIterative.Active(GridLogIterative.isActive());
    HadronsLogDebug.Active(GridLogDebug.isActive());
    LOG(Message) << "Grid initialized" << std::endl;
    
    // run setup
    Application application;
    
    // global parameters
    Application::GlobalPar globalPar;
    globalPar.trajCounter.start    = 1500;
    globalPar.trajCounter.end      = 1520;
    globalPar.trajCounter.step     = 20;
    globalPar.runId                = "test";
    globalPar.genetic.maxGen       = 1000;
    globalPar.genetic.maxCstGen    = 200;
    globalPar.genetic.popSize      = 20;
    globalPar.genetic.mutationRate = 0.1;
    application.setPar(globalPar);
    
    // gauge field
    application.createModule<MGauge::Random>("gauge");

    // DWF action and solver at flow time zero. Both condensate paths use
    // propagators obtained by inverting this same Dirac operator.
    std::string boundary = "1 1 1 -1";
    std::string twist = "0. 0. 0. 0.";

    MAction::DWF::Par actionPar;
    actionPar.gauge = "gauge";
    actionPar.Ls = 12;
    actionPar.M5 = 1.8;
    actionPar.mass = 0.25;
    actionPar.boundary = boundary;
    actionPar.twist = twist;
    application.createModule<MAction::DWF>("DWF", actionPar);

    MSolver::RBPrecCG::Par solverPar;
    solverPar.action = "DWF";
    solverPar.residual = 1.0e-8;
    solverPar.maxIteration = 10000;
    application.createModule<MSolver::RBPrecCG>("CG", solverPar);
    
    // Generate independent Z2 stochastic PropagatorField sources.
    MSource::Z2::Par noisePar;
    noisePar.tA = 0;
    noisePar.tB = GridDefaultLatt()[Tp] - 1;  // full time extent
    application.createModule<MSource::Z2>("eta", noisePar);
    application.createModule<MSource::Z2>("chi", noisePar);
    
    // Pre-compute gauge trajectory (forward flow)
    // This saves gauge fields at all flow times for the adjoint flow to use
    MGradientFlow::WilsonFlow::Par gaugeFlowPar;
    gaugeFlowPar.gauge = "gauge";
    gaugeFlowPar.steps = 1;
    gaugeFlowPar.step_size = 0.01;
    gaugeFlowPar.meas_interval = 1;  // save every step
    gaugeFlowPar.save_history = true;
    gaugeFlowPar.save_rk_stages = true;
    application.createModule<MGradientFlow::WilsonFlow>("WilsonFlow", gaugeFlowPar);
    // Outputs include WilsonFlow_U_t0.00, WilsonFlow_U_t0.01, and the
    // W1/W2 stages labelled by the start time of their RK step.
    
    // Adjoint-flow the noise backward: eta(tau=0.01) -> xi(tau=0.00)
    // Reuse the saved RK stages instead of reconstructing them from U(0).
    MGradientFlow::AdjointFermionFlow::Par adjointPar;
    adjointPar.gauge = "WilsonFlow_U_t0.00";     // U at earlier time (tau=0.00)
    adjointPar.stage1 = "WilsonFlow_W1_t0.00";
    adjointPar.stage2 = "WilsonFlow_W2_t0.00";
    adjointPar.sources = {"eta"};                // xi at s+eps (treating eta as if at tau=0.01)
    adjointPar.sourceTypes = {"PropagatorField"};
    adjointPar.outSources = {"xi_t0.00"};        // xi at s (after flowing to tau=0.00)
    adjointPar.steps = 1;                        // single step
    adjointPar.step_size = 0.01;
    adjointPar.bc = -1;                          // antiperiodic time
    application.createModule<MGradientFlow::AdjointFermionFlow>("AdjointFlow", adjointPar);

    // Keep the reconstruction path as a direct diagnostic that saved raw
    // stages reproduce the existing adjoint-flow implementation.
    MGradientFlow::AdjointFermionFlow::Par reconstructedAdjointPar = adjointPar;
    reconstructedAdjointPar.stage1.clear();
    reconstructedAdjointPar.stage2.clear();
    reconstructedAdjointPar.outSources = {"xi_reconstructed_t0.00"};
    application.createModule<MGradientFlow::AdjointFermionFlow>(
        "AdjointFlowReconstructed", reconstructedAdjointPar);

    // Invert the flow-zero Dirac operator on the back-flowed source.
    MFermion::GaugeProp::Par adjointPropagatorPar;
    adjointPropagatorPar.solver = "CG";
    adjointPropagatorPar.source = "xi_t0.00";
    application.createModule<MFermion::GaugeProp>("xi_propagator_t0.00",
                                                   adjointPropagatorPar);
    //
    // First invert D chi_propagator = chi, then flow both the stochastic
    // source and its propagator with the shared gauge trajectory.
    MFermion::GaugeProp::Par positivePropagatorPar;
    positivePropagatorPar.solver = "CG";
    positivePropagatorPar.source = "chi";
    application.createModule<MFermion::GaugeProp>("chi_propagator",positivePropagatorPar);

    // Forward-flow chi from tau=0 to tau=0.01. Together with the adjoint
    // path, this supplies the two inner products in the discrete identity:
    // (xi_0, chi_0) = (eta_tau, chi_tau).
    
    MGradientFlow::FermionFlow::Par forwardPar;
    forwardPar.gauge = "gauge";
    forwardPar.steps = 1;
    forwardPar.step_size = 0.01;
    forwardPar.meas_interval = 1;
    forwardPar.props = {"chi", "chi_propagator"};
    forwardPar.defaultType = "PropagatorField";
    forwardPar.outProps = {"chi_t0.01", "chi_propagator_t0.01"};
    forwardPar.bc = -1;
    application.createModule<MGradientFlow::FermionFlow>("ForwardFlow", forwardPar);

    // Record the norms and inner products entering the discrete-adjoint test.
    MUtilities::NormCheckPropagator::Par backwardNormPar;
    backwardNormPar.field = "xi_t0.00";
    backwardNormPar.reference = "chi";
    application.createModule<MUtilities::NormCheckPropagator>("BackwardNorm", backwardNormPar);

    MUtilities::NormCheckPropagator::Par forwardNormPar;
    forwardNormPar.field = "eta";
    forwardNormPar.reference = "chi_t0.01";
    application.createModule<MUtilities::NormCheckPropagator>("ForwardNorm", forwardNormPar);

    MUtilities::NormCheckPropagator::Par stageNormPar;
    stageNormPar.field = "xi_t0.00";
    stageNormPar.reference = "xi_reconstructed_t0.00";
    application.createModule<MUtilities::NormCheckPropagator>("StageNorm", stageNormPar);

    // Condensate from the adjoint-flow path: D xi_propagator = xi at flow
    // time zero. With Identity and c_fl=0 this is -xi^dagger D^-1 xi.
    MContraction::StochasticCondensatePropagator::Par backwardCondensatePar;
    backwardCondensatePar.eta = "xi_t0.00";
    backwardCondensatePar.phi = "xi_propagator_t0.00";
    backwardCondensatePar.gamma = Gamma::Algebra::Identity;
    backwardCondensatePar.c_fl = 0.0;
    application.createModule<MContraction::StochasticCondensatePropagator>(
        "BackwardCondensate", backwardCondensatePar);

    // Positive-flow condensate: flow the stochastic source and its solution
    // together, then contract them at the positive flow time.
    MContraction::StochasticCondensatePropagator::Par forwardCondensatePar;
    forwardCondensatePar.eta = "chi_t0.01";
    forwardCondensatePar.phi = "chi_propagator_t0.01";
    forwardCondensatePar.gamma = Gamma::Algebra::Identity;
    forwardCondensatePar.c_fl = 0.0;
    application.createModule<MContraction::StochasticCondensatePropagator>(
        "ForwardCondensate", forwardCondensatePar);
    
    // Results to save
    std::vector<std::string> results = {
        "ForwardFlow", "BackwardNorm", "ForwardNorm", "StageNorm",
        "BackwardCondensate", "ForwardCondensate"
    };
    
    // save data
    MIO::WriteResultGroup::Par wPar;
    wPar.results = results;
    wPar.output = "results";
    application.createModule<MIO::WriteResultGroup>("WriteToFile", wPar);
    
    // execution
    application.saveParameterFile("TestAdjointFlow.xml");
    application.run();
    
    // epilogue
    LOG(Message) << "Grid is finalizing now" << std::endl;
    Grid_finalize();
    
    return EXIT_SUCCESS;
}
