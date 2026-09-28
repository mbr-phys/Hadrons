#include <Hadrons/Application.hpp>
#include <Hadrons/Modules.hpp>

using namespace Grid;
using namespace Hadrons;

int main(int argc, char *argv[])
{
    // initialization //////////////////////////////////////////////////////////
    Grid_init(&argc, &argv);
    HadronsLogError.Active(GridLogError.isActive());
    HadronsLogWarning.Active(GridLogWarning.isActive());
    HadronsLogMessage.Active(GridLogMessage.isActive());
    HadronsLogIterative.Active(GridLogIterative.isActive());
    HadronsLogDebug.Active(GridLogDebug.isActive());
    LOG(Message) << "Grid initialized" << std::endl;
    
    // run setup ///////////////////////////////////////////////////////////////
    Application              application;
    
    // global parameters
    Application::GlobalPar globalPar;
    globalPar.trajCounter.start    = 1500;
    globalPar.trajCounter.end      = 1520;
    globalPar.trajCounter.step     = 20;
    globalPar.runId                = "test_stochastic_meson";
    globalPar.genetic.maxGen       = 1000;
    globalPar.genetic.maxCstGen    = 200;
    globalPar.genetic.popSize      = 20;
    globalPar.genetic.mutationRate = 0.1;
    application.setPar(globalPar);
    
    // gauge field
    application.createModule<MGauge::Random>("gauge");

    // set fermion boundary conditions to be periodic space, antiperiodic time.
    std::string boundary = "1 1 1 -1";
    std::string twist = "0. 0. 0. 0.";

    // point source for standard propagator (for comparison)
    MSource::Point::Par ptPar;
    ptPar.position = "0 0 0 0";
    std::string srcName  = "pt_0";
    application.createModule<MSource::Point>(srcName, ptPar);

    // point sink
    MSink::Point::Par sinkPar;
    sinkPar.mom = "0 0 0";
    application.createModule<MSink::ScalarPoint>("sink", sinkPar);
    
    // action
    MAction::DWF::Par actionPar;
    actionPar.gauge = "gauge";
    actionPar.Ls    = 12;
    actionPar.M5    = 1.8;
    actionPar.mass  = 0.25;
    actionPar.boundary = boundary;
    actionPar.twist = twist;
    application.createModule<MAction::DWF>("DWF", actionPar);
    
    // solver
    MSolver::RBPrecCG::Par solverPar;
    solverPar.action       = "DWF";
    solverPar.residual     = 1.0e-8;
    solverPar.maxIteration = 10000;
    application.createModule<MSolver::RBPrecCG>("CG",solverPar);

    // standard propagator (for comparison)
    MFermion::GaugeProp::Par quarkPar;
    quarkPar.solver = "CG";
    quarkPar.source = srcName;
    std::vector<std::string> qName = {"Qpt_0"};
    application.createModule<MFermion::GaugeProp>(qName[0], quarkPar);

    // zero-flow meson from standard propagator (for comparison)
    MContraction::Meson::Par mesPar;
    mesPar.q1     = qName[0];
    mesPar.q2     = qName[0];
    mesPar.gammas = "(Gamma5 Gamma5)";
    mesPar.sink   = "sink";
    application.createModule<MContraction::Meson>("meson_std",mesPar);

    // Create two independent Z2 wall sources (A and B) and their propagators.
    // Contract them with StochasticMeson module.
    
    unsigned int sourceTime = 0;  // wall source time slice
    
    // Generate Z2 wall noise for source A at t = sourceTime
    MSource::Z2::Par z2ParA;
    z2ParA.tA = sourceTime;
    z2ParA.tB = sourceTime;
    std::string etaA = "eta_A";
    application.createModule<MSource::Z2>(etaA, z2ParA);
    
    // Generate Z2 wall noise for source B at t = sourceTime 
    MSource::Z2::Par z2ParB;
    z2ParB.tA = sourceTime;
    z2ParB.tB = sourceTime;
    std::string etaB = "eta_B";
    application.createModule<MSource::Z2>(etaB, z2ParB);
    
    // Solve for phi_A = D^-1 eta_A
    MFermion::GaugeProp::Par phiParA;
    phiParA.solver = "CG";
    phiParA.source = etaA;
    std::string phiA = "phi_A";
    application.createModule<MFermion::GaugeProp>(phiA, phiParA);
    
    // Solve for phi_B = D^-1 eta_B
    MFermion::GaugeProp::Par phiParB;
    phiParB.solver = "CG";
    phiParB.source = etaB;
    std::string phiB = "phi_B";
    application.createModule<MFermion::GaugeProp>(phiB, phiParB);
    
    // Stochastic meson contraction using StochasticMeson module
    // Use noise pair (A=0, B=0) at source time 0, zero flow time
    MContraction::StochasticMeson::Par stochMesPar;
    stochMesPar.phi1A      = phiA;
    stochMesPar.etaA       = etaA;
    stochMesPar.phi2B      = phiB;
    stochMesPar.etaB       = etaB;
    stochMesPar.gammas     = "(Gamma5 Gamma5)";  
    stochMesPar.sink       = "sink";
    stochMesPar.sourceTime = sourceTime;
    stochMesPar.noiseA     = 0;
    stochMesPar.noiseB     = 0;
    stochMesPar.flowTime   = 0.0;
    stochMesPar.output     = "stoch_meson_ps";
    application.createModule<MContraction::StochasticMeson>("stoch_meson_t0", stochMesPar);

    // Collect names of results to be written to file
    std::vector<std::string> results = {
        "meson_std",
        "stoch_meson_t0"
    };
    
    // save data
    MIO::WriteResultGroup::Par wPar;
    wPar.results = results;
    wPar.output = "results";
    application.createModule<MIO::WriteResultGroup>("WriteToFile", wPar);

    // execution
    application.saveParameterFile("StochasticMeson.xml");
    application.run();
    
    // epilogue
    LOG(Message) << "Grid is finalizing now" << std::endl;
    Grid_finalize();
    
    return EXIT_SUCCESS;
}
