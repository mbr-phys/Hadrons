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
    double        mass    = .25;
    
    // global parameters
    Application::GlobalPar globalPar;
    globalPar.trajCounter.start    = 1500;
    globalPar.trajCounter.end      = 1520;
    globalPar.trajCounter.step     = 20;
    globalPar.runId                = "test";
    globalPar.genetic.maxGen       = 1000;
    globalPar.genetic.maxCstGen    = 200;
    globalPar.genetic.popSize      = 20;
    globalPar.genetic.mutationRate = .1;
    application.setPar(globalPar);
    
    // gauge field
    application.createModule<MGauge::Random>("gauge");

    // set fermion boundary conditions to be periodic space, antiperiodic time.
    std::string boundary = "1 1 1 -1";
    std::string twist = "0. 0. 0. 0.";

    // point source
    MSource::Point::Par ptPar;
    ptPar.position = "0 0 0 0";
    std:: string srcName  = "pt_0";
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
    actionPar.mass  = mass;
    actionPar.boundary = boundary;
    actionPar.twist = twist;
    application.createModule<MAction::DWF>("DWF", actionPar);
    
    // solver
    MSolver::RBPrecCG::Par solverPar;
    solverPar.action       = "DWF";
    solverPar.residual     = 1.0e-8;
    solverPar.maxIteration = 10000;
    application.createModule<MSolver::RBPrecCG>("CG",solverPar);

    // propagator
    MFermion::GaugeProp::Par quarkPar;
    quarkPar.solver = "CG";
    quarkPar.source = srcName;
    std::vector<std::string> qName = {"Qpt_0"};
    application.createModule<MFermion::GaugeProp>(qName[0], quarkPar);

    // zero flow contractions
    MContraction::Meson::Par mesPar;
    mesPar.q1     = qName[0];
    mesPar.q2     = qName[0];
    mesPar.gammas = "all";
    mesPar.sink   = "sink";
    application.createModule<MContraction::Meson>("meson_t0.00",mesPar);

    std::vector<std::string> results = {"meson_t0.00"}; // collect names of results to be written to file
    
    // ///////////////////////////////////////////////////////////////////////
    // POSITIVE FLOW: Flow both noise and solution 
    // together in a single module execution. This shares gauge evolution
    // across both field types, avoiding redundant gauge RK stage computation.
    //
    // For positive-flow bilinear estimators:
    //   - eta: stochastic noise field 
    //   - phi: solution D^-1 eta 
    //   - Both are flowed forward with the same gauge trajectory
    //   - Contraction: -eta_t^dagger A_t phi_t (scalar: A_t = 1)
    // ///////////////////////////////////////////////////////////////////////
    
    // Generate sparse Z2 stochastic noise 
    // creates spin-color diagonal noise with nsparse dilution
    unsigned int nsrc = 1;      // number of noise sources
    unsigned int nsparse = 2;   // sparse dilution factor (2^4 = 16 dilutions)
    std::string noiseBase = "eta";
    
    // create sparse spin-color diagonal noise
    MNoise::SparseSpinColorDiagonal::Par sparsePar;
    sparsePar.nsrc = nsrc;
    sparsePar.nsparse = nsparse;
    std::string sparseNoiseName = noiseBase + "_sparse";
    application.createModule<MNoise::SparseSpinColorDiagonal>(sparseNoiseName, sparsePar);
    
    // apply Z2 dilution 
    MSource::Z2Diluted::Par dilutedPar;
    dilutedPar.noise = sparseNoiseName;
    application.createModule<MSource::Z2Diluted>(noiseBase, dilutedPar);
    
    // unpack the diluted noise into individual sources
    // The output is a vector of PropagatorField: eta_0_0, eta_0_1, ..., eta_0_(N-1)
    unsigned int nDilutions = pow(nsparse, 4);  // 4D dilution
    std::vector<std::string> etaNames(nDilutions);
    MUtilities::PropagatorVectorUnpack::Par unpackPar;
    unpackPar.input = noiseBase;
    unpackPar.size = nDilutions;
    std::string unpackName = noiseBase + "_unpacked";
    application.createModule<MUtilities::PropagatorVectorUnpack>(unpackName, unpackPar);
    
    // Get the unpacked noise names
    for (unsigned int i = 0; i < nDilutions; i++) {
        etaNames[i] = unpackName + "_" + std::to_string(i);
    }
    
    // For this test, just use the first noise source
    std::string etaName = etaNames[0];
    
    // Solve for phi = D^-1 eta
    MFermion::GaugeProp::Par phiPar;
    phiPar.solver = "CG";
    phiPar.source = etaName;
    std::vector<std::string> phiName = {"phi_0"};
    application.createModule<MFermion::GaugeProp>(phiName[0], phiPar);
    
    // Positive flow: flow eta, phi, AND the standard propagator Qpt_0 together
    // This is the key efficiency win - gauge evolution happens ONCE for all fields
    MGradientFlow::FermionFlow::Par positiveFlowPar;
    positiveFlowPar.gauge = "gauge";
    positiveFlowPar.steps = 10;
    positiveFlowPar.step_size = 0.01;
    positiveFlowPar.meas_interval = 10;
    positiveFlowPar.props = {etaName, phiName[0], qName[0]};  // three fields together
    positiveFlowPar.defaultType = "PropagatorField";          // homogeneous list
    //positiveFlowPar.propTypes = {"PropagatorField", "PropagatorField", "PropagatorField"}; // alternatively define type of each object, e.g. if using FermionFields as well
    positiveFlowPar.bc = -1;
    application.createModule<MGradientFlow::FermionFlow>("PositiveFlow", positiveFlowPar);
    
    // Positive-flow contractions using new StochasticCondensate module
    // Scalar condensate at flow time t=0.10
    MContraction::StochasticCondensatePropagator::Par scalarPar;
    scalarPar.eta = etaName + "_t0.10";
    scalarPar.phi = phiName[0] + "_t0.10";
    scalarPar.gamma = "Identity";
    scalarPar.c_fl = 0.0;  // DWF action, chiral symmetry protects
    application.createModule<MContraction::StochasticCondensatePropagator>("scalar_t0.10", scalarPar);
    results.push_back("scalar_t0.10");
    
    // Derivative condensate for Z_chi (ringed scheme) using DslashField + StochasticCondensate
    // Step 1: Apply D-slash to flowed phi
    MContraction::DslashFieldPropagator::Par dslashPar;
    dslashPar.input = phiName[0] + "_t0.10";
    dslashPar.gauge = "PositiveFlow_U_t0.10";  // flowed gauge at same flow time
    application.createModule<MContraction::DslashFieldPropagator>("Dslash_phi_t0.10", dslashPar);
    
    // Step 2: Contract eta with Dslash_phi
    MContraction::StochasticCondensatePropagator::Par derivPar;
    derivPar.eta = etaName + "_t0.10";
    derivPar.phi = "Dslash_phi_t0.10";
    derivPar.gamma = "Identity";
    derivPar.c_fl = 0.0;  // no c_fl for derivative condensate
    application.createModule<MContraction::StochasticCondensatePropagator>("deriv_condensate_t0.10", derivPar);
    results.push_back("deriv_condensate_t0.10");
    
    // Flowed standard propagator contractions (for comparison with stochastic)
    MContraction::Meson::Par mesflowPar;
    mesflowPar.q1     = qName[0]+"_t0.10";  // flowed standard propagator
    mesflowPar.q2     = qName[0]+"_t0.10";  // flowed standard propagator
    mesflowPar.gammas = "all";
    mesflowPar.sink   = "sink";
    application.createModule<MContraction::Meson>("meson_std_t0.10",mesflowPar);
    results.push_back("meson_std_t0.10");
    
    // save data
    MIO::WriteResultGroup::Par wPar;
    wPar.results = results;
    wPar.output = "results";
    application.createModule<MIO::WriteResultGroup>("WriteToFile",wPar);

    // execution
    application.saveParameterFile("FermionFlow.xml");
    application.run();
    
    // epilogue
    LOG(Message) << "Grid is finalizing now" << std::endl;
    Grid_finalize();
    
    return EXIT_SUCCESS;
}
