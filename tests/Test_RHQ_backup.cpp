#include <Hadrons/Application.hpp>
#include <Hadrons/Modules.hpp>

using namespace Grid;
using namespace Hadrons;

int main(int argc, char *argv[])
{
    // *****************************
    // ****** initialization *******
    // *****************************

    Grid_init(&argc, &argv);
    HadronsLogError.Active(GridLogError.isActive());
    HadronsLogWarning.Active(GridLogWarning.isActive());
    HadronsLogMessage.Active(GridLogMessage.isActive());
    HadronsLogIterative.Active(GridLogIterative.isActive());
    HadronsLogDebug.Active(GridLogDebug.isActive());
    LOG(Message) << "Grid initialized" << std::endl;
    
    // run setup 
    Application              application;
    
    // global parameters
    Application::GlobalPar globalPar;
    globalPar.trajCounter.start             = 10000;
    globalPar.trajCounter.end               = 10020;
    globalPar.trajCounter.step              = 20;
    globalPar.runId                         = "test";
    application.setPar(globalPar);
    // gauge field
    MIO::LoadNersc::Par gauge;
    gauge.file = "../ckpoint_lat.IEEE64BIG";
    application.createModule<MIO::LoadNersc>("gauge",gauge);
    // sources
    MSource::Point::Par ptPar;
    ptPar.position = "0 0 0 0";
    application.createModule<MSource::Point>("pt", ptPar);
    // sink
    MSink::Point::Par sinkPar;
    sinkPar.mom = "0 0 0";
    application.createModule<MSink::ScalarPoint>("sink", sinkPar);
    
    // set fermion boundary conditions to be periodic space, antiperiodic time.
    std::string boundary = "1 1 1 -1";
    std::string twist = "0. 0. 0. 0.";

    // *****************************
    // *** strange quark modules ***
    // *****************************
    
    // action
    // MAction::DWF::Par actionPar;
    // actionPar.gauge = "gauge";
    // actionPar.Ls    = 16;
    // actionPar.M5    = 1.8;
    // actionPar.mass  = 0.32;
    // actionPar.boundary = boundary;
    // actionPar.twist = twist;
    // application.createModule<MAction::DWF>("DWF_s", actionPar);
    
    // solver
    // MSolver::RBPrecCG::Par solverPar;
    // solverPar.action       = "DWF_s";
    // solverPar.residual     = 1.0e-18;
    // solverPar.maxIteration = 8000;
    // application.createModule<MSolver::RBPrecCG>("CG_s",solverPar);
    
    // generate propagator
    // MFermion::GaugeProp::Par quarkPar;
    // quarkPar.solver = "CG_s";
    // quarkPar.source = "pt";
    // application.createModule<MFermion::GaugeProp>("Strange_Quark", quarkPar);

    // save propagator
    // MIO::SavePropagator::Par SaveProp;
    // SaveProp.name = "Strange_Quark";
    // SaveProp.fileStem = "../Strange_Quark_Propagator";
    // application.createModule<MIO::SavePropagator>("Save_Strange_Quark",SaveProp); 

    // load propagator
    MIO::LoadPropagator::Par LoadProp;
    LoadProp.name = "Strange_Quark";
    LoadProp.Ls = 1;
    LoadProp.fileStem = "../Strange_Quark_Propagator";
    application.createModule<MIO::LoadPropagator>("Load_Strange_Quark",LoadProp);

    // *****************************
    // ****** b quark modules ******
    // *****************************
    
    // RHQ action for b quark
    MAction::WilsonClover::Par Clover_b;
    Clover_b.gauge = "gauge";
    Clover_b.mass  = 7.47;
    Clover_b.csw_r = 4.92;
    Clover_b.csw_t = 4.92;
    WilsonAnisotropyCoefficients C_anisotropy;
    C_anisotropy.isAnisotropic = true;
    C_anisotropy.t_direction   = 3;
    C_anisotropy.xi_0          = 1;
    C_anisotropy.nu            = 2.93;
    Clover_b.clover_anisotropy = C_anisotropy;
    Clover_b.boundary = boundary;
    Clover_b.twist = twist;
    application.createModule<MAction::WilsonClover>("Clover_b", Clover_b);

    // Solver for RHQ Action
    MSolver::RBPrecCG::Par CG_b;
    CG_b.action       = "Clover_b";
    CG_b.residual     = 1e-45;
    CG_b.maxIteration = 1000;
    CG_b.eigenPack    = "";
    application.createModule<MSolver::RBPrecCG>("CG_b", CG_b);

    // b quark propagator
    MFermion::GaugeProp::Par B_Quark;
    B_Quark.solver = "CG_b";
    B_Quark.source = "pt";
    application.createModule<MFermion::GaugeProp>("B_Quark", B_Quark);

    // *****************************
    // ** begin sequential source **
    // *****************************

    MSource::SeqGamma::Par SeqStrangeSource;
    SeqStrangeSource.q = "Strange_Quark";
    SeqStrangeSource.tA = 14;
    SeqStrangeSource.tB = 14;
    SeqStrangeSource.gamma = Gamma::Algebra::Gamma5;
    SeqStrangeSource.mom = "0. 0. 0. 0.";
    application.createModule<MSource::SeqGamma>("SeqStrangeSource", SeqStrangeSource);

    // seq source with a strange quark
    MFermion::GaugeProp::Par SeqStrangeB;
    SeqStrangeB.solver = "CG_b";
    SeqStrangeB.source = "SeqStrangeSource";
    application.createModule<MFermion::GaugeProp>("SeqStrange-B", SeqStrangeB);

    // *****************************
    // ******* RHQ insertion *******
    // *****************************
    
    // RHQ Improvement II
    MRHQ::RHQInsertionI::Par RHQImprII_T;
    RHQImprII_T.q = "SeqStrange-B";
    RHQImprII_T.gauge = "gauge";
    RHQImprII_T.gauge_index = 3;
    application.createModule<MRHQ::RHQInsertionI>("RHQImprII_T", RHQImprII_T);

    // RHQ Improvement I
    MRHQ::RHQInsertionI::Par RHQImprI_T;
    RHQImprI_T.q = "Strange_Quark";
    RHQImprI_T.gauge = "gauge";
    RHQImprI_T.gauge_index = 3;
    application.createModule<MRHQ::RHQInsertionI>("RHQImprI_T", RHQImprI_T);

    // RHQ Improvement IV
    MRHQ::RHQInsertionIII::Par RHQImprIV;
    RHQImprIV.q = "SeqStrange-B";
    RHQImprIV.gauge = "gauge";
    application.createModule<MRHQ::RHQInsertionIII>("RHQImprIV", RHQImprIV);

    // RHQ Improvement III
    MRHQ::RHQInsertionIII::Par RHQImprIII;
    RHQImprIII.q = "Strange_Quark";
    RHQImprIII.gauge = "gauge";
    application.createModule<MRHQ::RHQInsertionIII>("RHQImprIII", RHQImprIII);

    // RHQ Improvement VI
    MRHQ::RHQInsertionV::Par RHQImprVI;
    RHQImprVI.q = "SeqStrange-B";
    RHQImprVI.gauge = "gauge";
    application.createModule<MRHQ::RHQInsertionV>("RHQImprVI", RHQImprVI);

    // RHQ Improvement V
    MRHQ::RHQInsertionV::Par RHQImprV;
    RHQImprV.q = "Strange_Quark";
    RHQImprV.gauge = "gauge";
    application.createModule<MRHQ::RHQInsertionV>("RHQImprV", RHQImprV);

    // *****************************
    // ***** Meson Contraction *****
    // *****************************

    // Light-Light Contraction (Eta_s 2-pt function)
    MContraction::Meson::Par pt_ll;
    pt_ll.q1 = "Strange_Quark";
    pt_ll.q2 = "Strange_Quark";
    pt_ll.gammas = "(Gamma5 Gamma5)";
    pt_ll.sink = "sink";
    pt_ll.output = "./eta_s_2pt";
    application.createModule<MContraction::Meson>("meson_pt_ll", pt_ll);

    // Light-Heavy Contraction (Bs meson 2-pt function)
    MContraction::Meson::Par pt_hl;
    pt_hl.q1 = "Strange_Quark";
    pt_hl.q2 = "B_Quark";
    pt_hl.gammas = "(Gamma5 Gamma5)";
    pt_hl.sink = "sink";
    pt_hl.output = "./b_meson_2pt";
    application.createModule<MContraction::Meson>("meson_pt_hl", pt_hl);

    // Heavy-Heavy Contraction (Bs meson 2-pt function)
    MContraction::Meson::Par pt_hh;
    pt_hh.q1 = "B_Quark";
    pt_hh.q2 = "B_Quark";
    pt_hh.gammas = "(Gamma5 Gamma5)";
    pt_hh.sink = "sink";
    pt_hh.output = "./bbar_2pt";
    application.createModule<MContraction::Meson>("meson_pt_hh", pt_hh);

    // *****************************
    // ****** 3pt Contraction ******
    // *****************************
    
    // Sequential Prop & Prop Contraction (Bs->Eta_s 3-pt function)
    MContraction::Meson::Par BstoEtas_3pt;
    BstoEtas_3pt.q1 = "Strange_Quark";
    BstoEtas_3pt.q2 = "B_Quark";
    BstoEtas_3pt.gammas = "(GammaT Gamma5)";
    BstoEtas_3pt.sink = "sink";
    BstoEtas_3pt.output = "./3pt_BstoEta_s_T_000";
    application.createModule<MContraction::Meson>("3pt_BstoEta_s_T", BstoEtas_3pt);
    
    // *****************************
    // ****** RHQ Contraction ******
    // *****************************
    
    // Bs->Eta_s 3-pt function, RHQII
    MContraction::Meson::Par BstoEtas_RHQII_3pt;
    BstoEtas_RHQII_3pt.q1 = "Strange_Quark";
    BstoEtas_RHQII_3pt.q2 = "RHQImprII_T";
    BstoEtas_RHQII_3pt.gammas = "(Identity Gamma5)";
    BstoEtas_RHQII_3pt.sink = "sink";
    BstoEtas_RHQII_3pt.output = "./3pt_BstoEta_s_T_RHQII_000";
    application.createModule<MContraction::Meson>("3pt_BstoEta_s_T_RHQII", BstoEtas_RHQII_3pt);

    // Bs->Eta_s 3-pt function, RHQI
    MContraction::Meson::Par BstoEtas_RHQI_3pt;
    BstoEtas_RHQI_3pt.q1 = "RHQImprI_T";
    BstoEtas_RHQI_3pt.q2 = "SeqStrange-B";
    BstoEtas_RHQI_3pt.gammas = "(Identity Gamma5)";
    BstoEtas_RHQI_3pt.sink = "sink";
    BstoEtas_RHQI_3pt.output = "./3pt_BstoEta_s_T_RHQI_000";
    application.createModule<MContraction::Meson>("3pt_BstoEta_s_T_RHQI", BstoEtas_RHQI_3pt);

    // Bs->Eta_s 3-pt function, RHQIV
    MContraction::Meson::Par BstoEtas_RHQIV_3pt;
    BstoEtas_RHQIV_3pt.q1 = "Strange_Quark";
    BstoEtas_RHQIV_3pt.q2 = "RHQImprIV";
    BstoEtas_RHQIV_3pt.gammas = "(GammaT Gamma5)";
    BstoEtas_RHQIV_3pt.sink = "sink";
    BstoEtas_RHQIV_3pt.output = "./3pt_BstoEta_s_T_RHQIV_000";
    application.createModule<MContraction::Meson>("3pt_BstoEta_s_T_RHQIV", BstoEtas_RHQIV_3pt);

    // Bs->Eta_s 3-pt function, RHQIII
    MContraction::Meson::Par BstoEtas_RHQIII_3pt;
    BstoEtas_RHQIII_3pt.q1 = "RHQImprIII";
    BstoEtas_RHQIII_3pt.q2 = "SeqStrange-B";
    BstoEtas_RHQIII_3pt.gammas = "(GammaT Gamma5)";
    BstoEtas_RHQIII_3pt.sink = "sink";
    BstoEtas_RHQIII_3pt.output = "./3pt_BstoEta_s_T_RHQIII_000";
    application.createModule<MContraction::Meson>("3pt_BstoEta_s_T_RHQIII", BstoEtas_RHQIII_3pt);

    // Bs->Eta_s 3-pt function, RHQVI
    MContraction::Meson::Par BstoEtas_RHQVI_3pt;
    BstoEtas_RHQVI_3pt.q1 = "Strange_Quark";
    BstoEtas_RHQVI_3pt.q2 = "RHQImprVI";
    BstoEtas_RHQVI_3pt.gammas = "(GammaT Gamma5)";
    BstoEtas_RHQVI_3pt.sink = "sink";
    BstoEtas_RHQVI_3pt.output = "./3pt_BstoEta_s_T_RHQVI_000";
    application.createModule<MContraction::Meson>("3pt_BstoEta_s_T_RHQVI", BstoEtas_RHQVI_3pt);

    // Bs->Eta_s 3-pt function, RHQV
    MContraction::Meson::Par BstoEtas_RHQV_3pt;
    BstoEtas_RHQV_3pt.q1 = "RHQImprV";
    BstoEtas_RHQV_3pt.q2 = "SeqStrange-B";
    BstoEtas_RHQV_3pt.gammas = "(GammaT Gamma5)";
    BstoEtas_RHQV_3pt.sink = "sink";
    BstoEtas_RHQV_3pt.output = "./3pt_BstoEta_s_T_RHQV_000";
    application.createModule<MContraction::Meson>("3pt_BstoEta_s_T_RHQV", BstoEtas_RHQV_3pt);

    // execution
    application.run();
    
    // epilogue
    LOG(Message) << "Grid is finalizing now" << std::endl;
    Grid_finalize();
    
    return EXIT_SUCCESS;
}
