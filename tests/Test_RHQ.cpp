#include <Hadrons/Application.hpp>
#include <Hadrons/Modules.hpp>

using namespace Grid;
using namespace Hadrons;

namespace TestInputs
{
    class RunPar : Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(RunPar, 
                                        std::string, runId);
    };

    class ConfigPar : Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(ConfigPar,
                                        std::string,  fileStem,
                                        unsigned int, begin,
                                        unsigned int, end,
                                        unsigned int, step);
    };
    
    class DwfPar : Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(DwfPar,
                                        std::string, quark,
                                        double, mass,
                                        double, M5,
                                        unsigned int, Ls,
                                        double, residual,
                                        unsigned int, maxiter);
    };

    class RhqPar : Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(RhqPar,
                                        std::string, quark,
                                        double, mass,
                                        double, cswr,
                                        double, cswt,
                                        unsigned int, tdir,
                                        unsigned int, xi0,
                                        double, nu,
                                        double, residual,
                                        unsigned int, maxiter);
    };

    class MesonPar : Serializable 
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(MesonPar,
                                        std::string, cont_name,
                                        std::string, gamma_src,
                                        std::string, gamma_snk);
    };
}

struct TestPar
{
    TestInputs::RunPar      runPar;
    TestInputs::ConfigPar   configPar;
    TestInputs::DwfPar      dwfPar;
    TestInputs::RhqPar      rhqPar;
    TestInputs::MesonPar    mesonPar;
};

std::vector<std::string> split_gammas(std::string str, char del){
    std::string temp = "";
    std::vector<std::string> strings;
    for (int i = 0; i < (int)str.size(); i++) {
        if (str[i] != del) {
            temp += str[i];
        }
        else {
            strings.push_back(temp);
            temp = "";
        }
    }
    if (temp != "") {
        strings.push_back(temp);
    }
    return strings;
}

// general contraction 
std::string make_contraction(Application &application, std::string meson, std::string q1, std::string q2, std::array<std::string, 2> gammas, std::string sink, std::string fileSt){

    MContraction::Meson::Par contraction;
    std::string gamma_snk = gammas[0];
    std::string gamma_src = gammas[1];
    std::string contraction_name = "";

    contraction.q1 = q1;
    contraction.q2 = q2;
    if ((gamma_snk == "all") || (gamma_src == "all")) {
        contraction.gammas = "all";
        contraction_name = "all";
    }
    else {
        contraction.gammas = "(" + gamma_snk + " " + gamma_src + ")";
        contraction_name = gamma_snk + "_" + gamma_src; 
                                                       
    }
    contraction.sink   = sink;

    contraction.output = fileSt + meson + "_" + contraction_name;

    application.createModule<MContraction::Meson>(meson + "_" + contraction_name, contraction);
    return contraction_name;
}

int main(int argc, char *argv[])
{
    //program expects parameter file as first argument of command line
    if (argc < 2)
    {
        std::cerr << "usage: " << argv[0] << " <parameter file>";
        std::cerr << std::endl;

        return EXIT_FAILURE;
    }
    std::string parFilename;
    parFilename = argv[1];

    TestPar testPar;
    XmlReader reader(parFilename);

    read(reader, "runPar",     testPar.runPar);
    read(reader, "configPar",  testPar.configPar);
    read(reader, "dwfPar",     testPar.dwfPar);
    read(reader, "rhqPar",     testPar.rhqPar);
    read(reader, "mesonPar",   testPar.mesonPar);

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
    globalPar.trajCounter.start             = testPar.configPar.begin;
    globalPar.trajCounter.end               = testPar.configPar.end;
    globalPar.trajCounter.step              = testPar.configPar.step;
    globalPar.runId                         = testPar.runPar.runId;
    application.setPar(globalPar);
    // gauge field
    MIO::LoadNersc::Par gauge;
    gauge.file = testPar.configPar.fileStem;
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

    // make vectors of src and snk gamma structures to loop over if more than 1 set
    std::vector<std::string> gamma_src_vec = split_gammas(testPar.mesonPar.gamma_src,' ');                     
    std::vector<std::string> gamma_snk_vec = split_gammas(testPar.mesonPar.gamma_snk,' ');
    // program expects same number of snk and src gammas
    if (gamma_src_vec.size() != gamma_snk_vec.size()) {
        std::cerr << "gamma_src and gamma_snk must have same number of entries:" << std::endl
                  << "    gamma_src has " << gamma_src_vec.size() << " entries," << std::endl
                  << "    gamma_snk has " << gamma_snk_vec.size() << " entries." << std::endl;

        return EXIT_FAILURE;                                                                                   
    }       

    // *****************************
    // *** strange quark modules ***
    // *****************************
    
    // action
    // MAction::DWF::Par actionPar;
    // actionPar.gauge = "gauge";
    // actionPar.Ls    = testPar.dwfPar.Ls;
    // actionPar.M5    = testPar.dwfPar.M5;
    // actionPar.mass  = testPar.dwfPar.mass;
    // actionPar.boundary = boundary;
    // actionPar.twist = twist;
    // application.createModule<MAction::DWF>("DWF_s", actionPar);
    
    // solver
    // MSolver::RBPrecCG::Par solverPar;
    // solverPar.action       = "DWF_s";
    // solverPar.residual     = testPar.dwfPar.residual;
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
    Clover_b.mass  = testPar.rhqPar.mass;
    Clover_b.csw_r = testPar.rhqPar.cswr;
    Clover_b.csw_t = testPar.rhqPar.cswt;
    WilsonAnisotropyCoefficients C_anisotropy;
    C_anisotropy.isAnisotropic = true;
    C_anisotropy.t_direction   = testPar.rhqPar.tdir;
    C_anisotropy.xi_0          = testPar.rhqPar.xi0;
    C_anisotropy.nu            = testPar.rhqPar.nu;
    Clover_b.clover_anisotropy = C_anisotropy;
    Clover_b.boundary = boundary;
    Clover_b.twist = twist;
    application.createModule<MAction::WilsonClover>("Clover_b", Clover_b);

    // Solver for RHQ Action
    MSolver::RBPrecCG::Par CG_b;
    CG_b.action       = "Clover_b";
    CG_b.residual     = testPar.rhqPar.residual;
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

//    MSource::SeqGamma::Par SeqStrangeSource;
//    SeqStrangeSource.q = "Strange_Quark";
//    SeqStrangeSource.tA = 14;
//    SeqStrangeSource.tB = 14;
//    SeqStrangeSource.gamma = Gamma::Algebra::Gamma5;
//    SeqStrangeSource.mom = "0. 0. 0. 0.";
//    application.createModule<MSource::SeqGamma>("SeqStrangeSource", SeqStrangeSource);
//
//    // seq source with a strange quark
//    MFermion::GaugeProp::Par SeqStrangeB;
//    SeqStrangeB.solver = "CG_b";
//    SeqStrangeB.source = "SeqStrangeSource";
//    application.createModule<MFermion::GaugeProp>("SeqStrange-B", SeqStrangeB);

    // *****************************
    // ******* RHQ insertion *******
    // *****************************
    
    // RHQ Improvement II 
    MRHQ::RHQInsertionI::Par RHQImprII_T;
    RHQImprII_T.q = "B_Quark";
    RHQImprII_T.gauge = "gauge";
    RHQImprII_T.index = 3;
    RHQImprII_T.gamma5 = Gamma::Algebra::Gamma5;
    application.createModule<MRHQ::RHQInsertionI>("RHQImprII", RHQImprII_T);

    // RHQ Improvement I 
    MRHQ::RHQInsertionI::Par RHQImprI_T;
    RHQImprI_T.q = "Strange_Quark";
    RHQImprI_T.gauge = "gauge";
    RHQImprI_T.index = 3;
    RHQImprI_T.gamma5 = Gamma::Algebra::Gamma5;
    application.createModule<MRHQ::RHQInsertionI>("RHQImprI", RHQImprI_T);

    // RHQ Improvement IV 
    MRHQ::RHQInsertionIV::Par RHQImprIV;
    RHQImprIV.q = "B_Quark";
    RHQImprIV.gauge = "gauge";
    RHQImprIV.index = 3;
    RHQImprIV.gamma5 = Gamma::Algebra::Gamma5;
    RHQImprIV.flag = 0;
    application.createModule<MRHQ::RHQInsertionIV>("RHQImprIV", RHQImprIV);

    // RHQ Improvement III 
    MRHQ::RHQInsertionIII::Par RHQImprIII;
    RHQImprIII.q = "Strange_Quark";
    RHQImprIII.gauge = "gauge";
    RHQImprIII.index = 3;
    RHQImprIII.gamma5 = Gamma::Algebra::Gamma5;
    RHQImprIII.flag = 0;
    application.createModule<MRHQ::RHQInsertionIII>("RHQImprIII", RHQImprIII);

    // RHQ Improvement VI 
    MRHQ::RHQInsertionVI::Par RHQImprVI;
    RHQImprVI.q = "B_Quark";
    RHQImprVI.gauge = "gauge";
    RHQImprVI.index = 3;
    RHQImprVI.gamma5 = Gamma::Algebra::Gamma5;
    application.createModule<MRHQ::RHQInsertionVI>("RHQImprVI", RHQImprVI);

    // RHQ Improvement V 
    MRHQ::RHQInsertionV::Par RHQImprV;
    RHQImprV.q = "Strange_Quark";
    RHQImprV.gauge = "gauge";
    RHQImprV.index = 3;
    RHQImprV.gamma5 = Gamma::Algebra::Gamma5;
    application.createModule<MRHQ::RHQInsertionV>("RHQImprV", RHQImprV);

    // *****************************
    // ***** Meson Contraction *****
    // *****************************

    for (int i = 0; i < (int)gamma_src_vec.size(); i++) {
        std::array<std::string, 2> Gammas = {gamma_snk_vec[i], gamma_src_vec[i]};
//        make_contraction(application, "eta_s_2pt", "Strange_Quark", "Strange_Quark", Gammas, "sink", "./");
//        make_contraction(application, "bbar_2pt", "B_Quark", "B_Quark", Gammas, "sink", "./");
        make_contraction(application, "sb_O0", "Strange_Quark", "B_Quark", Gammas, "sink", "./");
    }
    std::array<std::string, 2> Gammas = {"Identity", "Gamma5"};
    make_contraction(application, "sb_O1", "RHQImprI", "B_Quark", Gammas, "sink", "./");
    make_contraction(application, "sb_O2", "Strange_Quark", "RHQImprII", Gammas, "sink", "./");
    make_contraction(application, "sb_O3", "RHQImprIII", "B_Quark", Gammas, "sink", "./");
    make_contraction(application, "sb_O4", "Strange_Quark", "RHQImprIV", Gammas, "sink", "./");
    make_contraction(application, "sb_O5", "RHQImprV", "B_Quark", Gammas, "sink", "./");
    make_contraction(application, "sb_O6", "Strange_Quark", "RHQImprVI", Gammas, "sink", "./");

    // Light-Light Contraction (Eta_s 2-pt function)
//    MContraction::Meson::Par pt_ll;
//    pt_ll.q1 = "Strange_Quark";
//    pt_ll.q2 = "Strange_Quark";
//    pt_ll.gammas = "(Gamma5 Gamma5)";
//    pt_ll.sink = "sink";
//    pt_ll.output = "./eta_s_2pt";
//    application.createModule<MContraction::Meson>("meson_pt_ll", pt_ll);

//    // Light-Heavy Contraction (Bs meson 2-pt function)
//    MContraction::Meson::Par pt_hl;
//    pt_hl.q1 = "Strange_Quark";
//    pt_hl.q2 = "B_Quark";
//    pt_hl.gammas = "(Gamma5 Gamma5)";
//    pt_hl.sink = "sink";
//    pt_hl.output = "./b_meson_2pt";
//    application.createModule<MContraction::Meson>("meson_pt_hl", pt_hl);

//    // Heavy-Heavy Contraction (Bs meson 2-pt function)
//    MContraction::Meson::Par pt_hh;
//    pt_hh.q1 = "B_Quark";
//    pt_hh.q2 = "B_Quark";
//    pt_hh.gammas = "(Gamma5 Gamma5)";
//    pt_hh.sink = "sink";
//    pt_hh.output = "./bbar_2pt";
//    application.createModule<MContraction::Meson>("meson_pt_hh", pt_hh);

    // *****************************
    // ****** 3pt Contraction ******
    // *****************************
    
    // Sequential Prop & Prop Contraction (Bs->Eta_s 3-pt function)
//    MContraction::Meson::Par BstoEtas_3pt;
//    BstoEtas_3pt.q1 = "Strange_Quark";
//    BstoEtas_3pt.q2 = "B_Quark";
//    BstoEtas_3pt.gammas = "(GammaT Gamma5)";
//    BstoEtas_3pt.sink = "sink";
//    BstoEtas_3pt.output = "./3pt_BstoEta_s_T_000";
//    application.createModule<MContraction::Meson>("3pt_BstoEta_s_T", BstoEtas_3pt);
    
    // *****************************
    // ****** RHQ Contraction ******
    // *****************************
    
    // Bs->Eta_s 3-pt function, RHQII
//    MContraction::Meson::Par BstoEtas_RHQII_3pt;
//    BstoEtas_RHQII_3pt.q1 = "Strange_Quark";
//    BstoEtas_RHQII_3pt.q2 = "RHQImprII_T";
//    BstoEtas_RHQII_3pt.gammas = "(Identity Gamma5)";
//    BstoEtas_RHQII_3pt.sink = "sink";
//    BstoEtas_RHQII_3pt.output = "./3pt_BstoEta_s_T_RHQII_000";
//    application.createModule<MContraction::Meson>("3pt_BstoEta_s_T_RHQII", BstoEtas_RHQII_3pt);
//
//    // Bs->Eta_s 3-pt function, RHQI
//    MContraction::Meson::Par BstoEtas_RHQI_3pt;
//    BstoEtas_RHQI_3pt.q1 = "RHQImprI_T";
//    BstoEtas_RHQI_3pt.q2 = "SeqStrange-B";
//    BstoEtas_RHQI_3pt.gammas = "(Identity Gamma5)";
//    BstoEtas_RHQI_3pt.sink = "sink";
//    BstoEtas_RHQI_3pt.output = "./3pt_BstoEta_s_T_RHQI_000";
//    application.createModule<MContraction::Meson>("3pt_BstoEta_s_T_RHQI", BstoEtas_RHQI_3pt);
//
//    // Bs->Eta_s 3-pt function, RHQIV
//    MContraction::Meson::Par BstoEtas_RHQIV_3pt;
//    BstoEtas_RHQIV_3pt.q1 = "Strange_Quark";
//    BstoEtas_RHQIV_3pt.q2 = "RHQImprIV";
//    BstoEtas_RHQIV_3pt.gammas = "(GammaT Gamma5)";
//    BstoEtas_RHQIV_3pt.sink = "sink";
//    BstoEtas_RHQIV_3pt.output = "./3pt_BstoEta_s_T_RHQIV_000";
//    application.createModule<MContraction::Meson>("3pt_BstoEta_s_T_RHQIV", BstoEtas_RHQIV_3pt);
//
//    // Bs->Eta_s 3-pt function, RHQIII
//    MContraction::Meson::Par BstoEtas_RHQIII_3pt;
//    BstoEtas_RHQIII_3pt.q1 = "RHQImprIII";
//    BstoEtas_RHQIII_3pt.q2 = "SeqStrange-B";
//    BstoEtas_RHQIII_3pt.gammas = "(GammaT Gamma5)";
//    BstoEtas_RHQIII_3pt.sink = "sink";
//    BstoEtas_RHQIII_3pt.output = "./3pt_BstoEta_s_T_RHQIII_000";
//    application.createModule<MContraction::Meson>("3pt_BstoEta_s_T_RHQIII", BstoEtas_RHQIII_3pt);
//
//    // Bs->Eta_s 3-pt function, RHQVI
//    MContraction::Meson::Par BstoEtas_RHQVI_3pt;
//    BstoEtas_RHQVI_3pt.q1 = "Strange_Quark";
//    BstoEtas_RHQVI_3pt.q2 = "RHQImprVI";
//    BstoEtas_RHQVI_3pt.gammas = "(GammaT Gamma5)";
//    BstoEtas_RHQVI_3pt.sink = "sink";
//    BstoEtas_RHQVI_3pt.output = "./3pt_BstoEta_s_T_RHQVI_000";
//    application.createModule<MContraction::Meson>("3pt_BstoEta_s_T_RHQVI", BstoEtas_RHQVI_3pt);
//
//    // Bs->Eta_s 3-pt function, RHQV
//    MContraction::Meson::Par BstoEtas_RHQV_3pt;
//    BstoEtas_RHQV_3pt.q1 = "RHQImprV";
//    BstoEtas_RHQV_3pt.q2 = "SeqStrange-B";
//    BstoEtas_RHQV_3pt.gammas = "(GammaT Gamma5)";
//    BstoEtas_RHQV_3pt.sink = "sink";
//    BstoEtas_RHQV_3pt.output = "./3pt_BstoEta_s_T_RHQV_000";
//    application.createModule<MContraction::Meson>("3pt_BstoEta_s_T_RHQV", BstoEtas_RHQV_3pt);

    // execution
    // application.saveParameterFile("spectrum.xml");
    application.run();
    
    // epilogue
    LOG(Message) << "Grid is finalizing now" << std::endl;
    Grid_finalize();
    
    return EXIT_SUCCESS;
}
