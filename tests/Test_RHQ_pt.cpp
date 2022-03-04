#include <Hadrons/Application.hpp>
#include <Hadrons/Modules.hpp>

using namespace Grid;
using namespace Hadrons;

// set fermion boundary conditions to be periodic space, antiperiodic time.
// these will likely need to be generalised, e.g. in xml, at some point
std::string boundary = "1 1 1 -1";
std::string twist = "0. 0. 0. 0.";

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
                                        int, maxiter,
                                        std::string, mom);
    };

    class MesonPar : Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(MesonPar,
                                        std::string, cont_name,
                                        std::string, gamma_src,
                                        std::string, gamma_snk,
                                        std::string, snk_mom);
    };
}

struct TestPar
{
    TestInputs::RunPar      runPar;
    TestInputs::ConfigPar   configPar;
    TestInputs::DwfPar      dwfPar1;
    TestInputs::DwfPar      dwfPar2;
    TestInputs::MesonPar    mesonPar;
};

std::string removeSpaces(std::string word) {
    std::string newWord;
    for (int i = 0; i < word.length(); i++) {
        if (word[i] != ' ') {
            newWord += word[i];}
    }
    return newWord;
}

// dwf propagator --> add flag for loading from / saving to binary?
std::string make_dwfpropagator(Application &application, TestInputs::DwfPar &dwfPar, std::string src){
    // action
    MAction::DWF::Par actionPar;
    actionPar.gauge    = "gauge"; // maybe generalise this somehow, if at all necessary
    actionPar.Ls       = dwfPar.Ls;
    actionPar.M5       = dwfPar.M5;
    actionPar.mass     = dwfPar.mass;
    actionPar.boundary = boundary;
    actionPar.twist    = twist;

    std::string act_name = "DWF_action_" + dwfPar.quark;
    application.createModule<MAction::DWF>(act_name, actionPar);

    // solver
    MSolver::RBPrecCG::Par solverPar;
    solverPar.action       = act_name;
    solverPar.residual     = dwfPar.residual;
    solverPar.maxIteration = dwfPar.maxiter;

    std::string solve_name = "CG_" + dwfPar.quark;
    application.createModule<MSolver::RBPrecCG>(solve_name,solverPar);

    // generate propagator
    MFermion::GaugeProp::Par quarkPar;
    quarkPar.solver = solve_name;
    quarkPar.source = src;

    std::string prop_name = dwfPar.quark + "_quark";
    application.createModule<MFermion::GaugeProp>(prop_name, quarkPar);
    return prop_name;
}

// (phased) pt wall source
std::string make_ptsrc(Application &application, TestInputs::DwfPar &dwfPar){
    // start with standard pt source
    MSource::Point::Par ptPar;
    ptPar.position = "0 0 0 0"; // maybe put this as xml input too
    std::string pt_mom = "pt_" + dwfPar.quark + "_p000_t0";
    application.createModule<MSource::Point>(pt_mom, ptPar);
    // then add momentum if not 0
    if (dwfPar.mom != "0 0 0") {
        MSource::MomentumPhase::Par momPar;
        momPar.src = pt_mom;
        momPar.mom = dwfPar.mom + " 0";
        std::string mom_str = dwfPar.mom;
        pt_mom = "pt_" + dwfPar.quark + "_p" + removeSpaces(mom_str) + "_t0";
        application.createModule<MSource::MomentumPhase>(pt_mom,momPar);
    }
    return pt_mom;
}

// general contraction --> from Alessandro
std::string make_contraction(Application &application, std::string contraction_name, std::string q1, std::string q2, std::array<std::string, 2> gammas, std::string sink){

    MContraction::Meson::Par contraction;
    std::string gamma_snk = gammas[0];
    std::string gamma_src = gammas[1];

    contraction.q1 = q1;
    contraction.q2 = q2;
    if ((gamma_snk == "all") || (gamma_src == "all")) {
        contraction.gammas = "all";
        contraction_name += "_all";
    }
    else {
        contraction.gammas = "(" + gamma_snk + " " + gamma_src + ")";
        contraction_name += "_" + gamma_snk + "_" + gamma_src;
    }
    contraction.sink   = sink;

    contraction.output = contraction_name;

    application.createModule<MContraction::Meson>(contraction_name, contraction);
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
    read(reader, "dwfPar1",    testPar.dwfPar1);
    read(reader, "dwfPar2",    testPar.dwfPar2);
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
    Application application;
    
    // global parameters
    Application::GlobalPar globalPar;
    globalPar.trajCounter.start    = testPar.configPar.begin;
    globalPar.trajCounter.end      = testPar.configPar.end;
    globalPar.trajCounter.step     = testPar.configPar.step;
    globalPar.runId                = testPar.runPar.runId;
    application.setPar(globalPar);
    
    // gauge field
    MIO::LoadNersc::Par gauge;
    gauge.file = testPar.configPar.fileStem;
    application.createModule<MIO::LoadNersc>("gauge",gauge);

    // sources
    std::string src1 = make_ptsrc(application, testPar.dwfPar1);

    std::string src2 = make_ptsrc(application, testPar.dwfPar2);

    // sink 
    MSink::Point::Par sinkPar;
    sinkPar.mom = testPar.mesonPar.snk_mom;
    application.createModule<MSink::ScalarPoint>("sink", sinkPar);
    
    // *****************************
    // ******** dwf modules ********
    // *****************************
    
    std::string name_q1 = make_dwfpropagator(application, testPar.dwfPar1, src1); 

    std::string name_q2 = make_dwfpropagator(application, testPar.dwfPar2, src2);

    // *****************************
    // ***** Meson Contraction *****
    // *****************************

    // generalise to loop over multiple gamma structure but not all of them, if this could be useful
    // --> if either gamma_snk, gamma_src == 'all', will give all 256 entries 
    std::vector<std::string> gamma_src_vec = split_gammas(testPar.mesonPar.gamma_src,' ');
    std::vector<std::string> gamma_snk_vec = split_gammas(testPar.mesonPar.gamma_snk,' ');
    // program expects same number of snk and src gammas
    if (gamma_src_vec.size() != gamma_snk_vec.size()) {
        std::cerr << "gamma_src and gamma_snk must have same number of entries:" << std::endl
                  << "    gamma_src has " << gamma_src_vec.size() << " entries," << std::endl
                  << "    gamma_snk has " << gamma_snk_vec.size() << " entries." << std::endl;

        return EXIT_FAILURE;
    }
    for (int i = 0; i < (int)gamma_src_vec.size(); i++) {
        std::array<std::string, 2> Gammas = {gamma_snk_vec[i], gamma_src_vec[i]};
        make_contraction(application, testPar.mesonPar.cont_name, name_q1, name_q2, Gammas, "sink");
    }

    // execution
    application.run();
    
    // epilogue
    LOG(Message) << "Grid is finalizing now" << std::endl;
    Grid_finalize();
    
    return EXIT_SUCCESS;
}
