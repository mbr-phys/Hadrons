#include <Hadrons/Application.hpp>
#include <Hadrons/Modules.hpp>

using namespace Grid;
using namespace Hadrons;

// set fermion boundary conditions to be periodic space, antiperiodic time.
// these will likely need to be generalised, e.g. in xml, at some point
std::string boundary = "1 1 1 -1";
std::string twist = "0. 0. 0. 0.";

// momenta
std::vector<std::string> moms = {
    "0 0 0",                                                // n2 = 0
    "1 0 0","-1 0 0","0 1 0","0 -1 0","0 0 1","0 0 -1",     // n2 = 1
    "1 1 0","-1 -1 0","0 1 1","0 -1 -1","1 0 1","-1 0 -1",  // n2 = 2 
    "1 -1 0","-1 1 0","0 1 -1","0 -1 1","1 0 -1","-1 0 1",
    "1 1 1","-1 -1 -1","1 1 -1","-1 -1 1","1 -1 -1",        // n2 = 3
    "-1 1 1","1 -1 1","-1 1 -1",
    "2 0 0","-2 0 0","0 2 0","0 -2 0","0 0 2","0 0 -2"      // n2 = 4
};

int find_ind(std::string Mom) {
    int ind;
    for (int i = 0; i < (int)moms.size(); i++) {
        if (moms[i] == Mom) {
            ind = i;
            break;
        }
    }
    return ind;
}

std::string squareMom(std::string Mom) {
    int ind = find_ind(Mom);
    std::string sqM = "n2_";
    if (ind == 0) {
        sqM += "0";}
    else if (ind > 0 && ind <= 6) {
        sqM += "1";}
    else if (ind > 6 && ind <= 18) {
        sqM += "2";}
    else if (ind > 18 && ind <= 26) {
        sqM += "3";}
    else if (ind > 26) {
        sqM += "4";}
    return sqM;
}

int count_char(std::string str, char chr) {
    int count = 0;
    for (int i = 0; i < (int)str.size(); i++) {
        if (str[i] == chr) {
            count++;
        }
    }
    return count;
}

std::vector<std::string> pick_moms(unsigned int min_n2, unsigned int max_n2, std::string mom_set) {
    std::vector<std::string> momenta;
    std::vector<int> min_ind = {0,1,7,19,27};
    std::vector<int> max_ind = {1,7,19,27,(int)moms.size()};
    for (int i = min_ind[min_n2]; i < max_ind[max_n2]; i++) {
        if (mom_set == "all") { // all momenta -- only way to get terms with +ve and -ve values in the vector
            momenta.push_back(moms[i]);
        }
        else if (mom_set == "+ve") { // only all positive momenta
            int negs = count_char(moms[i],'-');
            if (negs == 0) {
                momenta.push_back(moms[i]);
            }
        }
        else if (mom_set == "-ve") { // only all negative momenta
            int ones = count_char(moms[i],'1');
            int negs = count_char(moms[i],'-');
            int zeros = count_char(moms[i],'0');
            if (ones == negs || zeros == 3) {
                momenta.push_back(moms[i]);
            }
        }
    }
    return momenta;
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
                                        unsigned int, step,
                                        std::string, source_loc);
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
                                        std::string, gamma_snk,
                                        std::string, rhq_impr,
                                        unsigned int, n2_min,
                                        unsigned int, n2_max,
                                        std::string, moms);
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

std::array<std::string, 2> split_pos(std::string pos_4d){
    std::string space = "";
    std::string time = "";
    int count = 0;
    char del = ' ';
    for (int i = 0; i< (int)pos_4d.size(); i++) {
        if (pos_4d[i] == del) {
            count++;
        }
        if (count < 3) {
            space += pos_4d[i];
        }
        else {
            if (pos_4d[i] != del) {
                time += pos_4d[i];
            }
        }
    }
    std::array<std::string, 2> spacetime = {space, time};
    return spacetime;
}

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

std::string removeSpaces(std::string word) {
    std::string newWord;
    for (int i = 0; i < word.length(); i++) {
        if (word[i] != ' ') {
            newWord += word[i];
        }
    }
    return newWord;
}

// dwf propagator --> add flag for loading from / saving to binary?
std::string make_dwfpropagator(Application &application, TestInputs::DwfPar &dwfPar, std::string src, std::string mom, int i){
    std::string solve_name = "CG_" + dwfPar.quark;
    if (i == 0) {
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

        application.createModule<MSolver::RBPrecCG>(solve_name,solverPar);
    }

    // generate propagator
    MFermion::GaugeProp::Par quarkPar;
    quarkPar.solver = solve_name;
    quarkPar.source = src;

    std::string prop_name = dwfPar.quark + "_p" + mom;
    application.createModule<MFermion::GaugeProp>(prop_name, quarkPar);
    return prop_name;
}

// rhq action + solver
std::string make_rhqsolver(Application &application, TestInputs::RhqPar &rhqPar){
    // action
    MAction::WilsonClover::Par actionPar;
    actionPar.gauge = "gauge"; // maybe generalise this somehow, if at all necessary
    actionPar.mass  = rhqPar.mass;
    actionPar.csw_r = rhqPar.cswr;
    actionPar.csw_t = rhqPar.cswt;
    WilsonAnisotropyCoefficients C_anisotropy;
    C_anisotropy.isAnisotropic = true;
    C_anisotropy.t_direction   = rhqPar.tdir;
    C_anisotropy.xi_0          = rhqPar.xi0;
    C_anisotropy.nu            = rhqPar.nu;
    actionPar.clover_anisotropy = C_anisotropy;
    actionPar.boundary = boundary;
    actionPar.twist    = twist;

    std::string act_name = "RHQ_action_" + rhqPar.quark;
    application.createModule<MAction::WilsonClover>(act_name, actionPar);

    // solver
    MSolver::RBPrecCG::Par solverPar;
    solverPar.action       = act_name;
    solverPar.residual     = rhqPar.residual;
    solverPar.maxIteration = rhqPar.maxiter;

    std::string solve_name = "CG_" + rhqPar.quark;
    application.createModule<MSolver::RBPrecCG>(solve_name,solverPar);
    return solve_name;
}

// rhq propagator 
std::string make_rhqpropagator(Application &application, TestInputs::RhqPar &rhqPar, std::string src, std::string solve_name, std::string add){
    // generate propagator
    MFermion::GaugeProp::Par quarkPar;
    quarkPar.solver = solve_name;
    quarkPar.source = src;

    std::string prop_name = rhqPar.quark + add; 
    application.createModule<MFermion::GaugeProp>(prop_name, quarkPar);
    return prop_name;
}

// rhq insertion I/II
std::string insert_rhqI(Application &application, std::string quark, std::string gamma, int index) {
    std::string insert_str = "RHQImprI_" + quark;
    MRHQ::RHQInsertionI::Par RHQImprI;
    RHQImprI.q = quark;
    RHQImprI.gauge = "gauge";
    RHQImprI.index = index;
    if (gamma == "Identity") {
        RHQImprI.gamma5 = Gamma::Algebra::Identity;
        insert_str += "_V";
    }
    else if (gamma == "Gamma5") {
        RHQImprI.gamma5 = Gamma::Algebra::Gamma5;
        insert_str += "_A";
    }

    application.createModule<MRHQ::RHQInsertionI>(insert_str, RHQImprI);
    return insert_str;
}

// rhq insertion III
std::string insert_rhqIII(Application &application, std::string quark, std::string gamma, int index) {
    std::string insert_str = "RHQImprIII_" + quark;
    MRHQ::RHQInsertionIII::Par RHQImprIII;
    RHQImprIII.q = quark;
    RHQImprIII.gauge = "gauge";
    RHQImprIII.index = index;
    RHQImprIII.flag = 0;
    if (gamma == "Identity") {
        RHQImprIII.gamma5 = Gamma::Algebra::Identity;
        insert_str += "_V";
    }
    else if (gamma == "Gamma5") {
        RHQImprIII.gamma5 = Gamma::Algebra::Gamma5;
        insert_str += "_A";
    }

    application.createModule<MRHQ::RHQInsertionIII>(insert_str, RHQImprIII);
    return insert_str;
}

// rhq insertion IV
std::string insert_rhqIV(Application &application, std::string quark, std::string gamma, int index) {
    std::string insert_str = "RHQImprIV_" + quark;
    MRHQ::RHQInsertionIV::Par RHQImprIV;
    RHQImprIV.q = quark;
    RHQImprIV.gauge = "gauge";
    RHQImprIV.index = index;
    RHQImprIV.flag = 0;
    if (gamma == "Identity") {
        RHQImprIV.gamma5 = Gamma::Algebra::Identity;
        insert_str += "_V";
    }
    else if (gamma == "Gamma5") {
        RHQImprIV.gamma5 = Gamma::Algebra::Gamma5;
        insert_str += "_A";
    }

    application.createModule<MRHQ::RHQInsertionIV>(insert_str, RHQImprIV);
    return insert_str;
}

// rhq insertion V
std::string insert_rhqV(Application &application, std::string quark, std::string gamma, int index) {
    std::string insert_str = "RHQImprV_" + quark;
    MRHQ::RHQInsertionV::Par RHQImprV;
    RHQImprV.q = quark;
    RHQImprV.gauge = "gauge";
    RHQImprV.index = index;
    if (gamma == "Identity") {
        RHQImprV.gamma5 = Gamma::Algebra::Identity;
        insert_str += "_V";
    }
    else if (gamma == "Gamma5") {
        RHQImprV.gamma5 = Gamma::Algebra::Gamma5;
        insert_str += "_A";
    }

    application.createModule<MRHQ::RHQInsertionV>(insert_str, RHQImprV);
    return insert_str;
}

// rhq insertion VI
std::string insert_rhqVI(Application &application, std::string quark, std::string gamma, int index) {
    std::string insert_str = "RHQImprVI_" + quark;
    MRHQ::RHQInsertionVI::Par RHQImprVI;
    RHQImprVI.q = quark;
    RHQImprVI.gauge = "gauge";
    RHQImprVI.index = index;
    if (gamma == "Identity") {
        RHQImprVI.gamma5 = Gamma::Algebra::Identity;
        insert_str += "_V";
    }
    else if (gamma == "Gamma5") {
        RHQImprVI.gamma5 = Gamma::Algebra::Gamma5;
        insert_str += "_A";
    }

    application.createModule<MRHQ::RHQInsertionVI>(insert_str, RHQImprVI);
    return insert_str;
}

// (phased) Z2 wall source
std::string make_Z2src(Application &application, std::string quark, int time){
    // start with standard Z2 source
    std::string z2_mom = "";
    MSource::Z2::Par z2Par;
    z2Par.tA = time; 
    z2Par.tB = time;
    z2_mom = "z2_p000_t0";
    application.createModule<MSource::Z2>(z2_mom, z2Par);
    return z2_mom;
}

std::string make_PhasedZ2(Application &application, std::string quark, std::string mom, int time){
    std::string z2_mom = "";
    MSource::MomentumPhase::Par momPar;
    momPar.src = "z2_p000_t0";
    momPar.mom = mom + " 0";
    std::string mom_str = mom;
    z2_mom = "z2_" + quark + "_p" + removeSpaces(mom_str) + "_t0";
    application.createModule<MSource::MomentumPhase>(z2_mom,momPar);
    return z2_mom;
}

// general contraction 
std::string make_contraction(Application &application, std::string meson, std::string q1, std::string q2, std::string gammas, std::string sink, std::string fileSt){

    MContraction::Meson::Par contraction;
    contraction.q1 = q1;
    contraction.q2 = q2;
    contraction.gammas = gammas;
    contraction.sink   = sink;

    contraction.output = fileSt + meson;

    application.createModule<MContraction::Meson>(meson, contraction);
    return fileSt + meson;
}

int main(int argc, char *argv[])
{
    //program expects parameter file as first argument of command line
    if (argc < 2) {
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

    std::vector<std::string> gamma_src_vec = split_gammas(testPar.mesonPar.gamma_src,' ');
    std::vector<std::string> gamma_snk_vec = split_gammas(testPar.mesonPar.gamma_snk,' ');
    std::vector<std::string> rhq_impr_vec = split_gammas(testPar.mesonPar.rhq_impr,' ');
    // program expects same number of snk and src gammas
    if (gamma_src_vec.size() != gamma_snk_vec.size()) {
        std::cerr << "gamma_src and gamma_snk must have same number of entries:" << std::endl
                  << "    gamma_src has " << gamma_src_vec.size() << " entries," << std::endl
                  << "    gamma_snk has " << gamma_snk_vec.size() << " entries." << std::endl;

        return EXIT_FAILURE;
    }
    std::string gamma_pairs = "";
    std::string gamma_ax = "";
    std::string gamma_vt = "";
    for (int i = 0; i < gamma_src_vec.size(); i++) {
        std::string pair = "(" + gamma_snk_vec[i] + " " + gamma_src_vec[i] + ")";
        gamma_pairs += pair;
        if (rhq_impr_vec[i] == "yes") {
            std::string pair_rhq = "(Identity " + gamma_src_vec[i] + ")";
            if (gamma_snk_vec[i].size() == 12) {
                gamma_ax += pair_rhq;}
            else {
                gamma_vt += pair_rhq;}
        }
    }

    // point source
    // MSource::Point::Par ptPar;
    // ptPar.position = testPar.configPar.source_loc;
    // application.createModule<MSource::Point>("pt_src",ptPar);

    std::array<std::string, 2> spacetime = split_pos(testPar.configPar.source_loc);
    int src_time = std::stoi(spacetime[1]);
    std::string src0 = make_Z2src(application, testPar.rhqPar.quark, src_time); 
    
    // *****************************
    // ***** RHQ quark modules *****
    // *****************************
        
    std::string b_solver = make_rhqsolver(application, testPar.rhqPar);
    std::string name_q2 = make_rhqpropagator(application, testPar.rhqPar, src0, b_solver, "");

    // *****************************
    // ***** RHQ RHQ insertion *****
    // *****************************

    // RHQ Improvement 2
    std::string O2_V = insert_rhqI(application, name_q2, "Identity", 3);
    std::string O2_A = insert_rhqI(application, name_q2, "Gamma5", 3);
    
    // RHQ Improvement 4
    std::string O4_V = insert_rhqIV(application, name_q2, "Identity", 3);
    std::string O4_A = insert_rhqIV(application, name_q2, "Gamma5", 3);
    
    // loop over all momenta up to max n2
    std::vector<std::string> momenta = pick_moms(testPar.mesonPar.n2_min, testPar.mesonPar.n2_max, testPar.mesonPar.moms);
    for (int i = 0; i < momenta.size(); i++) {
        // dwf source
        std::string src_mom = momenta[i];
        std::string src1 = "";
        if (src_mom == "0 0 0") {
            src1 = src0;
        }
        else {
            src1 = make_PhasedZ2(application, testPar.dwfPar.quark, src_mom, src_time);
        }

        std::string mes = testPar.rhqPar.quark + testPar.dwfPar.quark;
        std::string fileSt = mes + "_" + squareMom(src_mom) + "/";

        std::string snk_mom;
        // snk mom must be opposite sign of Z2 src mom
        if (src_mom != moms[0]) {
            int i = find_ind(src_mom);
            int j;
            if (i % 2 == 0) {
                j = i - 1;
            }
            else {
                j = i + 1;
            }
            snk_mom = moms[j];
        }
        else {
            snk_mom = moms[0];
        }
//        mes += "_p" + removeSpaces(snk_mom);
        MSink::Point::Par sinkPar;
        sinkPar.mom = snk_mom;
        std::string snk_str = "sink_p" + removeSpaces(snk_mom);
        application.createModule<MSink::ScalarPoint>(snk_str, sinkPar);
        
        // *****************************
        // ******* quark modules *******
        // *****************************
        
        // add for loop over different quark masses? --> same as gammas below
        std::string name_q1 = make_dwfpropagator(application, testPar.dwfPar, src1, removeSpaces(src_mom), i);

        // *****************************
        // ***** RHQ dwf insertion *****
        // *****************************

        // RHQ Improvement 1
        std::string O1_V = insert_rhqI(application, name_q1, "Identity", 3);
        std::string O1_A = insert_rhqI(application, name_q1, "Gamma5", 3);
    
        // RHQ Improvement 3
        std::string O3_V = insert_rhqIII(application, name_q1, "Identity", 3);
        std::string O3_A = insert_rhqIII(application, name_q1, "Gamma5", 3);

        std::string O1_Vmes = name_q2 + O1_V;
        std::string O2_Vmes = O2_V + name_q1;
        std::string O3_Vmes = name_q2 + O3_V;
        std::string O4_Vmes = O4_V + name_q1;
        std::string O1_Ames = name_q2 + O1_A;
        std::string O2_Ames = O2_A + name_q1;
        std::string O3_Ames = name_q2 + O3_A;
        std::string O4_Ames = O4_A + name_q1;

        // *****************************
        // ***** Meson Contraction *****
        // *****************************

        // O(0) operators
        make_contraction(application, name_q2 + name_q1 + "_Z2_PT", name_q1, name_q2, gamma_pairs, snk_str, fileSt);

        // O(a) operators Axial
        if (gamma_ax != "") {
            make_contraction(application, O1_Ames + "_Z2_PT", O1_A, name_q2, gamma_ax, snk_str, fileSt);
            make_contraction(application, O2_Ames + "_Z2_PT", name_q1, O2_A, gamma_ax, snk_str, fileSt);
            make_contraction(application, O3_Ames + "_Z2_PT", O3_A, name_q2, gamma_ax, snk_str, fileSt);
            make_contraction(application, O4_Ames + "_Z2_PT", name_q1, O4_A, gamma_ax, snk_str, fileSt);
        }

        // O(a) operators Vector
        if (gamma_vt != "") {
            make_contraction(application, O1_Vmes + "_Z2_PT", O1_V, name_q2, gamma_vt, snk_str, fileSt);
            make_contraction(application, O2_Vmes + "_Z2_PT", name_q1, O2_V, gamma_vt, snk_str, fileSt);
            make_contraction(application, O3_Vmes + "_Z2_PT", O3_V, name_q2, gamma_vt, snk_str, fileSt);
            make_contraction(application, O4_Vmes + "_Z2_PT", name_q1, O4_V, gamma_vt, snk_str, fileSt);
        }
    }

    // execution
    application.run();
    
    // epilogue
    LOG(Message) << "Grid is finalizing now" << std::endl;
    Grid_finalize();
    
    return EXIT_SUCCESS;
}
