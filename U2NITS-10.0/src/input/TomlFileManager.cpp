#include "TomlFileManager.h"
#include "../global/GlobalPara.h"
#include <iostream>
#include "../output/LogWriter.h"
#include "AirParameterConverter.h"
#include "../math/PhysicalKernel.h"
#include "../global/StringProcessor.h"
#include "../global/CExit.h"
#include "../drivers/RoutineController.h"
#include "../output/FieldWriter.h"
#include "FieldInitializer.h"

// 类指针
TomlFileManager* TomlFileManager::classPointer = nullptr;

TomlFileManager* TomlFileManager::getInstance() {
    if (classPointer == nullptr) {
        classPointer = new TomlFileManager();
    }
    return classPointer;
}

// CPP读取TOML文件格式的方法 https://blog.csdn.net/jiaostyle/article/details/125695972
void TomlFileManager::readTomlFile(std::string fullFilePath) {
    // 对于OptionalValue，不应该用异常处理，因为throw会中断try后面的代码，
    // 导致后面的变量没有被初始化
    try {
        m_tree = cpptoml::parse_file(fullFilePath);// 若文件路径错误，或者内容拼写错误，则抛出异常
        // 转globalpara
        treeToGlobalParameter();
        //handleConflictingInputs();
    }
    catch (cpptoml::parse_exception par) {
        std::stringstream ss;
        ss << "parse_exception: " << par.what() << " @TomlFileManager::readTomlFile\n";
        LogWriter::logAndPrintError(ss.str());
        CExit::saveAndExit(114514);
    }

    ifFailedThenExit();

}


inline double get_velocity_magnitude(myfloat* ruvp) {
    myfloat u = ruvp[1];
    myfloat v = ruvp[2];
    return sqrt(u * u + v * v);
}

// 获取参考来流参数
void get_reference_value_according_to_initialType(myfloat& rho_ref, myfloat& U_ref, myfloat& L_ref, int initial_type) {
    /*
    对于激波管和双马赫反射，取inlet和outlet的速度大者作为来流参数
    其他情况，取inf作为来流参数
    */
    using namespace GlobalPara::boundaryCondition::_2D;
    if (initial_type == type_shock_tube || initial_type == type_double_mach_reflection) {
        myfloat U_inlet = get_velocity_magnitude(inlet::ruvp);
        myfloat U_outlet = get_velocity_magnitude(outlet::ruvp);
        if (U_inlet > U_outlet) {
            rho_ref = inlet::ruvp[0];
            U_ref = U_inlet;
        }
        else {
            rho_ref = outlet::ruvp[0];
            U_ref = U_outlet;
        }
    }
    else {
        rho_ref = inf::ruvp[0];
        U_ref = get_velocity_magnitude(inf::ruvp);
    }
    L_ref = GlobalPara::constant::referenceArea;// 二维，无需开根号
}

// 从toml读取mu0和Re。读取1个，另一个根据来流参数自动确定。依赖于GlobalPara::referencearea
void read_mu0_Re_from_config() {
    TomlFileManager* t = TomlFileManager::getInstance();

    myfloat rho_ref{ 1.0 }, U_ref{ 1.0 }, L_ref{ 1.0 };
    get_reference_value_according_to_initialType(rho_ref, U_ref, L_ref, FieldInitializer::getInstance()->get_initial_type());// 获取参考来流参数

    bool calculate_mu0_by_Re = false;
    t->getValueIfExists("constant.calculate_mu0_by_Re", calculate_mu0_by_Re);
    if (calculate_mu0_by_Re) {
        t->getValue("constant.Re", GlobalPara::constant::Re);
        GlobalPara::constant::mu0 = rho_ref * U_ref * L_ref / GlobalPara::constant::Re;
    }
    else {
        t->getValueIfExists("constant.mu0", GlobalPara::constant::mu0);
        GlobalPara::constant::Re = rho_ref * U_ref * L_ref / GlobalPara::constant::mu0;
    }

    std::stringstream ss;
    ss << std::scientific
        << "modify_mu0_by_Re_using_config: "
        << "calculate_mu0_by_Re=" << calculate_mu0_by_Re << ", "
        << "Re=" << GlobalPara::constant::Re << ", "
        << "mu0=" << GlobalPara::constant::mu0 << "\n";
    LogWriter::logAndPrint(ss.str());
    //CExit::pressAnyKeyToExit();
}

void read_initialCondition_from_config() {
    FieldInitializer::getInstance()->initialize_using_config();
}

void TomlFileManager::treeToGlobalParameter() {
    getValue("basic.continue", GlobalPara::basic::_continue);
    getValue("basic.dimension", GlobalPara::basic::dimension);
    getValue("basic.filename", GlobalPara::basic::filename);
    getValue("basic.meshFileType", GlobalPara::basic::meshFileType);
    getValue("basic.useGPU", GlobalPara::basic::useGPU);
    getValue("basic.isDebugMode", GlobalPara::basic::isDebugMode);

    getValue("constant.T0", GlobalPara::constant::T0);
    getValue("constant.p0", GlobalPara::constant::p0);
    getValue("constant.Pr", GlobalPara::constant::Pr);
    getValue("constant.gamma", GlobalPara::constant::gamma);
    getValue("constant.referenceArea", GlobalPara::constant::referenceArea);
    getValueIfExists("constant.mesh_scale_factor", GlobalPara::constant::mesh_scale_factor);
    // 应在读取referenceArea和mesh_scale_factor后修正referenceArea
    GlobalPara::constant::referenceArea *= GlobalPara::constant::mesh_scale_factor;

    getValue_boundaryCondition2D("boundaryCondition.2D.inf", GlobalPara::boundaryCondition::_2D::inf::ruvp);
    getValue_boundaryCondition2D("boundaryCondition.2D.inlet", GlobalPara::boundaryCondition::_2D::inlet::ruvp);
    getValue_boundaryCondition2D("boundaryCondition.2D.outlet", GlobalPara::boundaryCondition::_2D::outlet::ruvp);

    read_initialCondition_from_config();
    read_mu0_Re_from_config();// 应放在referenceArea、initial_type后面

    getValue("output.step_per_print", GlobalPara::output::step_per_print);
    getValue("output.step_per_output_field", GlobalPara::output::step_per_output_field);
    getValue("output.step_per_output_hist", GlobalPara::output::step_per_output_hist);
    getValueIfExists("output.start_output_field", GlobalPara::output::start_output_field);
    if (GlobalPara::output::start_output_field != 0) {
        LogWriter::log("start_output_field = " + std::to_string(GlobalPara::output::start_output_field) + "\n");
    }
    getValue("output.maxIteration", GlobalPara::output::maxIteration);
    getValue("output.tolerace_residual", GlobalPara::output::tolerace_residual);

    FieldWriter::getInstance()->initialize_outputScheme_usingConfig();

    getValue("physicsModel.equation", GlobalPara::physicsModel::equation);

    int time_strategy = 0;
    getValueIfExists("time.strategy", time_strategy);
    RoutineController::getInstance()->setStrategy(time_strategy);

    getValue("time.is_steady", GlobalPara::time::is_steady);
    getValueOnCondition("time.CFL_steady", GlobalPara::time::CFL_steady, GlobalPara::time::is_steady);
    getValueOnCondition("time.CFL", GlobalPara::time::CFL, !GlobalPara::time::is_steady);

    getValue("time.max_physical_time", GlobalPara::time::max_physical_time);
    getValue("time.time_advance", GlobalPara::time::time_advance);

    getValue("inviscid_flux_method.flag_reconstruct", GlobalPara::inviscid_flux_method::flag_reconstruct);
    getValue("inviscid_flux_method.flag_gradient", GlobalPara::inviscid_flux_method::flag_gradient);
    getValue("inviscid_flux_method.flux_conservation_scheme", GlobalPara::inviscid_flux_method::flux_conservation_scheme);
    getValue("inviscid_flux_method.flux_limiter", GlobalPara::inviscid_flux_method::flux_limiter);

}

void TomlFileManager::treeToSGlobalPara(SGlobalPara::SGlobalPara& spara) {
    // [未完成]建议放在CInput，因为后面SetInitialCondition时会修改GlobalPara
    LogWriter::logAndPrintError("Unimplemented error @TomlFileManager::treeToSGlobalPara");
}

void TomlFileManager::handleConflictingInputs() {
    // 处理部分输入参数，根据Ma和AOA计算ruvp
    
    // ruvp 
    using namespace GlobalPara::boundaryCondition::_2D;
    if (!inf::input_mode) {
        //std::cout << "inf: Ma = " << inf::Ma << ", "
        //    << "AOA = " << inf::AOA << "\n";
        AirParameterConverter::get_ruvp_by_Ma_AoA_T0_p0_gamma_R(
            inf::ruvp, inf::Ma, inf::AOA, GlobalPara::constant::T0, GlobalPara::constant::p0, GlobalPara::constant::gamma, GlobalPara::constant::R);
    }
    if (!inlet::input_mode) {
        AirParameterConverter::get_ruvp_by_Ma_AoA_T0_p0_gamma_R(
            inlet::ruvp, inlet::Ma, inlet::AOA, GlobalPara::constant::T0, GlobalPara::constant::p0, GlobalPara::constant::gamma, GlobalPara::constant::R);
    }
    if (!outlet::input_mode) {
        AirParameterConverter::get_ruvp_by_Ma_AoA_T0_p0_gamma_R(
            outlet::ruvp, outlet::Ma, outlet::AOA, GlobalPara::constant::T0, GlobalPara::constant::p0, GlobalPara::constant::gamma, GlobalPara::constant::R);
    }


}

void TomlFileManager::getValue_boundaryCondition2D(std::string parent, myfloat ruvp[4]) {
    /*
    事实上input_mode Ma AoA无需作为全局变量，因为仅仅在输入时用到
    */
    // eg. parent = "boundaryCondition.2D.inf"
    int input_mode = -1;
    myfloat Ma = 0;
    myfloat AoA = 0;
    myfloat U[4]{};
    myfloat gamma = 1.4;
    if (!treeContainsKey(parent + ".input_mode")) {
        LogWriter::logWarning("No input_mode at " + parent + ". Will use default.\n");
        return;
    }
    getValueIfExists(parent + ".input_mode", input_mode);
    const std::string description = "0- Ma,AoA \n1- rho,u,v,p \n2- rho, rhou, rhov, rhoE \n3- rho, u, angle_degree, p\n";
    switch (input_mode) {

    case 0:
        if (treeContainsKey(parent + ".AOA")) {
            getValue(parent + ".AOA", AoA);
        }
        else {
            getValue(parent + ".AoA", AoA);
        }
        getValue(parent + ".Ma", Ma);
        AirParameterConverter::get_ruvp_by_Ma_AoA_T0_p0_gamma_R(ruvp, Ma, AoA, GlobalPara::constant::T0, GlobalPara::constant::p0, GlobalPara::constant::gamma, GlobalPara::constant::R);
        
        break;
    case 1:
        getValue(parent + ".rho", ruvp[0]);
        getValue(parent + ".u", ruvp[1]);
        getValue(parent + ".v", ruvp[2]);
        getValue(parent + ".p", ruvp[3]);

        break;
    case 2:
        getValue(parent + ".rho", U[0]);
        getValue(parent + ".rhou", U[1]);
        getValue(parent + ".rhov", U[2]);
        getValue(parent + ".rhoE", U[3]);
        getValueIfExists("constant.gamma", gamma);
        U2NITS::Math::U2ruvp_host(U, ruvp, gamma);

        break;
    case 3:
        {
        myfloat u = 0.0;
        myfloat angle = 0.0;
        getValue(parent + ".rho", ruvp[0]);
        getValue(parent + ".u", u);
        getValue(parent + ".angle_degree", angle);
        getValue(parent + ".p", ruvp[3]);
        angle = angle * U2NITS::Math::PI / 180.0;
        ruvp[1] = u * cos(angle);
        ruvp[2] = u * sin(angle);
        }
        break;

    default:
        LogWriter::logAndPrintError("Invalid input mode at " + parent + ".\n" + description);
        has_getValueFailed = true;
        break;
    }
   
}

void TomlFileManager::printTreeContainsKeyOrNot(std::string key) {
    if (treeContainsKey(key)) {
        std::cout << "Tree contains key \"" << key << "\"\n";
    }
    else {
        std::cout << "Tree doesn't contain key \"" << key << "\"\n";
    }
}

void TomlFileManager::ifFailedThenExit() {
    if (has_getValueFailed) {
        std::stringstream ss;
        ss << "read_error_exception @TomlFileManager::readTomlFile. ";
        ss << "Please check \"input.toml\" for: " << "\n";
        ss << "1.Spelling mistake.\n";
        ss << "2.Invalid value type.\n";
        ss << " @TomlFileManager::readTomlFile \n";
        LogWriter::logAndPrintError(ss.str());
    }

    bool isDebugMode = GlobalPara::basic::isDebugMode;
    if (isDebugMode) {
        std::cout << "Final ruvp: \n";
        std::cout << "  inf: \t" << StringProcessor::doubleArray_2_string(GlobalPara::boundaryCondition::_2D::inf::ruvp, 4);
        std::cout << "\n  inlet: \t" << StringProcessor::doubleArray_2_string(GlobalPara::boundaryCondition::_2D::inlet::ruvp, 4);
        std::cout << "\n  outlet: \t" << StringProcessor::doubleArray_2_string(GlobalPara::boundaryCondition::_2D::outlet::ruvp, 4);

        if (has_getValueFailed) {
            std::cout << "has_getValueFailed. Will exit.\n";
        }
        printTree();
    }

    if (isDebugMode) {
        
        CExit::pressAnyKeyToExit();
    }

    if (has_getValueFailed) {
        CExit::saveAndExit(114514);
    }
}

