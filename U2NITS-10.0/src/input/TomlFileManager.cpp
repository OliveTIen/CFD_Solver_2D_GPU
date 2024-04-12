#include "TomlFileManager.h"
#include "../global/GlobalPara.h"
#include <iostream>
#include "../output/LogWriter.h"
#include "AirParameterConverter.h"
#include "../math/PhysicalKernel.h"
#include "../global/StringProcessor.h"

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
        exit(114514);
    }

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
        system("pause");
        exit(-1);
    }

    if (has_getValueFailed) {
        exit(114514);
    }

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
    getValue("constant.Re", GlobalPara::constant::Re);
    getValue("constant.Pr", GlobalPara::constant::Pr);
    getValue("constant.gamma", GlobalPara::constant::gamma);

    // 旧 和 handleConflictingInputs();配套使用
    //getValueIfExists("boundaryCondition.2D.inf.input_mode", GlobalPara::boundaryCondition::_2D::inf::input_mode);
    //getValueIfExists("boundaryCondition.2D.inf.rho", GlobalPara::boundaryCondition::_2D::inf::ruvp[0]);
    //getValueIfExists("boundaryCondition.2D.inf.u", GlobalPara::boundaryCondition::_2D::inf::ruvp[1]);
    //getValueIfExists("boundaryCondition.2D.inf.v", GlobalPara::boundaryCondition::_2D::inf::ruvp[2]);
    //getValueIfExists("boundaryCondition.2D.inf.p", GlobalPara::boundaryCondition::_2D::inf::ruvp[3]);
    //getValueIfExists("boundaryCondition.2D.inf.Ma", GlobalPara::boundaryCondition::_2D::inf::Ma);
    //getValueIfExists("boundaryCondition.2D.inf.AOA", GlobalPara::boundaryCondition::_2D::inf::AOA);

    //getValueIfExists("boundaryCondition.2D.inlet.input_mode", GlobalPara::boundaryCondition::_2D::inlet::input_mode);
    //getValueIfExists("boundaryCondition.2D.inlet.rho", GlobalPara::boundaryCondition::_2D::inlet::ruvp[0]);
    //getValueIfExists("boundaryCondition.2D.inlet.u", GlobalPara::boundaryCondition::_2D::inlet::ruvp[1]);
    //getValueIfExists("boundaryCondition.2D.inlet.v", GlobalPara::boundaryCondition::_2D::inlet::ruvp[2]);
    //getValueIfExists("boundaryCondition.2D.inlet.p", GlobalPara::boundaryCondition::_2D::inlet::ruvp[3]);
    //getValueIfExists("boundaryCondition.2D.inlet.Ma", GlobalPara::boundaryCondition::_2D::inlet::Ma);
    //getValueIfExists("boundaryCondition.2D.inlet.AOA", GlobalPara::boundaryCondition::_2D::inlet::AOA);

    //getValueIfExists("boundaryCondition.2D.outlet.input_mode", GlobalPara::boundaryCondition::_2D::outlet::input_mode);
    //getValueIfExists("boundaryCondition.2D.outlet.rho", GlobalPara::boundaryCondition::_2D::outlet::ruvp[0]);
    //getValueIfExists("boundaryCondition.2D.outlet.u", GlobalPara::boundaryCondition::_2D::outlet::ruvp[1]);
    //getValueIfExists("boundaryCondition.2D.outlet.v", GlobalPara::boundaryCondition::_2D::outlet::ruvp[2]);
    //getValueIfExists("boundaryCondition.2D.outlet.p", GlobalPara::boundaryCondition::_2D::outlet::ruvp[3]);
    //getValueIfExists("boundaryCondition.2D.outlet.Ma", GlobalPara::boundaryCondition::_2D::outlet::Ma);
    //getValueIfExists("boundaryCondition.2D.outlet.AOA", GlobalPara::boundaryCondition::_2D::outlet::AOA);

    // 新添加
    getValue_boundaryCondition2D("boundaryCondition.2D.inf", GlobalPara::boundaryCondition::_2D::inf::ruvp);
    getValue_boundaryCondition2D("boundaryCondition.2D.inlet", GlobalPara::boundaryCondition::_2D::inlet::ruvp);
    getValue_boundaryCondition2D("boundaryCondition.2D.outlet", GlobalPara::boundaryCondition::_2D::outlet::ruvp);

    getValue("initialCondition.type", GlobalPara::initialCondition::type);

    getValue("output.step_per_print", GlobalPara::output::step_per_print);
    getValue("output.step_per_output_field", GlobalPara::output::step_per_output_field);
    getValue("output.step_per_output_hist", GlobalPara::output::step_per_output_hist);
    getValue("output.maxIteration", GlobalPara::output::maxIteration);

    getValue("output.output_var.rho", GlobalPara::output::output_var_ruvp[0]);
    getValue("output.output_var.u", GlobalPara::output::output_var_ruvp[1]);
    getValue("output.output_var.v", GlobalPara::output::output_var_ruvp[2]);
    getValue("output.output_var.p", GlobalPara::output::output_var_ruvp[3]);

    getValue("physicsModel.equation", GlobalPara::physicsModel::equation);

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

void TomlFileManager::getValue_boundaryCondition2D(std::string parent, double ruvp[4]) {
    /*
    事实上input_mode Ma AoA无需作为全局变量，因为仅仅在输入时用到
    */
    // eg. parent = "boundaryCondition.2D.inf"
    int input_mode = -1;
    double Ma = 0;
    double AoA = 0;
    double U[4];
    double gamma = 1.4;
    if (!treeContainsKey(parent + ".input_mode")) {
        LogWriter::logAndPrintWarning("No input_mode at " + parent + ". Will use default.\n");
        return;
    }
    getValueIfExists(parent + ".input_mode", input_mode);
    const std::string description = "0- Ma,AoA \n1- rho,u,v,p \n2- rho, rhou, rhov, rhoE \n";
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

