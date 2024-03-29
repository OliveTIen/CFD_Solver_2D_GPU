#include "TomlFileManager.h"
#include "../global/GlobalPara.h"
#include <iostream>
#include "../output/LogWriter.h"
#include "AirParameterConverter.h"

// 类指针
TomlFileManager* TomlFileManager::classPointer = nullptr;

TomlFileManager* TomlFileManager::getInstance() {
    if (classPointer == nullptr) {
        classPointer = new TomlFileManager();
    }
    return classPointer;
}

// CPP读取TOML文件格式的方法 https://blog.csdn.net/jiaostyle/article/details/125695972
void TomlFileManager::readFileAndParseFile(std::string fullFilePath) {
    // 对于OptionalValue，不应该用异常处理，因为throw会中断try后面的代码，
    // 导致后面的变量没有被初始化
    try {
        m_parsedFile = cpptoml::parse_file(fullFilePath);// 若文件路径错误，或者内容拼写错误，则抛出异常
        modifyGlobalParametersAccordingToParsedFile();
        initialize_ruvp();
    }
    catch (cpptoml::parse_exception par) {
        LogWriter::logAndPrintError("Parse Toml file failed. Will exit. \n");
        std::cout << par.what() << std::endl;
        exit(1919810);
    }


    if (has_getOptionalValueFailed) {
        std::cout << "Warning: Get optional value failed, and they'll be ignored. Look into log for details.\n";
    }
    if (has_getValueFailed) {
        LogWriter::logAndPrintError("Error in reading input parameter. Will exit. \n");
        std::cout << "Please check the file \"input.toml\": " << "\n";
        std::cout << "1.Spelling mistake.\n";
        std::cout << "2.For boolean values, you should write \"true\" or \"false\", rather than \"1\" or \"0\".\n";
        std::cout << "(Code: " << __FILE__ << ")" << "\n";
        exit(114514);
    }

}

void TomlFileManager::printParsedFile() {
    std::cout << "File: \n" << *m_parsedFile << std::endl;
}

void TomlFileManager::modifyGlobalParametersAccordingToParsedFile() {
    getValue("basic.continue", GlobalPara::basic::_continue);
    getValue("basic.dimension", GlobalPara::basic::dimension);
    getValue("basic.filename", GlobalPara::basic::filename);
    getValue("basic.meshFileType", GlobalPara::basic::meshFileType);
    getValue("basic.useGPU", GlobalPara::basic::useGPU);
    //meshFileType = "inp"

    getValue("constant.T0", GlobalPara::constant::T0);
    getValue("constant.p0", GlobalPara::constant::p0);
    getValue("constant.Re", GlobalPara::constant::Re);
    getValue("constant.Pr", GlobalPara::constant::Pr);
    getValue("constant.gamma", GlobalPara::constant::gamma);

    getOptionalValueIfExists("boundaryCondition.2D.inf.AOA", GlobalPara::boundaryCondition::_2D::inf::AOA);
    getOptionalValueIfExists("boundaryCondition.2D.inf.Ma", GlobalPara::boundaryCondition::_2D::inf::Ma);
    getOptionalValueIfExists("boundaryCondition.2D.inf.p", GlobalPara::boundaryCondition::_2D::inf::ruvp[3]);
    getOptionalValueIfExists("boundaryCondition.2D.inf.rho", GlobalPara::boundaryCondition::_2D::inf::ruvp[0]);
    getOptionalValueIfExists("boundaryCondition.2D.inf.u", GlobalPara::boundaryCondition::_2D::inf::ruvp[1]);
    getOptionalValueIfExists("boundaryCondition.2D.inf.use_ruvp", GlobalPara::boundaryCondition::_2D::inf::use_ruvp);
    getOptionalValueIfExists("boundaryCondition.2D.inf.v", GlobalPara::boundaryCondition::_2D::inf::ruvp[2]);

    getOptionalValueIfExists("boundaryCondition.2D.inlet.AOA", GlobalPara::boundaryCondition::_2D::inlet::AOA);
    getOptionalValueIfExists("boundaryCondition.2D.inlet.Ma", GlobalPara::boundaryCondition::_2D::inlet::Ma);
    getOptionalValueIfExists("boundaryCondition.2D.inlet.p", GlobalPara::boundaryCondition::_2D::inlet::ruvp[3]);
    getOptionalValueIfExists("boundaryCondition.2D.inlet.rho", GlobalPara::boundaryCondition::_2D::inlet::ruvp[0]);
    getOptionalValueIfExists("boundaryCondition.2D.inlet.u", GlobalPara::boundaryCondition::_2D::inlet::ruvp[1]);
    getOptionalValueIfExists("boundaryCondition.2D.inlet.use_ruvp", GlobalPara::boundaryCondition::_2D::inlet::use_ruvp);
    getOptionalValueIfExists("boundaryCondition.2D.inlet.v", GlobalPara::boundaryCondition::_2D::inlet::ruvp[2]);

    getOptionalValueIfExists("boundaryCondition.2D.outlet.AOA", GlobalPara::boundaryCondition::_2D::outlet::AOA);
    getOptionalValueIfExists("boundaryCondition.2D.outlet.Ma", GlobalPara::boundaryCondition::_2D::outlet::Ma);
    getOptionalValueIfExists("boundaryCondition.2D.outlet.p", GlobalPara::boundaryCondition::_2D::outlet::ruvp[3]);
    getOptionalValueIfExists("boundaryCondition.2D.outlet.rho", GlobalPara::boundaryCondition::_2D::outlet::ruvp[0]);
    getOptionalValueIfExists("boundaryCondition.2D.outlet.u", GlobalPara::boundaryCondition::_2D::outlet::ruvp[1]);
    getOptionalValueIfExists("boundaryCondition.2D.outlet.use_ruvp", GlobalPara::boundaryCondition::_2D::outlet::use_ruvp);
    getOptionalValueIfExists("boundaryCondition.2D.outlet.v", GlobalPara::boundaryCondition::_2D::outlet::ruvp[2]);

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

    getValue("space.space_1D.nElement", GlobalPara::space::_1D::nElement);
    getValue("space.space_1D.x1", GlobalPara::space::_1D::x1);
    getValue("space.space_1D.x2", GlobalPara::space::_1D::x2);

    getValue("time.CFL", GlobalPara::time::CFL);
    getValue("time.T", GlobalPara::time::T);
    getValue("time.time_advance", GlobalPara::time::time_advance);
    getValue("time.residual", GlobalPara::time::residual);

    getValue("inviscid_flux_method.flag_reconstruct", GlobalPara::space::flag_reconstruct);
    getValue("inviscid_flux_method.flag_gradient", GlobalPara::space::flag_gradient);
    getValue("inviscid_flux_method.flux_conservation_scheme", GlobalPara::inviscid_flux_method::flux_conservation_scheme);
    getValue("inviscid_flux_method.flux_limiter", GlobalPara::inviscid_flux_method::flux_limiter);

    //if (has_getValueFailed) {
    //    std::cout << "Please check the file \"input.toml\": " << "\n";
    //    std::cout << "1.Spelling mistake.\n";
    //    std::cout << "2.For boolean values, you should write \"true\" or \"false\", rather than \"1\" or \"0\".\n";
    //    std::cout << "(Code: " << __FILE__ << ")" << "\n";
    //    LogWriter::logAndPrint("Error in reading input parameter. Will exit. \n");
    //    exit(114514);
    //}
}

void TomlFileManager::initialize_ruvp() {
    // 处理部分输入参数，根据Ma和AOA计算ruvp
    //GlobalPara::boundaryCondition::_2D::ini_ruvp_by_Ma_AOA();
        //若不给定ruvp，则根据Ma和AOA计算ruvp
    using namespace GlobalPara::boundaryCondition::_2D;
    if (!inf::use_ruvp) {
        //std::cout << "inf: Ma = " << inf::Ma << ", "
        //    << "AOA = " << inf::AOA << "\n";
        AirParameterConverter::get_ruvp_by_Ma_AoA_T0_p0_gamma_R(
            inf::ruvp, inf::Ma, inf::AOA, GlobalPara::constant::T0, GlobalPara::constant::p0, GlobalPara::constant::gamma, GlobalPara::constant::R);
    }
    if (!inlet::use_ruvp) {
        AirParameterConverter::get_ruvp_by_Ma_AoA_T0_p0_gamma_R(
            inlet::ruvp, inlet::Ma, inlet::AOA, GlobalPara::constant::T0, GlobalPara::constant::p0, GlobalPara::constant::gamma, GlobalPara::constant::R);
    }
    if (!outlet::use_ruvp) {
        AirParameterConverter::get_ruvp_by_Ma_AoA_T0_p0_gamma_R(
            outlet::ruvp, outlet::Ma, outlet::AOA, GlobalPara::constant::T0, GlobalPara::constant::p0, GlobalPara::constant::gamma, GlobalPara::constant::R);
    }


}

void TomlFileManager::logErrorMessage(ErrorMessage errorMessage) {
    LogWriter::log(errorMessage.getMessage() + "\n");
}

