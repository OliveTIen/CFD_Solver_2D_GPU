#include "TomlFileManager.h"
#include "../GlobalPara.h"
#include <iostream>
#include "../output/LogWriter.h"
#include "AirParameterConverter.h"

// CPP读取TOML文件格式的方法 https://blog.csdn.net/jiaostyle/article/details/125695972

void TomlFileManager::readFileAndParseFile(std::string fullFilePath) {
    try {
        m_parsedFile = cpptoml::parse_file(fullFilePath);// 若文件路径错误，或者内容拼写错误，则抛出异常
        modifyGlobalParametersAccordingToParsedFile();
        initialize_ruvp();
    }
    catch (cpptoml::parse_exception par) {
        // 异常处理
        std::cout << "Toml File Parse Error: ";
        std::cout << par.what() << std::endl;
        std::cout << "(File path: " << fullFilePath << ")" << "\n";
        std::cout << "(Code: " << __FILE__ << ")" << "\n";
    }
}

void TomlFileManager::printParsedFile() {
    std::cout << "File: \n" << *m_parsedFile << std::endl;
}

void TomlFileManager::modifyGlobalParametersAccordingToParsedFile() {
    getValueIfExists("basic.continue", GlobalPara::basic::_continue);
    getValueIfExists("basic.dimension", GlobalPara::basic::dimension);
    getValueIfExists("basic.filename", GlobalPara::basic::filename);

    getValueIfExists("constant.T0", Constant::T0);
    getValueIfExists("constant.p0", Constant::p0);

    getValueIfExists("boundaryCondition.2D.inf.AOA", GlobalPara::boundaryCondition::_2D::inf::AOA);
    getValueIfExists("boundaryCondition.2D.inf.Ma", GlobalPara::boundaryCondition::_2D::inf::Ma);
    getValueIfExists("boundaryCondition.2D.inf.p", GlobalPara::boundaryCondition::_2D::inf::ruvp[3]);
    getValueIfExists("boundaryCondition.2D.inf.rho", GlobalPara::boundaryCondition::_2D::inf::ruvp[0]);
    getValueIfExists("boundaryCondition.2D.inf.u", GlobalPara::boundaryCondition::_2D::inf::ruvp[1]);
    getValueIfExists("boundaryCondition.2D.inf.use_ruvp", GlobalPara::boundaryCondition::_2D::inf::use_ruvp);
    getValueIfExists("boundaryCondition.2D.inf.v", GlobalPara::boundaryCondition::_2D::inf::ruvp[2]);

    getValueIfExists("boundaryCondition.2D.inlet.AOA", GlobalPara::boundaryCondition::_2D::inlet::AOA);
    getValueIfExists("boundaryCondition.2D.inlet.Ma", GlobalPara::boundaryCondition::_2D::inlet::Ma);
    getValueIfExists("boundaryCondition.2D.inlet.p", GlobalPara::boundaryCondition::_2D::inlet::ruvp[3]);
    getValueIfExists("boundaryCondition.2D.inlet.rho", GlobalPara::boundaryCondition::_2D::inlet::ruvp[0]);
    getValueIfExists("boundaryCondition.2D.inlet.u", GlobalPara::boundaryCondition::_2D::inlet::ruvp[1]);
    getValueIfExists("boundaryCondition.2D.inlet.use_ruvp", GlobalPara::boundaryCondition::_2D::inlet::use_ruvp);
    getValueIfExists("boundaryCondition.2D.inlet.v", GlobalPara::boundaryCondition::_2D::inlet::ruvp[2]);

    getValueIfExists("boundaryCondition.2D.outlet.AOA", GlobalPara::boundaryCondition::_2D::outlet::AOA);
    getValueIfExists("boundaryCondition.2D.outlet.Ma", GlobalPara::boundaryCondition::_2D::outlet::Ma);
    getValueIfExists("boundaryCondition.2D.outlet.p", GlobalPara::boundaryCondition::_2D::outlet::ruvp[3]);
    getValueIfExists("boundaryCondition.2D.outlet.rho", GlobalPara::boundaryCondition::_2D::outlet::ruvp[0]);
    getValueIfExists("boundaryCondition.2D.outlet.u", GlobalPara::boundaryCondition::_2D::outlet::ruvp[1]);
    getValueIfExists("boundaryCondition.2D.outlet.use_ruvp", GlobalPara::boundaryCondition::_2D::outlet::use_ruvp);
    getValueIfExists("boundaryCondition.2D.outlet.v", GlobalPara::boundaryCondition::_2D::outlet::ruvp[2]);

    getValueIfExists("initialCondition.type", GlobalPara::initialCondition::type);

    getValueIfExists("output.step_per_print", GlobalPara::output::step_per_print);
    getValueIfExists("output.step_per_output", GlobalPara::output::step_per_output);
    getValueIfExists("output.step_per_output_hist", GlobalPara::output::step_per_output_hist);

    getValueIfExists("output.output_var.rho", GlobalPara::output::output_var_ruvp[0]);
    getValueIfExists("output.output_var.u", GlobalPara::output::output_var_ruvp[1]);
    getValueIfExists("output.output_var.v", GlobalPara::output::output_var_ruvp[2]);
    getValueIfExists("output.output_var.p", GlobalPara::output::output_var_ruvp[3]);

    getValueIfExists("physicsModel.equation", GlobalPara::physicsModel::equation);
    getValueIfExists("physicsModel.gamma", Constant::gamma); // ! [todo]未统一

    getValueIfExists("space.space_1D.nElement", GlobalPara::space::_1D::nElement);
    getValueIfExists("space.space_1D.x1", GlobalPara::space::_1D::x1);
    getValueIfExists("space.space_1D.x2", GlobalPara::space::_1D::x2);

    getValueIfExists("time.CFL", GlobalPara::time::CFL);
    getValueIfExists("time.T", GlobalPara::time::T);

    if (has_error_on_modifying) {
        std::cout << "Please check the file \"input.toml\": " << "\n";
        std::cout << "1.Spelling mistake.\n";
        std::cout << "2.For boolean values, you should write \"true\" or \"false\", rather than \"1\" or \"0\".\n";
        std::cout << "(Code: " << __FILE__ << ")" << "\n";
        LogWriter::writeLogAndCout("Error in reading input parameter. Will exit. \n");
        exit(114514);
    }
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
            inf::ruvp, inf::Ma, inf::AOA, Constant::T0, Constant::p0, Constant::gamma, Constant::R);
    }
    if (!inlet::use_ruvp) {
        AirParameterConverter::get_ruvp_by_Ma_AoA_T0_p0_gamma_R(
            inlet::ruvp, inlet::Ma, inlet::AOA, Constant::T0, Constant::p0, Constant::gamma, Constant::R);
    }
    if (!outlet::use_ruvp) {
        AirParameterConverter::get_ruvp_by_Ma_AoA_T0_p0_gamma_R(
            outlet::ruvp, outlet::Ma, outlet::AOA, Constant::T0, Constant::p0, Constant::gamma, Constant::R);
    }


}

