#include "CInput.h"
#include "FieldInitializer.h"
#include "InpMeshReader.h"
#include "SU2MeshReader.h"
#include "TomlFileManager.h"
#include "../legacy/FVM_2D.h"
#include "../global/FilePathManager.h"
#include "../global/GlobalPara.h"
#include "../global/SystemInfo.h"
#include "../output/LogWriter.h"
#include "../output/ConsolePrinter.h"
#include "ContinueFileReader.h"
#include "../global/CExit.h"
#include "../global/StringProcessor.h"
#include <sstream>
#include <map>
#include "../drivers/RoutineController.h"


/*
检查输入值是否属于列表
value: 待检测值
category: 值所属类型
table: 列表，包含已知值和描述
b_print: 是否输出到屏幕
*/
void check_if_value_belongs_to_table(int value, std::string category, std::map<int,std::string> table, bool b_print = false) {
	/*
	2024-05-01
	检查输入值是否属于列表
	若属于，则输出对应的描述，否则报错退出
	列表为空时，相当于“不属于”
	*/
	std::stringstream ss;
	auto iter = table.find(value);
	if (iter != table.end()) {
		ss << category << ": " << iter->second << "\n";
		LogWriter::log(ss.str());
		if (b_print)LogWriter::print(ss.str());
	}
	else {
		ss << "invalid " << category << ": " << value << "\n";
		LogWriter::logError(ss.str());
		exit(-1);
	}
	
}

void U2NITS::CInput::readConfig() {
	LogWriter::log(SystemInfo::getCurrentDateTime() + "\n");
	TomlFileManager::getInstance()->readTomlFile(FilePathManager::getInstance()->getTomlFilePath());

	check_if_value_belongs_to_table(RoutineController::getInstance()->getStrategy(), "time strategy",
		{ {0, "null"}, {1, "adaptive CFL"} }, true);
	check_if_value_belongs_to_table(GlobalPara::basic::dimension, "dimension", { {2, "2D"} });
	check_if_value_belongs_to_table(GlobalPara::physicsModel::equation, "equation", { {_EQ_euler, "Euler"}, {_EQ_NS, "NS"} }, true);
	check_if_value_belongs_to_table(GlobalPara::basic::useGPU, "platform", { {0,"CPU(host)"},{1,"GPU(device)"},{2,"half GPU"},{3,"development mode"} }, true);
}


// 记录边界参数
void log_boundary_condition() {
	using namespace GlobalPara::boundaryCondition::_2D;
	myfloat* inlet_ruvp = inlet::ruvp;
	myfloat* outlet_ruvp = outlet::ruvp;
	myfloat* inf_ruvp = inf::ruvp;
	const int num_ruvp = 4;
	// 边界参数
	std::string str;
	str += "BoundaryCondition:\n";
	str += "inlet::ruvp\t" + StringProcessor::doubleArray_2_string(inlet_ruvp, num_ruvp)
		+ "\noutlet::ruvp\t" + StringProcessor::doubleArray_2_string(outlet_ruvp, num_ruvp)
		+ "\ninf::ruvp\t" + StringProcessor::doubleArray_2_string(inf_ruvp, num_ruvp)
		+ "\n";
	LogWriter::log(str);
	// 输出雷诺数参见FieldInitialize，因为要根据initial_type确定用inf, inlet还是outlet


}

void read_mesh_file_and_initialize() {
	/*
	从头开始。读取网格，初始化流场和边界
	*/
	// 读取网格
	std::string dir = FilePathManager::getInstance()->getInputDirectory();// 目录
	std::string basename = GlobalPara::basic::filename;// 文件名主体
	std::string type = GlobalPara::basic::meshFileType;// 后缀
	if (type == "su2") {
		std::string filename = basename + "." + type;
		LogWriter::log("mesh file: " + filename + "\n");
		SU2MeshReader::read_mesh_and_process(dir + filename, true);
	}
	else {
		LogWriter::logAndPrint("Invalid mesh file type: " + type + " @CInput::readField_1()\n");
		exit(-1);
	}
	// 初始化流场和边界
	FieldInitializer::getInstance()->setInitialAndBoundaryCondition();
}

void U2NITS::CInput::readField_1() {

	/*
	若Config要求续算，则直接读取续算文件，否则从头开始
	*/
	if (GlobalPara::basic::_continue) {
		// 若找不到续算文件，则从头开始
		if (ContinueFileReader::readContinueFile_1() == -1) {
			LogWriter::logAndPrint("failed to read continue file(pause_*.dat), will read mesh file.\n");
			read_mesh_file_and_initialize();
			GlobalPara::basic::_continue = false;// 保证_continue为false，防止后面histWriter不写文件头
		}
	}
	else {
		read_mesh_file_and_initialize();
	}

	// 日志记录边界参数
	log_boundary_condition();
}

void U2NITS::CInput::readField_2_unused() {
	/*
	读取网格文件，
	法一：读到单元数量那一行后，就可以开辟solver中单元的内存空间。缺点是你无法保证这个数字是否正确
	法二：先用std::vector一个一个push。完成后，开辟solver数据的内存并初始化，然后删除临时的std::vector
	*/
	LogWriter::logAndPrintError("unimplemented error @CInput::readField_2_unused\n");
	exit(-1);
	
}
