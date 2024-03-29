#include "CInput.h"
#include "FieldInitializer.h"
#include "InpMeshReader.h"
#include "SU2MeshReader.h"
#include "TomlFileManager.h"
#include "../FVM_2D.h"
#include "../global/FilePathManager.h"
#include "../global/GlobalPara.h"
#include "../global/SystemInfo.h"
#include "../output/LogWriter.h"
#include "../output/ConsolePrinter.h"
#include "ContinueFileReader.h"


void U2NITS::CInput::readConfig() {
	LogWriter::log(SystemInfo::getCurrentDateTime() + "\n");
	TomlFileManager::getInstance()->readFileAndParseFile(FilePathManager::getInstance()->getTomlFilePath());


	if (GlobalPara::basic::dimension != 2) {
		LogWriter::logAndPrintError("invalid dimension type.\n");
		exit(2);
	}
	if (GlobalPara::physicsModel::equation != _EQ_euler) {
		LogWriter::logAndPrintError("Invalid equation type.\n");
		exit(_EQ_euler);
	}
}

void U2NITS::CInput::readField_old() {

	//// 初始化
	// 尝试读取续算文件。若读取失败或者_continue==false则从头开始
	bool startFromZero = false;
	if (GlobalPara::basic::_continue) {
		int flag_readContinue = ContinueFileReader::readContinueFile_1();
		if (flag_readContinue == -1) {
			LogWriter::logAndPrint("Fail to read previous mesh(pause_*.dat). Will try to start from zero again.\n");
			startFromZero = true;
			// 防止后面histWriter不写文件头
			GlobalPara::basic::_continue = false;
		}
	}
	else {
		startFromZero = true;
	}
	// 从头开始。读取网格、初始化
	if (startFromZero) {
		const std::string& type = GlobalPara::basic::meshFileType;
		if (type == "inp") {
			std::string dir = FilePathManager::getInstance()->getInputDirectory();
			int flag_readMesh = InpMeshReader::readGmeshFile(dir + GlobalPara::basic::filename + ".inp");
			if (flag_readMesh == -1) {
				LogWriter::logAndPrintError("Fail to read mesh. Program will exit.\n");
				return;//退出
			}
			FieldInitializer::setInitialAndBoundaryCondition();
		}
		else if (type == "su2") {
			std::string dir = FilePathManager::getInstance()->getInputDirectory();
			int flag_readMesh = SU2MeshReader::readFile_2(dir + GlobalPara::basic::filename + ".su2", true);
			if (flag_readMesh == -1) {
				LogWriter::logAndPrint("Fail to read mesh. Program will exit.\n");
				return;//退出
			}
			FieldInitializer::setInitialAndBoundaryCondition();
		}
		else {
			LogWriter::logAndPrint("Invalid mesh file type: " + type + ". Program will exit.\n");
			return;
		}

	}

	// 日志记录边界参数
	using namespace GlobalPara::boundaryCondition::_2D;
	LogWriter::writeBoundaryCondition(inlet::ruvp, outlet::ruvp, inf::ruvp, 4);
}

void U2NITS::CInput::readField_2() {
	bool startFromZero = false;
	if (GlobalPara::basic::_continue) {
		int flag_readContinue = ContinueFileReader::readContinueFile_2();
		if (flag_readContinue != 0) {
			LogWriter::logAndPrint("Fail to read previous mesh(pause_*.dat).\n");
			startFromZero = true;
			GlobalPara::basic::_continue = false;// 防止后面histWriter不写文件头
		}
	}
	else {
		startFromZero = true;
	}
	// 从头开始。读取网格、初始化
	if (startFromZero) {
		const std::string& type = GlobalPara::basic::meshFileType;
		std::string dir = FilePathManager::getInstance()->getInputDirectory();
		int flag_readMesh = 502;
		if (type == "inp") {
			flag_readMesh = InpMeshReader::readGmeshFile_2(dir + GlobalPara::basic::filename + ".inp");
		}
		else if (type == "su2") {
			flag_readMesh = SU2MeshReader::readFile_2(dir + GlobalPara::basic::filename + ".su2", true);
		}
		else {
			LogWriter::logAndPrintError("Invalid mesh file type: " + type + "\n");
			flag_readMesh = 404;
		}
		if (flag_readMesh != 0) {
			LogWriter::logAndPrintError("Fail to read mesh. Program will exit.\n");
			exit(-3);//退出
		}
		FieldInitializer::setInitialAndBoundaryCondition();

	}
	// 日志记录边界参数
	using namespace GlobalPara::boundaryCondition::_2D;
	LogWriter::writeBoundaryCondition(inlet::ruvp, outlet::ruvp, inf::ruvp, 4);
}
