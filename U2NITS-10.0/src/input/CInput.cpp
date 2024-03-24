#include "CInput.h"
#include "Initializer.h"
#include "InpMeshReader.h"
#include "SU2MeshReader.h"
#include "TomlFileManager.h"
#include "../FVM_2D.h"
#include "../global/FilePathManager.h"
#include "../global/GlobalPara.h"
#include "../global/SystemInfo.h"
#include "../output/LogWriter.h"
#include "../output/ConsolePrinter.h"


void U2NITS::CInput::readConfig() {
	LogWriter::writeLogAndCout("Using VS Profiler. Please remove '/Profile' option on release.\n");
	ConsolePrinter::printHeader();
	LogWriter::writeLog(SystemInfo::getCurrentDateTime() + "\n");
	TomlFileManager::getInstance()->readFileAndParseFile(FilePathManager::getInstance()->getTomlFilePath());


	if (GlobalPara::basic::dimension != 2) {
		LogWriter::writeLogAndCout("Error: invalid dimension type.\n", LogWriter::Error, LogWriter::Error);
		exit(2);
	}
	if (GlobalPara::physicsModel::equation != _EQ_euler) {
		LogWriter::writeLogAndCout("Error: Invalid equation type.\n", LogWriter::Error, LogWriter::Error);
		exit(_EQ_euler);
	}
}

void U2NITS::CInput::readField() {

	FVM_2D* pFVM2D = FVM_2D::getInstance();
	//// 初始化
	// 尝试读取续算文件。若读取失败或者_continue==false则从头开始
	bool startFromZero = false;
	if (GlobalPara::basic::_continue) {
		int flag_readContinue = pFVM2D->readContinueFile();
		if (flag_readContinue == -1) {
			std::string error_msg = "Info: Fail to read previous mesh(pause_*.dat). ";
			error_msg += "Will try to start from zero again.\n";
			LogWriter::writeLogAndCout(error_msg, LogWriter::Info);
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
				std::string error_msg = "Error: Fail to read mesh. Program will exit.\n";
				LogWriter::writeLogAndCout(error_msg);
				return;//退出
			}
			Initializer::setInitialAndBoundaryCondition();
		}
		else if (type == "su2") {
			std::string dir = FilePathManager::getInstance()->getInputDirectory();
			int flag_readMesh = SU2MeshReader::readFile(dir + GlobalPara::basic::filename + ".su2", true);
			if (flag_readMesh == -1) {
				std::string error_msg = "Error: Fail to read mesh. Program will exit.\n";
				LogWriter::writeLogAndCout(error_msg);
				return;//退出
			}
			Initializer::setInitialAndBoundaryCondition();
		}
		else {
			std::string error_msg = "Invalid mesh file type: " + type + ". Program will exit.\n";
			LogWriter::writeLogAndCout(error_msg);
			return;
		}

	}

	// 日志记录边界参数
	using namespace GlobalPara::boundaryCondition::_2D;
	LogWriter::writeBoundaryCondition(inlet::ruvp, outlet::ruvp, inf::ruvp, 4);
}
