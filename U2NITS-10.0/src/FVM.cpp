#include "FVM.h"
#include "FVM_2D.h"
#include "JsonFileManager.h"
#include "output/ConsolePrinter.h"
#include "global/FilePathManager.h"
#include "global/SystemInfo.h"
#include "global/Config.h"
#include "output/LogWriter.h"

void FVM::run() {
	FilePathManager* filePathManager = FilePathManager::getInstance();
	filePathManager->initialzeVariables();
	JsonFileManager jsonFileManager;
	ConsolePrinter::printHeader();
	// 获取路径
	GlobalStatic::iniExePath();
	// 创建文件夹
	jsonFileManager.createFolderIfDoesntExist("input");
	jsonFileManager.createFolderIfDoesntExist("output");
	jsonFileManager.createFolderIfDoesntExist(".config");
	// 读取配置文件：读取上次文件名
	Config::readConfig();
	// 用户输入项目名称。必须在global::readConfig()之后，因为用到了lastFileName
	GlobalStatic::getFileName_UseUserInput();
	Config::writeConfig(GlobalStatic::filename);
	// 生成日志文件。须在filename初始化后使用
	LogWriter::writeLog(SystemInfo::getCurrentDateTime() + "\n", 0);
	// 读取初始参数，初始化GlobalPara
	jsonFileManager.readInputPara_initGlobalPara(GlobalStatic::exePath_withSlash + "input\\input.json");
	// 若未给定ruvp，则根据Ma和AOA计算ruvp
	GlobalPara::boundaryCondition::_2D::ini_ruvp_by_Ma_AOA();
	// 读网格，初始化初值，求解
	if (GlobalPara::basic::dimension == 1) {
		LogWriter::writeLogAndCout("1D has been deleted.\n");
	}
	else if (GlobalPara::basic::dimension == 2) {
		FVM_2D fvm2d;
		fvm2d.run();
	}
	// 收尾工作，写一些上次保存的配置。当前多余，但以后不一定多余
	Config::writeConfig(GlobalStatic::filename);
}
