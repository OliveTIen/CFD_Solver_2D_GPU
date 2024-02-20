#include "FVM.h"
#include "FVM_2D.h"
#include "JsonFileManager.h"
#include "output/ConsolePrinter.h"
#include "global/FilePathManager.h"
#include "global/SystemInfo.h"
#include "global/Config.h"
#include "output/LogWriter.h"
#include "input/TomlFileManager.h"

void FVM::run() {
	ConsolePrinter::printHeader();

	// 初始化路径，创建文件夹
	FilePathManager* filePathManager = FilePathManager::getInstance();
	filePathManager->createFolderIfDoesntExist("input");
	filePathManager->createFolderIfDoesntExist("output");
	// 生成日志
	LogWriter::writeLog(SystemInfo::getCurrentDateTime() + "\n", 0);
	// 读取input.toml输入参数，处理部分输入参数
	TomlFileManager tomlFileManager;
	tomlFileManager.readFileAndParseFile(filePathManager->getExePath_withSlash() + "input\\input.toml");
	//tomlFileManager.printParsedFile();
	// 读网格，初始化初值，求解
	if (GlobalPara::basic::dimension == 1) {
		std::cout << "Error: invalid dimension type. 1D has been deleted.\n";
	}
	else if (GlobalPara::basic::dimension == 2) {
		FVM_2D fvm2d;
		if (GlobalPara::basic::useGPU) {
			fvm2d.run_GPU();
		}
		else {
			fvm2d.run();
		}
	}
}
