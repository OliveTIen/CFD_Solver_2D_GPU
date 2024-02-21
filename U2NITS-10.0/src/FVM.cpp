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

	// ��ʼ��·���������ļ���
	FilePathManager* filePathManager = FilePathManager::getInstance();
	filePathManager->createFolderIfDoesntExist("input");
	filePathManager->createFolderIfDoesntExist("output");
	// ������־
	LogWriter::writeLog(SystemInfo::getCurrentDateTime() + "\n", 0);
	// ��ȡinput.toml����������������������
	TomlFileManager tomlFileManager;
	tomlFileManager.readFileAndParseFile(filePathManager->getExePath_withSlash() + "input\\input.toml");
	//tomlFileManager.printParsedFile();
	// �����񣬳�ʼ����ֵ�����
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
