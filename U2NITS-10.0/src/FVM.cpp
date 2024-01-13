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
	// ��ȡ·��
	GlobalStatic::iniExePath();
	// �����ļ���
	jsonFileManager.createFolderIfDoesntExist("input");
	jsonFileManager.createFolderIfDoesntExist("output");
	jsonFileManager.createFolderIfDoesntExist(".config");
	// ��ȡ�����ļ�����ȡ�ϴ��ļ���
	Config::readConfig();
	// �û�������Ŀ���ơ�������global::readConfig()֮����Ϊ�õ���lastFileName
	GlobalStatic::getFileName_UseUserInput();
	Config::writeConfig(GlobalStatic::filename);
	// ������־�ļ�������filename��ʼ����ʹ��
	LogWriter::writeLog(SystemInfo::getCurrentDateTime() + "\n", 0);
	// ��ȡ��ʼ��������ʼ��GlobalPara
	jsonFileManager.readInputPara_initGlobalPara(GlobalStatic::exePath_withSlash + "input\\input.json");
	// ��δ����ruvp�������Ma��AOA����ruvp
	GlobalPara::boundaryCondition::_2D::ini_ruvp_by_Ma_AOA();
	// �����񣬳�ʼ����ֵ�����
	if (GlobalPara::basic::dimension == 1) {
		LogWriter::writeLogAndCout("1D has been deleted.\n");
	}
	else if (GlobalPara::basic::dimension == 2) {
		FVM_2D fvm2d;
		fvm2d.run();
	}
	// ��β������дһЩ�ϴα�������á���ǰ���࣬���Ժ�һ������
	Config::writeConfig(GlobalStatic::filename);
}
