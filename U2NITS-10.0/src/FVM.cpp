#include "FVM.h"
#include "FVM_2D.h"
#include "FileManager.h"

void FVM::run() {
	FileManager myFileManager;
	cmd::printHeader();
	// ��ȡ·��
	global::iniExePath();
	// �����ļ���
	myFileManager.createFolderIfDoesntExist("input");
	myFileManager.createFolderIfDoesntExist("output");
	myFileManager.createFolderIfDoesntExist(".config");
	// ��ȡ�����ļ�����ȡ�ϴ��ļ���
	global::readConfig();
	// �û�������Ŀ���ơ�������global::readConfig()֮����Ϊ�õ���lastFileName
	global::getFileName_UseUserInput();
	// ������־�ļ�������filename��ʼ����ʹ��
	global::writeLog(global::currentDateTime() + "\n", 0);
	// ��ȡ��ʼ��������ʼ��GlobalPara
	myFileManager.readInputPara_initGlobalPara(global::exePath + "input\\input.json");
	// ��δ����ruvp�������Ma��AOA����ruvp
	GlobalPara::boundaryCondition::_2D::ini_ruvp_by_Ma_AOA();
	// �����񣬳�ʼ����ֵ�����
	if (GlobalPara::basic::dimension == 1) {
		global::writeLogAndCout("1D has been deleted.\n");
	}
	else if (GlobalPara::basic::dimension == 2) {
		FVM_2D fvm2d;
		fvm2d.run();
	}
	// ��β������дһЩ�ϴα�������á���ǰ���࣬���Ժ�һ������
	global::writeConfig();
}
