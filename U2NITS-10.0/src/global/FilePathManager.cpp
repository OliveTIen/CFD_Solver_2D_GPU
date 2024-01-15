#include "FilePathManager.h"
#include <direct.h>
#include <iostream>
#include <vector>
#include <io.h>
#include "Config.h"
#include "../GlobalPara.h"

FilePathManager* FilePathManager::getInstance() {
    if (!class_pointer) {
        class_pointer = new FilePathManager();
    }
    return class_pointer;
}

FilePathManager* FilePathManager::class_pointer = nullptr;

void FilePathManager::initialzeVariables() {

	// ��ȡ����Ŀ¼����·�� ���������"/"
	const int maxPathLength = 260;
	char buffer[maxPathLength];
	_getcwd(buffer, maxPathLength);
	workingDirectory = buffer;
	// ����Ŀ¼
	exePath_withSlash = workingDirectory + "\\";
	inputTomlFile_path = workingDirectory + "\\" + inputFolder_name + "\\" + "input.toml";
	inputJsonFile_path = workingDirectory + "\\" + inputFolder_name + "\\" + "input.json";
	outputFolder_path = workingDirectory + "\\" + outputFolder_name;

	// ���ñ�־������ʾ�ѳ�ʼ��
	if_initialized = true;
}

std::string FilePathManager::getExePath(int if_contain_slash) {
	if (!if_initialized) {
		std::cout << "Error: FilePathManager not initialized.\n";
	}
	std::string str = workingDirectory;
	if (if_contain_slash == 1) {
		str += std::string("\\");
	}
	return str;
}

std::string FilePathManager::getExePath_withSlash(int flag) {
	if (!if_initialized) {
		std::cout << "Error: FilePathManager not initialized.\n";
	}
	return exePath_withSlash;
}

std::vector<std::string> FilePathManager::ls(std::string path) {
	path += "*";
	_finddata64i32_t fileInfo;
	std::vector<std::string> files;
	intptr_t hFile = _findfirst(path.c_str(), &fileInfo);
	if (hFile == -1) {
		std::cout << "Error in FilePathManager::ls()" << std::endl;
		return files;
	}
	do {
		files.push_back(fileInfo.name);
	} while (_findnext(hFile, &fileInfo) == 0);
	return files;
}

void FilePathManager::getFileName_UseUserInput() {
	std::string lastf;
	if (Config::config.FindMember("lastfilename") != Config::config.MemberEnd()) {
		lastf = Config::config["lastfilename"].GetString();//��config��ֵ����lastf
	}
	if (1) {


		//������ʾ
		std::cout << "Please input filename (without suffix):" << std::endl;
		if (lastf != "[NULL]")
			//��lastf!="[NULL]"��˵������ֱ�Ӷ�ȡ�ϴ�filename
			std::cout << "(Press a single \"Enter\" to use last filename [" << lastf << "])\n";
		else {
			//��lastf=="[NULL]"��˵����ԭ��û��config�ļ�����˲���ʹ��last filename
			//�������inp�ļ�
		}
		//�����û����룬��Ϊ�ļ�����������Ϊ��(��ֱ�Ӱ�enter)����ʹ��lastf
		char a = 0;
		std::string str;
		while (a != '\n') {
			a = getchar();
			str.push_back(a);
		}
		str.pop_back();//ɾ������\n
		if (str == "")str = lastf;
		GlobalPara::basic::filename = str;
	}


}

void FilePathManager::createFolderIfDoesntExist(std::string foldername) {
	// ����Ŀ¼�¼���ļ�/Ŀ¼�Ƿ����
	if (_access(foldername.c_str(), 0) == 0) {// �ɹ�����0��ʧ�ܷ���-1
		// ����
		//std::cout << foldername << " exists.\n";
	}
	else {
		// ������
		//std::cout << foldername << " doesn't exist.\n";
		system(("mkdir " + workingDirectory + "\\" + foldername).c_str());
	}
}
