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

	createFolderIfDoesntExist(m_inputFolderName);
	createFolderIfDoesntExist(m_outputFolderName);

	m_initializeWorkingDirectory();
	m_inputDirectory = m_workingDirectory + m_inputFolderName + "\\";
	m_outputDirectory = m_workingDirectory + m_outputFolderName + "\\";
	
	m_tomlFilePath = m_inputDirectory + m_TomlFileName;
	m_jsonFilePath = m_inputDirectory + m_JsonFileName;

	// ���ñ�־������ʾ�ѳ�ʼ��
	if_initialized = true;
}


std::vector<std::string> FilePathManager::getFiles(std::string path) {
	path += "*";
	_finddata64i32_t fileInfo;
	std::vector<std::string> files;
	intptr_t hFile = _findfirst(path.c_str(), &fileInfo);
	if (hFile == -1) {
		std::cout << "Error in FilePathManager::getFiles()" << std::endl;
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

void FilePathManager::m_initializeWorkingDirectory() {
	// ��ȡ����Ŀ¼
	char* buffer = _getcwd(nullptr, 0);
	if (buffer == nullptr) {
		std::cout << "Error: getWorkingDirectory() failed.\n";
		m_workingDirectory = "";
		return;
	}
	std::string workingDirectory(buffer);
	workingDirectory += "\\";
	//std::cout << "Working Directory: " << workingDirectory << std::endl;
	free(buffer);// free����ɾ ��� https://baike.baidu.com/item/getcwd/4746955?fr=ge_ala
	m_workingDirectory = workingDirectory;
}

void FilePathManager::createFolderIfDoesntExist(std::string foldername) {
	if(foldername == "") return;
	// ����Ŀ¼�¼���ļ�/Ŀ¼�Ƿ����
	if (_access(foldername.c_str(), 0) == 0) {// �ɹ�����0��ʧ�ܷ���-1

	}
	else {
		system(("mkdir " + m_workingDirectory + foldername).c_str());
	}
}
