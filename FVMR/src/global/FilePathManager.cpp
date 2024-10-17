#include "FilePathManager.h"
#include <direct.h>
#include <iostream>
#include <vector>
#include <io.h>
#include "../global/GlobalPara.h"
#include "../output/LogWriter.h"

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
	if (m_inputFolderName != "") {
		m_inputDirectory = m_workingDirectory + m_inputFolderName + "\\";
	}
	else {
		m_inputDirectory = m_workingDirectory;
	}
	if (m_outputFolderName != "") {
		m_outputDirectory = m_workingDirectory + m_outputFolderName + "\\";
	}
	else {
		m_outputDirectory = m_workingDirectory;
	}

	m_tomlFilePath = m_inputDirectory + m_TomlFileName;

	// ���ñ�־������ʾ�ѳ�ʼ��
	if_initialized = true;
}


std::string FilePathManager::getTomlFilePath() {
	// �Ҳ���Ĭ��input.toml�����û�����
	if (_access(m_tomlFilePath.c_str(), 0) == -1) {// ����ļ�/Ŀ¼�Ƿ���� �ɹ�����0��ʧ�ܷ���-1
		std::cout << "Current working directory: " << m_workingDirectory << "\n";
		std::cout << "Cannot find toml file. Please input the file path: \n";
		std::cin >> m_tomlFilePath;
	}
	// ���Ҳ����������.toml��׺�ٳ���
	if (_access(m_tomlFilePath.c_str(), 0) == -1) {
		m_tomlFilePath = m_tomlFilePath + ".toml";
	}
	// ���Ҳ������򱨴��˳�
	if (_access(m_tomlFilePath.c_str(), 0) == -1) {
		LogWriter::logAndPrintError("cannot find file\n");
		exit(-1);
	}

	return m_tomlFilePath; 
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
