#include "FilePathManager.h"
#include <direct.h>
#include <iostream>
#include <vector>
#include <io.h>
#include "JsonConfig.h"
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
	m_jsonFilePath = m_inputDirectory + m_JsonFileName;

	// 设置标志符，表示已初始化
	if_initialized = true;
}


std::string FilePathManager::getTomlFilePath() {
	// 找不到默认input.toml，则用户输入
	if (_access(m_tomlFilePath.c_str(), 0) == -1) {// 检测文件/目录是否存在 成功返回0，失败返回-1
		std::cout << "Current working directory: " << m_workingDirectory << "\n";
		std::cout << "Cannot find toml file. Please input the file path: \n";
		std::cin >> m_tomlFilePath;
	}
	// 仍找不到，则加上.toml后缀再尝试
	if (_access(m_tomlFilePath.c_str(), 0) == -1) {
		m_tomlFilePath = m_tomlFilePath + ".toml";
	}
	// 仍找不到，则报错退出
	if (_access(m_tomlFilePath.c_str(), 0) == -1) {
		LogWriter::logAndPrintError("cannot find file\n");
		exit(-1);
	}

	return m_tomlFilePath; 
}

std::string FilePathManager::getJsonFilePath() {
	if (_access(m_jsonFilePath.c_str(), 0) == -1) {// 失败
		std::cout << "Current working directory: " << m_workingDirectory << "\n";
		std::cout << "Cannot find json file. Please input the file path: \n";
		std::cin >> m_jsonFilePath;
	}

	return m_jsonFilePath; 
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
	if (JsonConfig::config.FindMember("lastfilename") != JsonConfig::config.MemberEnd()) {
		lastf = JsonConfig::config["lastfilename"].GetString();//将config的值赋给lastf
	}
	if (1) {


		//输入提示
		std::cout << "Please input filename (without suffix):" << std::endl;
		if (lastf != "[NULL]")
			//若lastf!="[NULL]"则说明可以直接读取上次filename
			std::cout << "(Press a single \"Enter\" to use last filename [" << lastf << "])\n";
		else {
			//若lastf=="[NULL]"则说明是原本没有config文件，因此不能使用last filename
			//检测所有inp文件
		}
		//接受用户输入，作为文件名。若输入为空(即直接按enter)，则使用lastf
		char a = 0;
		std::string str;
		while (a != '\n') {
			a = getchar();
			str.push_back(a);
		}
		str.pop_back();//删除最后的\n
		if (str == "")str = lastf;
		GlobalPara::basic::filename = str;
	}


}

void FilePathManager::m_initializeWorkingDirectory() {
	// 获取工作目录
	char* buffer = _getcwd(nullptr, 0);
	if (buffer == nullptr) {
		std::cout << "Error: getWorkingDirectory() failed.\n";
		m_workingDirectory = "";
		return;
	}
	std::string workingDirectory(buffer);
	workingDirectory += "\\";
	//std::cout << "Working Directory: " << workingDirectory << std::endl;
	free(buffer);// free不可删 详见 https://baike.baidu.com/item/getcwd/4746955?fr=ge_ala
	m_workingDirectory = workingDirectory;
}

void FilePathManager::createFolderIfDoesntExist(std::string foldername) {
	if(foldername == "") return;
	// 工作目录下检测文件/目录是否存在
	if (_access(foldername.c_str(), 0) == 0) {// 成功返回0，失败返回-1

	}
	else {
		system(("mkdir " + m_workingDirectory + foldername).c_str());
	}
}
