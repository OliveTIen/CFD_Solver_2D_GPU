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

	// 获取工作目录所在路径 不含后面的"/"
	const int maxPathLength = 260;
	char buffer[maxPathLength];
	_getcwd(buffer, maxPathLength);
	workingDirectory = buffer;
	// 其他目录
	exePath_withSlash = workingDirectory + "\\";
	inputTomlFile_path = workingDirectory + "\\" + inputFolder_name + "\\" + "input.toml";
	inputJsonFile_path = workingDirectory + "\\" + inputFolder_name + "\\" + "input.json";
	outputFolder_path = workingDirectory + "\\" + outputFolder_name;

	// 设置标志符，表示已初始化
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
		lastf = Config::config["lastfilename"].GetString();//将config的值赋给lastf
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

void FilePathManager::createFolderIfDoesntExist(std::string foldername) {
	// 工作目录下检测文件/目录是否存在
	if (_access(foldername.c_str(), 0) == 0) {// 成功返回0，失败返回-1
		// 存在
		//std::cout << foldername << " exists.\n";
	}
	else {
		// 不存在
		//std::cout << foldername << " doesn't exist.\n";
		system(("mkdir " + workingDirectory + "\\" + foldername).c_str());
	}
}
