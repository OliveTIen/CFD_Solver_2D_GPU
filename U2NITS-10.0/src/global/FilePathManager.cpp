#include "FilePathManager.h"
#include <direct.h>
#include <iostream>
#include <vector>
#include <io.h>

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
	//
	exePath_withSlash = workingDirectory + "\\";

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

void FilePathManager::createFolderIfDoesntExist(std::string foldername) {
	std::string path = workingDirectory;
	std::vector<std::string> files = ls(path);
	for (int i = 0; i < files.size(); i++) {
		if (files[i] == foldername)return;//文件夹已经存在
	}
	system(("mkdir " + path + "\\" + foldername).c_str());
}
