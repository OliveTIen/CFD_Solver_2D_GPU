
#ifndef FILE_PATH_MANAGER_H
#define FILE_PATH_MANAGER_H

#include <string>
#include <vector>
class FilePathManager {
// -- 单例模式 start
// 单例模式使得每次得到的是同一个实例。它的作用类似静态成员变量
public:
	static FilePathManager* getInstance();

private:
	FilePathManager() {};
	static FilePathManager* class_pointer;
// -- 单例模式 end

private:
	bool if_initialized = false;
	std::string workingDirectory;// 当前执行文件所在的目录
	std::string exePath_withSlash;//有斜杠
	std::string residualFileName;

public:
	void initialzeVariables();
	std::string getExePath(int if_contain_slash = 1);
	// 此处flag无任何作用，只是为了重构代码时比较好过渡
	std::string getExePath_withSlash(int flag = 1);

private:
	std::vector<std::string> ls(std::string path);
	void createFolderIfDoesntExist(std::string foldername);
};

#endif