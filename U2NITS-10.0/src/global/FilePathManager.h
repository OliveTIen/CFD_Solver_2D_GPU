
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
	// 构造函数，在第一次获取实例时调用
	FilePathManager() { initialzeVariables(); };
	static FilePathManager* class_pointer;
// -- 单例模式 end

private:
	bool if_initialized = false;

public:
	std::string inputFolder_name = "input";// = "input"
	std::string outputFolder_name = "output";// = "output"
	std::string workingDirectory;// 当前执行文件所在的目录
	std::string exePath_withSlash;// 有斜杠
	std::string inputTomlFile_path;
	std::string inputJsonFile_path;
	std::string outputFolder_path;// 输出文件目录，无斜杠

public:
	void createFolderIfDoesntExist(std::string foldername);
	// 初始化，在构造函数中调用
	void initialzeVariables();
	std::string getExePath(int if_contain_slash = 1);
	// 此处flag无任何作用，只是为了重构代码时比较好过渡
	std::string getExePath_withSlash(int flag = 1);
	// 列举文件名称，path需加反斜杠[todo]该函数仅仅列出了目录下文件，没有列出子目录
	static std::vector<std::string> ls(std::string path);

	static void getFileName_UseUserInput();


private:
};

#endif