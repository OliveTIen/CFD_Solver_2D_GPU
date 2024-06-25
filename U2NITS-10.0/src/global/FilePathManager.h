
#ifndef FILE_PATH_MANAGER_H
#define FILE_PATH_MANAGER_H

#include <string>
#include <vector>
class FilePathManager {
// 单例模式使得每次得到的是同一个实例。它的作用类似静态成员变量
private:
	static FilePathManager* class_pointer;
	bool if_initialized = false;
	std::string m_inputFolderName = "";
	std::string m_outputFolderName = "";
	std::string m_TomlFileName = "input.toml";

	std::string m_workingDirectory;// 当前执行文件所在的目录 有斜杠
	std::string m_inputDirectory;// 输入文件目录
	std::string m_outputDirectory;// 输出文件目录

	std::string m_tomlFilePath;

private:
	FilePathManager() { initialzeVariables(); };// 构造函数，设为私有，在第一次获取实例时调用
	void initialzeVariables();
	void m_initializeWorkingDirectory();
	void createFolderIfDoesntExist(std::string foldername);


public:
	static FilePathManager* getInstance();// 接口，获取单例指针

	std::string getWorkingDirectory(){ return m_workingDirectory; }
	std::string getInputDirectory() { return m_inputDirectory; }
	std::string getOutputDirectory() { return m_outputDirectory; }
	std::string getTomlFilePath();
	// 列举文件名称，path需加反斜杠[todo]该函数仅仅列出了目录下文件，没有列出子目录
	std::vector<std::string> getFiles(std::string path);

};

#endif