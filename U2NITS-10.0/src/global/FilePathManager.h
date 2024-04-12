
#ifndef FILE_PATH_MANAGER_H
#define FILE_PATH_MANAGER_H

#include <string>
#include <vector>
class FilePathManager {
// ����ģʽʹ��ÿ�εõ�����ͬһ��ʵ���������������ƾ�̬��Ա����
private:
	static FilePathManager* class_pointer;
	bool if_initialized = false;
	std::string m_inputFolderName = "";
	std::string m_outputFolderName = "";
	std::string m_TomlFileName = "input.toml";
	std::string m_JsonFileName = "input.json";

	std::string m_workingDirectory;// ��ǰִ���ļ����ڵ�Ŀ¼ ��б��
	std::string m_inputDirectory;// �����ļ�Ŀ¼
	std::string m_outputDirectory;// ����ļ�Ŀ¼

	std::string m_tomlFilePath;
	std::string m_jsonFilePath;

private:
	FilePathManager() { initialzeVariables(); };// ���캯������Ϊ˽�У��ڵ�һ�λ�ȡʵ��ʱ����
	void initialzeVariables();
	void m_initializeWorkingDirectory();
	void createFolderIfDoesntExist(std::string foldername);


public:
	static FilePathManager* getInstance();// �ӿڣ���ȡ����ָ��

	std::string getWorkingDirectory(){ return m_workingDirectory; }
	std::string getInputDirectory() { return m_inputDirectory; }
	std::string getOutputDirectory() { return m_outputDirectory; }
	std::string getTomlFilePath();
	std::string getJsonFilePath();
	// �о��ļ����ƣ�path��ӷ�б��[todo]�ú��������г���Ŀ¼���ļ���û���г���Ŀ¼
	std::vector<std::string> getFiles(std::string path);
	void getFileName_UseUserInput();

};

#endif