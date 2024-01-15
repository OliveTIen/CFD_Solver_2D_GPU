
#ifndef FILE_PATH_MANAGER_H
#define FILE_PATH_MANAGER_H

#include <string>
#include <vector>
class FilePathManager {
// -- ����ģʽ start
// ����ģʽʹ��ÿ�εõ�����ͬһ��ʵ���������������ƾ�̬��Ա����
public:
	static FilePathManager* getInstance();

private:
	// ���캯�����ڵ�һ�λ�ȡʵ��ʱ����
	FilePathManager() { initialzeVariables(); };
	static FilePathManager* class_pointer;
// -- ����ģʽ end

private:
	bool if_initialized = false;

public:
	std::string inputFolder_name = "input";// = "input"
	std::string outputFolder_name = "output";// = "output"
	std::string workingDirectory;// ��ǰִ���ļ����ڵ�Ŀ¼
	std::string exePath_withSlash;// ��б��
	std::string inputTomlFile_path;
	std::string inputJsonFile_path;
	std::string outputFolder_path;// ����ļ�Ŀ¼����б��

public:
	void createFolderIfDoesntExist(std::string foldername);
	// ��ʼ�����ڹ��캯���е���
	void initialzeVariables();
	std::string getExePath(int if_contain_slash = 1);
	// �˴�flag���κ����ã�ֻ��Ϊ���ع�����ʱ�ȽϺù���
	std::string getExePath_withSlash(int flag = 1);
	// �о��ļ����ƣ�path��ӷ�б��[todo]�ú��������г���Ŀ¼���ļ���û���г���Ŀ¼
	static std::vector<std::string> ls(std::string path);

	static void getFileName_UseUserInput();


private:
};

#endif