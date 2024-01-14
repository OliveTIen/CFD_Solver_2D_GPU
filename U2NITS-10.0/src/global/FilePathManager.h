
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
	FilePathManager() {};
	static FilePathManager* class_pointer;
// -- ����ģʽ end

private:
	bool if_initialized = false;
	std::string workingDirectory;// ��ǰִ���ļ����ڵ�Ŀ¼
	std::string exePath_withSlash;//��б��
	std::string residualFileName;

public:
	void initialzeVariables();
	std::string getExePath(int if_contain_slash = 1);
	// �˴�flag���κ����ã�ֻ��Ϊ���ع�����ʱ�ȽϺù���
	std::string getExePath_withSlash(int flag = 1);

private:
	std::vector<std::string> ls(std::string path);
	void createFolderIfDoesntExist(std::string foldername);
};

#endif