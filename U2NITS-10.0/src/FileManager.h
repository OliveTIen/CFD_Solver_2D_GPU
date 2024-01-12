#ifndef FILEMANAGER_H
#define FILEMANAGER_H
#include "head.h"

class FileManager {
public:
	rapidjson::Document inputpara;

public:
	//json����Ĭ��ֵ�����Ϊ�ļ�
	void useDefaultInputPara_and_writeInputPara();
	//��json�����Ϊ�ļ�������global::exePath��ʼ�����Ѵ���folder��ʹ��
	void writeInputPara(std::string filepath_name);
	//��ȡ��ʼ��������ʼ��GlobalPara
	void readInputPara_initGlobalPara(std::string filepath_name);
	//[δ���]�����Դ�����ԡ�Ŀ���Ƿ�ֹjson�ļ����Ҹ�
	void checkIntegrity();
	//json������glboal:: data
	void useJsonTreeToUpdateGlobalPara();//����globalPara
	//��obj���ڸó�Ա���򽫶�Ӧֵ����variable�������׳��쳣������֧�����ͣ�Ҳ�׳��쳣
	template <typename T_type>
	void getMemberValue(rapidjson::Value& obj, const char* membername, T_type& variable);
	
	//�ļ�����ϵ��
	//�����ļ��С����Ѵ��ڣ��򲻴���������global::exePath��ʼ����ʹ��
	void createFolder(std::string foldername);
};

#endif

template<typename T_type>
inline void FileManager::getMemberValue(rapidjson::Value& obj, const char* membername, T_type& variable) {
	//�������ݺ�ע����ģ��������������е�ʵ�����������ǿ������ת��
	const char* err = "Error: Unsupported type, in function FileManager::getChildValue\n";
	const char* err2 = "Error: Unexisted key, in function FileManager::getChildValue\n";
	if (obj.HasMember(membername)) {
		if (typeid(variable) == typeid(1))variable = obj[membername].GetInt();
		else if (typeid(variable) == typeid(1.0))variable = obj[membername].GetDouble();
		else if (typeid(variable) == typeid(true))variable = obj[membername].GetInt();
		else throw err;
	}
	else throw err2;
}
