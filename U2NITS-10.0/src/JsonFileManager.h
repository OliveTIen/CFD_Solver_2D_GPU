#ifndef JSONFILEMANAGER_H
#define JSONFILEMANAGER_H
#include "head.h"

class JsonFileManager {
public:
	rapidjson::Document inputpara;

public:
	//json����Ĭ��ֵ�����Ϊ�ļ�
	void useDefaultInputPara_and_writeInputPara();
	//(����һ����������)��json�����Ϊ�ļ�������global::exePath��ʼ�����Ѵ���folder��ʹ��
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
	void getMemberString(rapidjson::Value& obj, const char* membername, std::string& variable);
	
	//�ļ�����ϵ��
	//�����ļ��С����Ѵ��ڣ��򲻴���������global::exePath��ʼ����ʹ��
	void createFolderIfDoesntExist(std::string foldername);
};

#endif

template<typename T_type>
inline void JsonFileManager::getMemberValue(rapidjson::Value& obj, const char* membername, T_type& variable) {
	//�������ݺ�ע����ģ���������������е�ʵ�����������ǿ������ת��
	const char* err = "Error: Unsupported type, in function JsonFileManager::getMemberValue\n";
	const char* err2 = "Error: Unexisted key, in function JsonFileManager::getMemberValue\n";
	if (obj.HasMember(membername)) {
		if (typeid(variable) == typeid(1))variable = obj[membername].GetInt();
		else if (typeid(variable) == typeid(1.0))variable = obj[membername].GetDouble();
		else if (typeid(variable) == typeid(true))variable = obj[membername].GetInt();
		else throw err;

	}
	else {
		std::cout << "Error in " << membername << " ";
		//throw err2;
	}
		
}

inline void JsonFileManager::getMemberString(rapidjson::Value& obj, const char* membername, std::string& variable) {
	const char* err2 = "Error: Unexisted key, in function JsonFileManager::getMemberString\n";
	if (obj.HasMember(membername)) {
		variable = obj[membername].GetString();
	}
	else throw err2;

}