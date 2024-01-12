#ifndef FILEMANAGER_H
#define FILEMANAGER_H
#include "head.h"

class FileManager {
public:
	rapidjson::Document inputpara;

public:
	//json树赋默认值，输出为文件
	void useDefaultInputPara_and_writeInputPara();
	//将json树输出为文件。须在global::exePath初始化、已存在folder后使用
	void writeInputPara(std::string filepath_name);
	//读取初始参数，初始化GlobalPara
	void readInputPara_initGlobalPara(std::string filepath_name);
	//[未完成]检查资源完整性。目的是防止json文件被乱改
	void checkIntegrity();
	//json树赋给glboal:: data
	void useJsonTreeToUpdateGlobalPara();//赋给globalPara
	//若obj存在该成员，则将对应值赋给variable，否则抛出异常。若不支持类型，也抛出异常
	template <typename T_type>
	void getMemberValue(rapidjson::Value& obj, const char* membername, T_type& variable);
	
	//文件操作系列
	//创建文件夹。若已存在，则不创建。须在global::exePath初始化后使用
	void createFolder(std::string foldername);
};

#endif

template<typename T_type>
inline void FileManager::getMemberValue(rapidjson::Value& obj, const char* membername, T_type& variable) {
	//更新内容后，注意在模板中添加所有现有的实例化，否则会强制类型转换
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
