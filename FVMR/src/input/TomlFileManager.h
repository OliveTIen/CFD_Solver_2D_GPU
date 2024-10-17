#ifndef TOML_FILE_MANAGER_H
#define TOML_FILE_MANAGER_H
#include <string>
#include <iostream>
#include <sstream>
//#include "../include/toml/cpptoml.h"
#include <toml/cpptoml.h>
#include "../global/SGlobalPara.h"
#include "../gpu/datatype/DefineType.h"

class TomlFileManager {
private:
	// 成员变量
	static TomlFileManager* classPointer;
	std::shared_ptr<cpptoml::table>m_tree = nullptr;
	bool has_getValueFailed = false;

public:
	static TomlFileManager* getInstance();
	void readTomlFile(std::string fullFilePath);
	
	// 接口，遇到错误会报错退出
	template <typename T_type>
	void getValueOrExit(std::string key, T_type& value_to_be_changed);// 必选

	// 以下函数，若遇到错误不会立刻退出，只是将has_getValueFailed置为true
	template <typename T_type>
	void getValueOnCondition(std::string key, T_type& value_to_be_changed, bool condition);// condition==true，则必选，否则可选
	template <typename T_type>
	void getValue(std::string key, T_type& value_to_be_changed);// 必选
	template <typename T_type>
	void getValueIfExists(std::string key, T_type& value_to_be_changed);// 可选
	

private:
	TomlFileManager() {};
	// m_tree转化为全局变量
	void treeToGlobalParameter();
	// [todo]
	void treeToSGlobalPara(SGlobalPara::SGlobalPara& spara);
	// 统一相互矛盾的输入
	void handleConflictingInputs();
	void getValue_boundaryCondition2D(std::string parent, myfloat ruvp[4]);

	// tree是否包含key。输入参数以"."分隔层级
	bool treeContainsKey(std::string key) { return m_tree->contains_qualified(key); }
	// 打印tree是否包含key
	void printTreeContainsKeyOrNot(std::string key);
	// 打印tree
	void printTree() { std::cout << *m_tree << std::endl; }
	// 检测has_getValueFailed是否为true并退出
	void ifFailedThenExit();
	
	template <typename T_type>
	void getValueOrigin(std::string key, T_type& value_to_be_changed);
};

template<typename T_type>
inline void TomlFileManager::getValueOrigin(std::string key, T_type& value_to_be_changed) {

	bool failed = false;
	// 如果是bool类型
	if (std::is_same<T_type, bool>::value) {
		// 尝试读取字符串'bool'（不带引号）
		cpptoml::option<bool> valueBool = m_tree->get_qualified_as<bool>(key);
		if (valueBool) {
			value_to_be_changed = *valueBool;
		}
		// 如果无法读取字符串'bool'，则按int类型读取数字
		else {
			cpptoml::option<int> valueInt = m_tree->get_qualified_as<int>(key);
			if (valueInt) {
				value_to_be_changed = (bool)*valueInt;// int转换为bool。即使在tomlTree中存储的valueInt是int类型的2或-1，但由于value_to_be_changed是bool类型，因此隐式转换为1
			}
			else {
				failed = true;
			}
		}
	}
	// cpptoml不识别float格式，应按照double读取
	else if (std::is_same<T_type, float>::value) {
		cpptoml::option<double> valueDouble = m_tree->get_qualified_as<double>(key);
		if (valueDouble) {
			value_to_be_changed = *valueDouble;
		}
		else {
			failed = true;
		}
	}
	// 如果是其他类型，正常读取
	else {
		cpptoml::option<T_type> valueOther = m_tree->get_qualified_as<T_type>(key);
		if (valueOther) {
			value_to_be_changed = *valueOther;
		}
		else {
			failed = true;
		}
	}

	if (failed) {
		std::stringstream ss;
		//ss << "invalid_value_exception @TomlFileManager::getValueOrigin\n";
		ss << "Parameter \"" << key << "\" gets a wrong value type.\n";
		LogWriter::logAndPrintError(ss.str());


		has_getValueFailed = true;
	}
}

template<typename T_type>
inline void TomlFileManager::getValueOnCondition(std::string key, T_type& value_to_be_changed, bool condition) {
	if (treeContainsKey(key)) {
		getValueOrigin<T_type>(key, value_to_be_changed);
	}
	else {
		if (condition) {// 如果有条件，则必选，因此会报错
			LogWriter::logAndPrintError("Required parameter \"" + key + "\" doesn't exist.\n");
			has_getValueFailed = true;
		}
		else {// 无条件，则可选，因此不会报错
			//LogWriter::log("Parameter \"" + key + "\" uses default value.\n", LogWriter::Info);

			std::stringstream ss;
			ss << "parameter \"" << key << "\" uses default value: " << value_to_be_changed << "\n";
			LogWriter::log(ss.str(), LogWriter::Info);
		}
	}
}

template<typename T_type>
inline void TomlFileManager::getValue(std::string key, T_type& value_to_be_changed) {
	getValueOnCondition(key, value_to_be_changed, true);
}

template<typename T_type>
inline void TomlFileManager::getValueOrExit(std::string key, T_type& value_to_be_changed) {
	getValue(key, value_to_be_changed);
	if (has_getValueFailed) {
		LogWriter::logAndPrintError("get value failed when reading key: " + key + " @TomlFileManager::getValueOrExit\n");
		exit(-1);
	}
}

// 可选参数。toml文件中可以有，可以没有
template<typename T_type>
inline void TomlFileManager::getValueIfExists(std::string key, T_type& value_to_be_changed) {
	getValueOnCondition(key, value_to_be_changed, false);
}

#endif // TOML_FILE_MANAGER_H


