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
	// ��Ա����
	static TomlFileManager* classPointer;
	std::shared_ptr<cpptoml::table>m_tree = nullptr;
	bool has_getValueFailed = false;

public:
	static TomlFileManager* getInstance();
	void readTomlFile(std::string fullFilePath);
	
	// �ӿڣ���������ᱨ���˳�
	template <typename T_type>
	void getValueOrExit(std::string key, T_type& value_to_be_changed);// ��ѡ

	// ���º��������������󲻻������˳���ֻ�ǽ�has_getValueFailed��Ϊtrue
	template <typename T_type>
	void getValueOnCondition(std::string key, T_type& value_to_be_changed, bool condition);// condition==true�����ѡ�������ѡ
	template <typename T_type>
	void getValue(std::string key, T_type& value_to_be_changed);// ��ѡ
	template <typename T_type>
	void getValueIfExists(std::string key, T_type& value_to_be_changed);// ��ѡ
	

private:
	TomlFileManager() {};
	// m_treeת��Ϊȫ�ֱ���
	void treeToGlobalParameter();
	// [todo]
	void treeToSGlobalPara(SGlobalPara::SGlobalPara& spara);
	// ͳһ�໥ì�ܵ�����
	void handleConflictingInputs();
	void getValue_boundaryCondition2D(std::string parent, myfloat ruvp[4]);

	// tree�Ƿ����key�����������"."�ָ��㼶
	bool treeContainsKey(std::string key) { return m_tree->contains_qualified(key); }
	// ��ӡtree�Ƿ����key
	void printTreeContainsKeyOrNot(std::string key);
	// ��ӡtree
	void printTree() { std::cout << *m_tree << std::endl; }
	// ���has_getValueFailed�Ƿ�Ϊtrue���˳�
	void ifFailedThenExit();
	
	template <typename T_type>
	void getValueOrigin(std::string key, T_type& value_to_be_changed);
};

template<typename T_type>
inline void TomlFileManager::getValueOrigin(std::string key, T_type& value_to_be_changed) {

	bool failed = false;
	// �����bool����
	if (std::is_same<T_type, bool>::value) {
		// ���Զ�ȡ�ַ���'bool'���������ţ�
		cpptoml::option<bool> valueBool = m_tree->get_qualified_as<bool>(key);
		if (valueBool) {
			value_to_be_changed = *valueBool;
		}
		// ����޷���ȡ�ַ���'bool'����int���Ͷ�ȡ����
		else {
			cpptoml::option<int> valueInt = m_tree->get_qualified_as<int>(key);
			if (valueInt) {
				value_to_be_changed = (bool)*valueInt;// intת��Ϊbool����ʹ��tomlTree�д洢��valueInt��int���͵�2��-1��������value_to_be_changed��bool���ͣ������ʽת��Ϊ1
			}
			else {
				failed = true;
			}
		}
	}
	// cpptoml��ʶ��float��ʽ��Ӧ����double��ȡ
	else if (std::is_same<T_type, float>::value) {
		cpptoml::option<double> valueDouble = m_tree->get_qualified_as<double>(key);
		if (valueDouble) {
			value_to_be_changed = *valueDouble;
		}
		else {
			failed = true;
		}
	}
	// ������������ͣ�������ȡ
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
		if (condition) {// ��������������ѡ����˻ᱨ��
			LogWriter::logAndPrintError("Required parameter \"" + key + "\" doesn't exist.\n");
			has_getValueFailed = true;
		}
		else {// �����������ѡ����˲��ᱨ��
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

// ��ѡ������toml�ļ��п����У�����û��
template<typename T_type>
inline void TomlFileManager::getValueIfExists(std::string key, T_type& value_to_be_changed) {
	getValueOnCondition(key, value_to_be_changed, false);
}

#endif // TOML_FILE_MANAGER_H


