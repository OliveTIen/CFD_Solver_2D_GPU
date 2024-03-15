#ifndef TOML_FILE_MANAGER_H
#define TOML_FILE_MANAGER_H
#include <string>
#include "../include/toml/cpptoml.h"

class TomlFileManager {
	enum ErrorCode {
		noError,
		getValueFailed,
		getOptionalValueFailed
	};
	class ErrorMessage {
		ErrorCode m_code;
		std::string m_message;
	public:
		ErrorMessage(ErrorCode code, std::string message) {
			m_code = code;
			m_message = message;
		}
		ErrorMessage(ErrorCode code, std::string message, int lineNumber) {
			m_code = code;
			m_message = message + "(line " + std::to_string(lineNumber) + ")";
		}
		ErrorCode getErrorCode() { return m_code; }
		std::string getMessage() { return m_message; }
	};

private:
	std::shared_ptr<cpptoml::table>m_parsedFile = nullptr;
	bool has_getValueFailed = false;
	bool has_getOptionalValueFailed = false;

public:
	void readFileAndParseFile(std::string fullFilePath);
	void printParsedFile();

private:
	// only used in readFileAndParseFile()
	void modifyGlobalParametersAccordingToParsedFile();
	// �޸�һЩ������ֵ
	void initialize_ruvp();
	void logErrorMessage(ErrorMessage errorMessage);

	// �����ڣ����value_to_be_changed��ֵ�����򲻱�
	// ע�ⲻʶ��float��ʽ
	template <typename T_type>
	void getValue(std::string key, T_type& value_to_be_changed);
	template <typename T_type>
	void getOptionalValueIfExists(std::string key, T_type& value_to_be_changed);
	

};

template<typename T_type>
inline void TomlFileManager::getValue(std::string key, T_type& value_to_be_changed) {
	cpptoml::option<T_type> value_optional = m_parsedFile->get_qualified_as<T_type>(key);
	// �����ڣ���ı�value_to_be_changed��ֵ�����򱨴��˳�
	if (value_optional) {
		//std::cout << "Debug: key=" << key << ", old_value=" << value_to_be_changed << ", new_value=";
		value_to_be_changed = *value_optional;
		//std::cout << value_to_be_changed << "\n";
	}
	else {
		has_getValueFailed = true;
		std::string message = "Get value failed, key = " + key;
		ErrorMessage m(ErrorCode::getValueFailed, message);
		logErrorMessage(m);
	}

}

template<typename T_type>
inline void TomlFileManager::getOptionalValueIfExists(std::string key, T_type& value_to_be_changed) {
	cpptoml::option<T_type> value_optional = m_parsedFile->get_qualified_as<T_type>(key);
	// �����ڣ���ı�value_to_be_changed��ֵ
	if (value_optional) {
		value_to_be_changed = *value_optional;
	}
	else {
		has_getOptionalValueFailed = true;
		std::string message = "Get optional value failed, key = " + key;
		ErrorMessage m(ErrorCode::getOptionalValueFailed, message);
		logErrorMessage(m);
	}

}


#endif // TOML_FILE_MANAGER_H


