#ifndef TOML_FILE_MANAGER_H
#define TOML_FILE_MANAGER_H
#include <string>
#include "../include/toml/cpptoml.h"

class TomlFileManager {
private:
	std::shared_ptr<cpptoml::table>m_parsedFile = nullptr;
	bool has_error_on_modifying = false;

public:
	void readFileAndParseFile(std::string fullFilePath);
	void printParsedFile();

private:
	// only used in readFileAndParseFile()
	void modifyGlobalParametersAccordingToParsedFile();
	// 修改一些变量的值
	void initialize_ruvp();

	// 若存在，则给value_to_be_changed赋值，否则不变
	// 注意不识别float格式
	template <typename T_type>
	void getValueIfExists(std::string key, T_type& value_to_be_changed);

};

#endif

template<typename T_type>
inline void TomlFileManager::getValueIfExists(std::string key, T_type& value_to_be_changed) {
	cpptoml::option<T_type> value_optional = m_parsedFile->get_qualified_as<T_type>(key);
	// 若存在，则改变value_to_be_changed的值
	if (value_optional) {
		//std::cout << "Debug: key=" << key << ", old_value=" << value_to_be_changed << ", new_value=";
		value_to_be_changed = *value_optional;
		//std::cout << value_to_be_changed << "\n";
	}
	else {
		std::cout << "Warning: Cannot read key \"" << key << "\" or its value.\n";
		has_error_on_modifying = true;
	}

}
