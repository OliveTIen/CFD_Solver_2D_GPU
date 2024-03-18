#ifndef _JSONCONFIG_H_
#define _JSONCONFIG_H_
#include <string>
#include "../include/rapidjson/document.h"
#include "../include/rapidjson/prettyWriter.h"
#include "../include/rapidjson/stringbuffer.h"

class JsonConfig {
public:
	static rapidjson::Document config;//json文件
	//配置文件
	static void writeConfig(std::string filename);
	//配置文件
	static void readConfig();
	//备用
	static void useBackupConfig();
};

#endif