#ifndef CONFIG_H
#define CONFIG_H
#include <string>
#include "../include/rapidjson/document.h"
#include "../include/rapidjson/prettyWriter.h"
#include "../include/rapidjson/stringbuffer.h"

class Config {
public:
	static rapidjson::Document config;//json�ļ�
	//�����ļ�
	static void writeConfig(std::string filename);
	//�����ļ�
	static void readConfig();
	//����
	static void useBackupConfig();
};

#endif