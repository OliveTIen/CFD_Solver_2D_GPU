#include "Config.h"
#include "FilePathManager.h"
#include <fstream>
#include <iostream>
#include "../output/LogWriter.h"

rapidjson::Document Config::config = rapidjson::Document();

void Config::writeConfig(std::string filename) {
	//变量赋给config
	config["lastfilename"].SetString(filename.c_str(), (rapidjson::SizeType)filename.size());//filename->"lastfilename"

	//config写入json https://blog.csdn.net/yang_aq/article/details/116934216
	rapidjson::StringBuffer buffer;
	rapidjson::PrettyWriter<rapidjson::StringBuffer> prettyWriter(buffer);//PrettyWriter的大括号格式化了更好看
	config.Accept(prettyWriter);
	std::string content = buffer.GetString();
	std::string exePath_withSlash = FilePathManager::getInstance()->getExePath_withSlash();
	std::ofstream outfile(exePath_withSlash + ".config\\config.json");
	if (outfile.is_open()) {
		outfile << content;
		outfile.close();
	}
	else std::cout << "Error: fail to open config.txt, in void Config::writeConfig()\n";

}

void Config::readConfig() {
	//读取json文件，赋给config
	std::string exePath_withSlash = FilePathManager::getInstance()->getExePath_withSlash();
	std::ifstream inf(exePath_withSlash + ".config\\config.json");
	if (inf.is_open()) {
		std::string json_content((std::istreambuf_iterator<char>(inf)), std::istreambuf_iterator<char>()); //将文件的数据流转为std::string类型
		inf.close();
		config.Parse(json_content.c_str());
	}
	else {
		//读取失败，则使用默认。注意，此时不必WriteConfig，因为Config要在获取文件名之后Write
		LogWriter::writeLog("Prompt: fail to open config.json, use backup config. (GlobalStatic::readConfig)\n");
		useBackupConfig();
	}
	//config的值 等后面需要的时候再读取
	//

}

void Config::useBackupConfig() {
	//ch_json赋给config
//由于是新生成的文件，用[NULL]指示
	const char* ch_json = R"({
		"lastfilename":"[NULL]",
		"version":"1.0"
		})";
	config.Parse(ch_json);
}
