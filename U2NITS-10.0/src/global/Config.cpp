#include "Config.h"
#include "FilePathManager.h"
#include <fstream>
#include <iostream>
#include "../output/LogWriter.h"

rapidjson::Document Config::config = rapidjson::Document();

void Config::writeConfig(std::string filename) {
	//��������config
	config["lastfilename"].SetString(filename.c_str(), (rapidjson::SizeType)filename.size());//filename->"lastfilename"

	//configд��json https://blog.csdn.net/yang_aq/article/details/116934216
	rapidjson::StringBuffer buffer;
	rapidjson::PrettyWriter<rapidjson::StringBuffer> prettyWriter(buffer);//PrettyWriter�Ĵ����Ÿ�ʽ���˸��ÿ�
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
	//��ȡjson�ļ�������config
	std::string exePath_withSlash = FilePathManager::getInstance()->getExePath_withSlash();
	std::ifstream inf(exePath_withSlash + ".config\\config.json");
	if (inf.is_open()) {
		std::string json_content((std::istreambuf_iterator<char>(inf)), std::istreambuf_iterator<char>()); //���ļ���������תΪstd::string����
		inf.close();
		config.Parse(json_content.c_str());
	}
	else {
		//��ȡʧ�ܣ���ʹ��Ĭ�ϡ�ע�⣬��ʱ����WriteConfig����ΪConfigҪ�ڻ�ȡ�ļ���֮��Write
		LogWriter::writeLog("Prompt: fail to open config.json, use backup config. (GlobalStatic::readConfig)\n");
		useBackupConfig();
	}
	//config��ֵ �Ⱥ�����Ҫ��ʱ���ٶ�ȡ
	//

}

void Config::useBackupConfig() {
	//ch_json����config
//�����������ɵ��ļ�����[NULL]ָʾ
	const char* ch_json = R"({
		"lastfilename":"[NULL]",
		"version":"1.0"
		})";
	config.Parse(ch_json);
}
