#include "head.h"
#include "global/Config.h"
#include "output/LogWriter.h"
using namespace std;


std::string GlobalStatic::filename = "";
std::string GlobalStatic::exePath_withSlash = "";//exe所在路径，有斜杠(\\)
int GlobalStatic::autosaveFileNum = 3;
int GlobalStatic::flag_reconstruct = _REC_constant;


void GlobalStatic::getFileName_UseUserInput() {
	std::string lastf;
	if (Config::config.FindMember("lastfilename") != Config::config.MemberEnd()) {
		lastf = Config::config["lastfilename"].GetString();//将config的值赋给lastf
	}
	if (filename.length() == 0) {


		//输入提示
		std::cout << "Please input filename (without suffix):" << std::endl;
		if(lastf!="[NULL]")
			//若lastf!="[NULL]"则说明可以直接读取上次filename
			std::cout << "(Press a single \"Enter\" to use last filename [" << lastf << "])\n";
		else {
			//若lastf=="[NULL]"则说明是原本没有config文件，因此不能使用last filename
			//检测所有inp文件
		}
		//接受用户输入，作为文件名。若输入为空(即直接按enter)，则使用lastf
		char a = 0;
		std::string str;
		while (a != '\n') {
			a = getchar();
			str.push_back(a);
		}
		str.pop_back();//删除最后的\n
		if (str == "")str = lastf;
		filename = str;
	}
	else {
		std::cout << "Filename:" << filename << std::endl;
	}

}

std::vector<std::string> GlobalStatic::splitString(std::string tLine) {
	tLine = tLine + ' ';//目的是防止最后一个单词没有被push_back进tWords
	std::vector<std::string> tWords;
	std::string tWord;
	bool isCurrentMeaningful = 0;
	bool isLastMeaningful = 0;
	const int length = (int)tLine.size();
	//std::cout << "tLine= " << tLine << std::endl;
	//std::cout << "size= " << length << std::endl; 成功返回tLine长度
	for (int j = 0; j < length; j++) {
		if (tLine[j] != ' ' && tLine[j] != '\t' && tLine[j] != ',')
			isCurrentMeaningful = 1;
		else
			isCurrentMeaningful = 0;
		if (isCurrentMeaningful)
			tWord.push_back(tLine[j]);
		else if (isLastMeaningful) {
			tWords.push_back(tWord);
			tWord.clear();
		}
		isLastMeaningful = isCurrentMeaningful;
	}
	return tWords;
}

void GlobalStatic::iniExePath() {
	exePath_withSlash = getExePath_withSlash(1);
}

std::string GlobalStatic::getExePath_withSlash(int flag) {
	if (exePath_withSlash != "")return exePath_withSlash;

#ifdef _WIN32
	//获取路径并返回字符串
	const int maxPathLength = 260;
	char buffer[maxPathLength];
	_getcwd(buffer, maxPathLength);
	std::string str = buffer;
	str += std::string("\\");
	return str;


#elif defined __linux__
	char* buffer;
	buffer = getcwd(NULL, 0);
	std::string str = buffer;
	free(buffer);
	str += std::string("\\");
	return str;
#endif 
}

std::vector<std::string> GlobalStatic::ls(std::string path) {
	path += "*";
	_finddata64i32_t fileInfo;
	std::vector<std::string> files;
	intptr_t hFile = _findfirst(path.c_str(), &fileInfo);
	if (hFile == -1) {
		std::cout << "Error in GlobalStatic::ls()" << std::endl;
		return files;
	}
	do {
		files.push_back(fileInfo.name);
	} while (_findnext(hFile, &fileInfo) == 0);
	return files;
}

void GlobalStatic::createFolderIfDoesntExist(std::string foldername) {
	std::string path = GlobalStatic::exePath_withSlash;
	std::vector<std::string> files = GlobalStatic::ls(path);
	for (int i = 0; i < files.size(); i++) {
		if (files[i] == foldername)return;//文件夹已经存在
	}
	foldername = "mkdir " + path + foldername;
	system(foldername.c_str());
}





