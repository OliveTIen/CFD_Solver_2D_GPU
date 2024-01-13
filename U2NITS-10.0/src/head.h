#ifndef HEAD_H
#define HEAD_H

#include <conio.h> 
#include <iostream>
#include <io.h>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#ifdef _WIN32
#include <direct.h>
#include <windows.h>
#elif defined __linux__
#include <sys/stat.h>
#include <unistd.h>
#endif
#include <stdlib.h> 
#include "./include/Eigen/Core"
#include "./include/Eigen/Dense"
#include <time.h>
#include <ctime>
#include <thread>
//#include "Math.h"
#include "./include/rapidjson/document.h"
#include "./include/rapidjson/prettyWriter.h"
#include "./include/rapidjson/stringbuffer.h"

#include "SignDefine.h"
#include "globalPara.h"

#ifdef _WIN32

#elif defined __linux__
typedef struct _COORD {
	short X;
	short Y;
} COORD, * PCOORD;
#endif

//namespace globalpara {}//命名空间定义的变量会导致重定义问题，而且仅声明也不可。但静态变量可以


class GlobalStatic {
public:
	static std::string filename;
	static std::string exePath_withSlash;//有斜杠。一开始已经初始化，放心使用
	static int autosaveFileNum;
	static int flag_reconstruct;//重构 暂时留着，因为不知道如何归类

	//接受用户输入
	static void getFileName_UseUserInput();
	//剔除' ''\t'','
	static std::vector<std::string> splitString(std::string);
	//初始化exe文件路径。在一开始便调用
	static void iniExePath();
	//获取exe文件路径，若已获得。始终有反斜杠
	static std::string getExePath_withSlash(int flag=1);
	//列举文件名称。path需加上反斜杠
	static std::vector<std::string> ls(std::string path);
	//创建文件夹。若已存在，则不创建。须在exePath初始化后使用
	static void createFolderIfDoesntExist(std::string foldername);



};



#endif
