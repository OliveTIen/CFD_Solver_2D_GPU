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
//#include "./rapidjson/document.h"
//#include "./rapidjson/prettyWriter.h"
//#include "./rapidjson/stringbuffer.h"
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


class global {
public:
	static rapidjson::Document config;//json文件
	static std::string filename;
	static std::string exePath;//一开始已经初始化，放心使用
	static int autosaveFileNum;

	static int flag_reconstruct;//重构 暂时留着，因为不知道如何归类

	//接受用户输入
	static void getFileName_UseUserInput();
	//剔除' ''\t'','
	static std::vector<std::string> splitString(std::string);
	//初始化exe文件路径。在一开始便调用
	static void iniExePath();
	//获取exe文件路径，若已获得。flag:1-加上斜杠(\\)0-不加
	static std::string getExePath(int flag=1);
	//列举文件名称。path需加上反斜杠
	static std::vector<std::string> ls(std::string path);
	//创建文件夹。若已存在，则不创建。须在exePath初始化后使用
	static void createFolder(std::string foldername);
	//strings转ints
	static std::vector<int> Words2Ints(std::vector<std::string> words);
	//数组转化为{}表示的形式
	static std::string DoubleArray_2_String(double* U, int length) {
		std::string str;
		str += "{";
		for (int i = 0; i < length; i++) {
			str += std::to_string(U[i]) + ",";
		}
		str[str.length() - 1] = '}';//将最后的","替换成"}"
		return str;
	}
	//app 1-续写 0-覆盖 须在filename初始化后使用
	static void writeLog(std::string, bool app = 1);
	static void writeLogAndCout(std::string);
	//当前日期时间
	static std::string currentDateTime();
	//配置文件
	static void writeConfig();
	//配置文件
	static void readConfig();
	//备用
	static void useBackupConfig();
};

namespace cmd {
	void printHeader();
	//获取控制台光标位置。x轴朝右，y轴朝下
	COORD getCursorPosition();
	//获取控制台窗口缓冲区大小
	COORD getScrnInfo();
	//设置控制台光标位置
	void setCursorPosition(COORD pstn);
	//清除p1, p2之间的控制台内容，并将光标定位于p1
	void clearDisplay(COORD p1, COORD p2);
	//绘制进度条。0<=percent<=100
	void drawProgressBar(double percent);
}



#endif
