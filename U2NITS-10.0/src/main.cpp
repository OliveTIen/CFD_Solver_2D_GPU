﻿#include "drivers/CDriver.h"
#include "global/GlobalPara.h"
#include "global/FilePathManager.h"
#include "global/SystemInfo.h"
#include "input/TomlFileManager.h"
#include "output/ConsolePrinter.h"
#include "output/LogWriter.h"

#include <stdlib.h>
using namespace std;

int main() {
	std::cout << "Using VS Profiler. Please remove '/Profile' option on release.\n";

	ConsolePrinter::printHeader();
	LogWriter::writeLog(SystemInfo::getCurrentDateTime() + "\n");
	// 读取input.toml输入参数，处理部分输入参数
	TomlFileManager::getInstance()->readFileAndParseFile(FilePathManager::getInstance()->getTomlFilePath());
	if (GlobalPara::basic::dimension == 2) {
		U2NITS::CDriver::run_GPU();
	}
	else {
		LogWriter::writeLogAndCout("Error: invalid dimension type.\n", LogWriter::Error, LogWriter::Error);
	}


	// 停止
#ifdef _WIN32
	system("pause");
#elif defined __linux__
	getchar();
#endif
	return 0;
}


//命名规范：
//ini_x初始化变量的值
//x_to_x计算并输出。输出方式可以为return、&
//set_x改变成员变量的值
//get_x输出成员变量的值
//cal_x计算并改变成员变量的值，不输出



//目前存在的问题
//1.[已解决]读取续算文件时，有时候并没有读取新的文件
//2.算不对。等熵涡出现发散问题。然而老师的程序是可以的。可能：1.时间推进不是rk3而是显式(已证否) 2.不是线性重构
// 可是等熵涡的边界只有周期边界，排除了边界的干扰。可能是计算黎曼通量有问题。
//3.新bug：在计算过程中，本来想按esc停止，但是误按了esc下面的"`"键，导致程序界面消失但仍然在计算。