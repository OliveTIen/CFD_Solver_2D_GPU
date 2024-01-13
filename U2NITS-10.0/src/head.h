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

//namespace globalpara {}//�����ռ䶨��ı����ᵼ���ض������⣬���ҽ�����Ҳ���ɡ�����̬��������


class GlobalStatic {
public:
	static std::string filename;
	static std::string exePath_withSlash;//��б�ܡ�һ��ʼ�Ѿ���ʼ��������ʹ��
	static int autosaveFileNum;
	static int flag_reconstruct;//�ع� ��ʱ���ţ���Ϊ��֪����ι���

	//�����û�����
	static void getFileName_UseUserInput();
	//�޳�' ''\t'','
	static std::vector<std::string> splitString(std::string);
	//��ʼ��exe�ļ�·������һ��ʼ�����
	static void iniExePath();
	//��ȡexe�ļ�·�������ѻ�á�ʼ���з�б��
	static std::string getExePath_withSlash(int flag=1);
	//�о��ļ����ơ�path����Ϸ�б��
	static std::vector<std::string> ls(std::string path);
	//�����ļ��С����Ѵ��ڣ��򲻴���������exePath��ʼ����ʹ��
	static void createFolderIfDoesntExist(std::string foldername);



};



#endif
