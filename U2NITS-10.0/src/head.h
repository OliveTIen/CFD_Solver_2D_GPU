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

//namespace globalpara {}//�����ռ䶨��ı����ᵼ���ض������⣬���ҽ�����Ҳ���ɡ�����̬��������


class global {
public:
	static rapidjson::Document config;//json�ļ�
	static std::string filename;
	static std::string exePath;//һ��ʼ�Ѿ���ʼ��������ʹ��
	static int autosaveFileNum;

	static int flag_reconstruct;//�ع� ��ʱ���ţ���Ϊ��֪����ι���

	//�����û�����
	static void getFileName_UseUserInput();
	//�޳�' ''\t'','
	static std::vector<std::string> splitString(std::string);
	//��ʼ��exe�ļ�·������һ��ʼ�����
	static void iniExePath();
	//��ȡexe�ļ�·�������ѻ�á�flag:1-����б��(\\)0-����
	static std::string getExePath(int flag=1);
	//�о��ļ����ơ�path����Ϸ�б��
	static std::vector<std::string> ls(std::string path);
	//�����ļ��С����Ѵ��ڣ��򲻴���������exePath��ʼ����ʹ��
	static void createFolder(std::string foldername);
	//stringsתints
	static std::vector<int> Words2Ints(std::vector<std::string> words);
	//����ת��Ϊ{}��ʾ����ʽ
	static std::string DoubleArray_2_String(double* U, int length) {
		std::string str;
		str += "{";
		for (int i = 0; i < length; i++) {
			str += std::to_string(U[i]) + ",";
		}
		str[str.length() - 1] = '}';//������","�滻��"}"
		return str;
	}
	//app 1-��д 0-���� ����filename��ʼ����ʹ��
	static void writeLog(std::string, bool app = 1);
	static void writeLogAndCout(std::string);
	//��ǰ����ʱ��
	static std::string currentDateTime();
	//�����ļ�
	static void writeConfig();
	//�����ļ�
	static void readConfig();
	//����
	static void useBackupConfig();
};

namespace cmd {
	void printHeader();
	//��ȡ����̨���λ�á�x�ᳯ�ң�y�ᳯ��
	COORD getCursorPosition();
	//��ȡ����̨���ڻ�������С
	COORD getScrnInfo();
	//���ÿ���̨���λ��
	void setCursorPosition(COORD pstn);
	//���p1, p2֮��Ŀ���̨���ݣ�������궨λ��p1
	void clearDisplay(COORD p1, COORD p2);
	//���ƽ�������0<=percent<=100
	void drawProgressBar(double percent);
}



#endif
