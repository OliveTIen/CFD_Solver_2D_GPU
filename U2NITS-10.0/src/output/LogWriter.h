
#ifndef LOG_WRITER_H
#define LOG_WRITER_H

#include <string>

class LogWriter {
	
private:
	//static std::string m_fullFilePath;
	static std::string dateAndTime;

public:
	//static void setLogPath(std::string fullFilePath) { m_fullFilePath = fullFilePath; }
	//app 1-��д 0-���� ����filename��ʼ����ʹ��
	static void writeLog(std::string content, bool app = 1);
	static void writeLogAndCout(std::string content);
};

#endif