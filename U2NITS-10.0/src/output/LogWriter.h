
#ifndef LOG_WRITER_H
#define LOG_WRITER_H

#include <string>

class LogWriter {
	
private:
	//static std::string m_fullFilePath;
	static std::string dateAndTime;

public:
	//static void setLogPath(std::string fullFilePath) { m_fullFilePath = fullFilePath; }
	//app 1-续写 0-覆盖 须在filename初始化后使用
	static void writeLog(std::string content, bool app = 1);
	static void writeLogAndCout(std::string content);
};

#endif