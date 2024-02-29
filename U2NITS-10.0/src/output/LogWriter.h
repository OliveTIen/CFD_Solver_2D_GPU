
#ifndef LOG_WRITER_H
#define LOG_WRITER_H

#include <string>
#include <fstream>
#include <iostream>

class LogWriter {
public:
	enum LogLevel {
		Fatal,
		Error,
		Warning,
		Info,
		Debug
	};
	
private:
	static std::ofstream m_logFile;
	static std::string m_logFilePath;
	static bool m_hasInitilized;
	static LogLevel m_logLevel;// 参考输出级别
	static LogLevel m_coutLevel;
	static std::string enumToString(LogLevel level);

public:
	//static void setLogPath(std::string fullFilePath) { m_fullFilePath = fullFilePath; }
	//app 1-续写 0-覆盖 须在filename初始化后使用
	static void writeLog(std::string content, LogLevel logLevel = Debug);
	static void writeLogAndCout(std::string content, LogLevel logLevel = Info, LogLevel coutLevel = Info);
};

#endif