
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
	static std::string enumToStringWithBracket(LogLevel level);

public:
	//static void setLogPath(std::string fullFilePath) { m_fullFilePath = fullFilePath; }
	//app 1-续写 0-覆盖 须在filename初始化后使用
	static void log(std::string content, LogLevel logLevel = Info);
	static void print(std::string content, LogLevel coutLevel = Info);
	static void logAndPrint(std::string content, LogLevel logLevel = Info, LogLevel coutLevel = Info);
	static void logError(std::string content);
	static void printError(std::string content);
	static void logAndPrintError(std::string content);

	static void writeBoundaryCondition(double* inlet_ruvp, double* outlet_ruvp,double* inf_ruvp, int num_ruvp);
	static void logAndPrintSignalInfo(int signal);
};

#endif