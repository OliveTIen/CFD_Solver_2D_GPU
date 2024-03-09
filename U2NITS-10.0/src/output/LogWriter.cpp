#include "LogWriter.h"

#include "../global/FilePathManager.h"
#include "../global/SystemInfo.h"

std::ofstream LogWriter::m_logFile = std::ofstream();
std::string LogWriter::m_logFilePath = std::string();
bool LogWriter::m_hasInitilized = false;
LogWriter::LogLevel LogWriter::m_logLevel = LogWriter::LogLevel::Debug;
LogWriter::LogLevel LogWriter::m_coutLevel = LogWriter::LogLevel::Info;

std::string LogWriter::enumToString(LogLevel level) {
	switch (level) {
	case LogLevel::Fatal:
		return "Fatal";
	case LogLevel::Error:
		return "Error";
	case LogLevel::Warning:
		return "Warning";
	case LogLevel::Info:
		return "Info";
	case LogLevel::Debug:
		return "Debug";
	default:
		return "Unknown";
	}
}

void LogWriter::writeLog(std::string content, LogLevel logLevel) {
	if (!m_hasInitilized) {
		std::string directory = FilePathManager::getInstance()->getOutputDirectory();
		std::string dateAndTime = SystemInfo::getCurrentDateTime_suitableForFileName();
		m_logFilePath = directory + dateAndTime + ".LOG";
		m_hasInitilized = true;
	}
	if (logLevel <= m_logLevel) {// logLevelС�ڲο���˵����������д����־
		m_logFile.open(m_logFilePath, std::ios::app);
		m_logFile << content;
		m_logFile.close();
	}

}

void LogWriter::writeLogAndCout(std::string content, LogLevel logLevel, LogLevel coutLevel) {
	writeLog(content);
	if (coutLevel <= m_coutLevel) {
		std::cout << content;
	}
}
