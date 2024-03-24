#include "LogWriter.h"

#include "../global/FilePathManager.h"
#include "../global/SystemInfo.h"
#include "../global/StringProcessor.h"

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
	if (logLevel <= m_logLevel) {// logLevel小于参考，说明更紧急，写入日志
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

void LogWriter::writeBoundaryCondition(double* inlet_ruvp, double* outlet_ruvp, double* inf_ruvp, int num_ruvp) {
	// 日志记录边界参数
	std::string str;
	str += "BoundaryCondition:\n";
	str += "inlet::ruvp\t" + StringProcessor::doubleArray_2_string(inlet_ruvp, num_ruvp)
		+ "\noutlet::ruvp\t" + StringProcessor::doubleArray_2_string(outlet_ruvp, num_ruvp)
		+ "\ninf::ruvp\t" + StringProcessor::doubleArray_2_string(inf_ruvp, num_ruvp)
		+ "\n";
	writeLog(str);
}
