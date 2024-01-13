#include "LogWriter.h"
#include <fstream>
#include <iostream>
#include "../global/FilePathManager.h"
#include "../global/SystemInfo.h"

std::string LogWriter::dateAndTime = std::string();

void LogWriter::writeLog(std::string content, bool app) {
	if (dateAndTime.empty()) {
		dateAndTime = SystemInfo::getCurrentDateTime_suitableForFileName();
	}

	//����filename��ʼ����ʹ�ã����������".LOG"
	std::ofstream f;
	// exePath_withSlash + "output\\" + ����ʱ�� + ".LOG"
	std::string exePath_withSlash = FilePathManager::getInstance()->getExePath_withSlash();
	std::string m_fullFilePath = exePath_withSlash + "output\\" + dateAndTime + ".LOG";
	if (app)f.open(m_fullFilePath, std::ios::app);
	else f.open(m_fullFilePath);
	f << content;
	f.close();
}

void LogWriter::writeLogAndCout(std::string content) {
	writeLog(content, 1);
	std::cout << content;
}
