#include "COutput.h"
#include "../global/FilePathManager.h"


void U2NITS::COutput::updateFileName(int istep) {
	if (!initialized) {
		outputPathWithSlash = FilePathManager::getInstance()->getOutputDirectory();
		basicFileName = GlobalPara::basic::filename;// 应当在读取CConfig后调用，防止filename未读取就使用
		initialized = true;
	}

	// 更新文件名
	char szBuffer[20];
	sprintf_s(szBuffer, _countof(szBuffer), "%04d", istep);//若istep小于4位数，则补0
	std::string str_istep_withBracket = "[" + std::string(szBuffer) + "]";
	tecplotFilePath = outputPathWithSlash + basicFileName + str_istep_withBracket + ".dat";
	continueFilePath = outputPathWithSlash + "pause_" + basicFileName + str_istep_withBracket + ".dat";
	continueFilePath_nan = outputPathWithSlash + "nan_" + basicFileName + str_istep_withBracket + ".dat";
}
