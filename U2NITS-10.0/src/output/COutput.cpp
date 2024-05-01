#include "COutput.h"
#include "../global/FilePathManager.h"


void U2NITS::COutput::updateFileName(int istep) {
	if (!initialized) {
		outputPathWithSlash = FilePathManager::getInstance()->getOutputDirectory();
		basicFileName = GlobalPara::basic::filename;// Ӧ���ڶ�ȡCConfig����ã���ֹfilenameδ��ȡ��ʹ��
		initialized = true;
	}

	// �����ļ���
	char szBuffer[20];
	sprintf_s(szBuffer, _countof(szBuffer), "%04d", istep);//��istepС��4λ������0
	std::string str_istep_withBracket = "[" + std::string(szBuffer) + "]";
	tecplotFilePath = outputPathWithSlash + basicFileName + str_istep_withBracket + ".dat";
	tecplotBoundaryFilePath = outputPathWithSlash + basicFileName + "_tec_boundary.dat";
	continueFilePath = outputPathWithSlash + "pause_" + basicFileName + str_istep_withBracket + ".dat";
	recoveryFilePath = outputPathWithSlash + "recovery_" + basicFileName + ".dat";
	continueFilePath_nan = outputPathWithSlash + "nan_" + basicFileName + str_istep_withBracket + ".dat";
}
