
#ifndef _COUTPUT_H_
#define _COUTPUT_H_

#include "ConsolePrinter.h"
#include "FieldWriter.h"
#include "HistWriter.h"
#include "LogWriter.h"
#include "ResidualCalculator.h"
#include "../global/CTimer.h"
#include "../global/FilePathManager.h"
namespace U2NITS {
	class COutput {
	private:
		std::string basicFileName;
		std::string outputPathWithSlash;
		bool initialized = false;// 上述两个string是否初始化

	public:
		std::string tecplotFilePath;
		std::string continueFilePath;
		std::string continueFilePath_nan;
		ConsolePrinter console;
		LogWriter log;
		HistWriter hist;
		FieldWriter field;
		CTimer timer;

		void updateFileName(int istep);
		std::string getDir() {
			return FilePathManager::getInstance()->getOutputDirectory();
		}
	};
}

#endif // _COUTPUT_H_