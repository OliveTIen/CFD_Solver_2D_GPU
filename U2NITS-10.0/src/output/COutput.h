
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
		std::string tecplotFilePath;// 流场输出
		std::string tecplotBoundaryFilePath;
		std::string continueFilePath;// 续算
		std::string continueFilePath_nan;
		std::string recoveryFilePath;
		std::string tecplot_hist_path;
		
		ConsolePrinter console;
		LogWriter log;
		HistWriter hist;
		CTimer timer;

		void updateFileName(int istep);
		void set_tecplot_hist_path(std::string path) { tecplot_hist_path = path; }
		std::string getDir() {
			return FilePathManager::getInstance()->getOutputDirectory();
		}

		FieldWriter* getFieldWriter() { return FieldWriter::getInstance(); }
	};
}

#endif // _COUTPUT_H_