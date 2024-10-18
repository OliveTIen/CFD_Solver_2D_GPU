#include "HistWriter.h"
#include <fstream>
#include <iostream>
#include "LogWriter.h"
#include "../global/GlobalPara.h"

void HistWriter::writeHistFile_2(int iteration, const double* residual, int length) {
	if (!hasBeenCalled_writeHistFile_2) {
		// 第一次调用时执行
		if (!GlobalPara::basic::_continue) {// 非续算，即第一次算。要写文件头
			writeHistFileHead_2();
		}
		hasBeenCalled_writeHistFile_2 = true;
	}

	if (length != 4) {
		LogWriter::logAndPrintError("Invalid vector size, @HistWriter::writeHistFile_2 \n");
		exit(4);
	}

	m_outfile.open(m_filePath, std::ios::app);
	if (m_outfile.fail()) {
		LogWriter::logAndPrintError("Write hist file data failed.\n");
		exit(-1);
	}
	m_outfile << std::scientific << iteration << " " << residual[0] << " " << residual[1] << " " << residual[2] << " " << residual[3] << "\n";
	m_outfile << std::defaultfloat;
	m_outfile.close();
}

void HistWriter::writeHistFileHead_2() {
	m_outfile.open(m_filePath);
	if (m_outfile.fail()) {
		LogWriter::logAndPrintError("Write hist file head failed.\n");
		exit(-1);
	}
	m_outfile << R"(TITLE=")" << GlobalPara::basic::filename << R"(")" << "\n";
	m_outfile << R"(VARIABLES="Iteration" "Residual_rho" "Residual_rhou" "Residual_rhov" "Residual_rhop")" << "\n";
	m_outfile << R"(ZONE  F=POINT)" << "\n";
	m_outfile.close();
}
