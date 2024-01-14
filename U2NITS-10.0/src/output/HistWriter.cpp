#include "HistWriter.h"
#include <fstream>
#include <iostream>
#include "LogWriter.h"
#include "../GlobalPara.h"

void HistWriter::writeHead() {
	// 如果是续算则不必添加文件头
	if (GlobalPara::basic::_continue) {
		hasWritten = true;
		return;
	}
	m_outfile.open(m_filePath);
	if (m_outfile.fail()) {
		std::cout << "Error: Cannot open file at " << m_filePath << "\n";
		std::cout << "(Code: " << __FILE__ << ")\n";
		return;
	}
	m_outfile << R"(TITLE=")"<< GlobalPara::basic::filename << R"(")" << "\n";
	m_outfile << R"(VARIABLES="Iteration" "Residual_rho" "Residual_rhou" "Residual_rhov" "Residual_rhop")" << "\n";
	m_outfile << R"(ZONE  F=POINT)" << "\n";
	m_outfile.close();
	hasWritten = true;
}

void HistWriter::writeData(int iteration, const double* residual, int length) {
	if (!hasWritten) {
		std::cout << "Error: has't written file head. \n";
		std::cout << "(Code: " << __FILE__ << ")\n";
		return;
	}
	if (length != 4) {
		std::cout << "Error: Invalid vector size. \n";
		std::cout << "(Code: " << __FILE__ << ")\n";
		return;
	}
	m_outfile.open(m_filePath, std::ios::app);
	if (m_outfile.fail()) {
		std::cout << "Error: Cannot open file at " << m_filePath << "\n";
		std::cout << "(Code: " << __FILE__ << ")\n";
		return;
	}
	m_outfile << std::scientific << iteration << " " << residual[0] << " " << residual[1] << " " << residual[2] << " " << residual[3] << "\n";
	m_outfile << std::defaultfloat;
	m_outfile.close();
	hasWritten = true;
}
