#ifndef HIST_WRITER_H
#define HIST_WRITER_H
#include <string>
#include <fstream>
#include <vector>

class HistWriter {
private:
	std::ofstream m_outfile;
	std::string m_filePath;
	bool hasWritten = false;

public:
	HistWriter(std::string fullFilePath) { m_filePath = fullFilePath; }
	void writeHead();
	void writeData(int iteration, const double* residual, int length);
};

#endif