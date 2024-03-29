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
	bool hasBeenCalled_writeHistFile_2 = false;

public:
	HistWriter() { m_filePath = "hist.dat"; }
	HistWriter(std::string fullFilePath) { m_filePath = fullFilePath; }
	void setFilePath(std::string fullFilePath){ m_filePath = fullFilePath; }
	void writeHistFileHead();
	void writeHistFileData(int iteration, const double* residual, int length);

	// 新函数，将writeHead集成在里面
	void writeHistFile_2(int iteration, const double* residual, int length);
	void writeHistFileHead_2();
};

#endif