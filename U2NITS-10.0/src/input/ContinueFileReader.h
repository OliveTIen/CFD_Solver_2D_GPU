#ifndef _CONTINUEFILEREADER_H_
#define _CONTINUEFILEREADER_H_

#include <string>
#include <fstream>
/// <summary>
///  ����FVM_2D::readContinueFile()����������ָ����ʼ���Ȳ������ܸ��ӣ�
/// Ŀǰ�Ȳ�æ��ɸ���
/// </summary>
class ContinueFileReader {
private:
	bool m_filepathHasInitialized = false;
	std::string m_filepath;
	std::ofstream m_outFile;

public:
	ContinueFileReader() {}
	ContinueFileReader(std::string filepath) {
		setFilePath(filepath);
	}

	void setFilePath(std::string filepath) {
		m_filepath = filepath;
		m_filepathHasInitialized = true;
	}

	//void read() {

	//}

};

#endif // _CONTINUEFILEREADER_H_