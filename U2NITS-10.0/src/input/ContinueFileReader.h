#ifndef _CONTINUEFILEREADER_H_
#define _CONTINUEFILEREADER_H_

#include <string>
/// <summary>
///  ����FVM_2D::readContinueFile()����������ָ����ʼ���Ȳ������ܸ��ӣ�
/// Ŀǰ�Ȳ�æ��ɸ���
/// </summary>
class ContinueFileReader {
private:
	bool filepathHasInitialized = false;
	std::string m_filepath;


public:
	ContinueFileReader(std::string filepath) {
		setFilePath(filepath);
	}

	void setFilePath(std::string filepath) {
		m_filepath = filepath;
		filepathHasInitialized = true;
	}

	//void read() {

	//}

};

#endif // _CONTINUEFILEREADER_H_