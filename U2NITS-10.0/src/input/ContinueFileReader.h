#ifndef _CONTINUEFILEREADER_H_
#define _CONTINUEFILEREADER_H_

#include <string>
/// <summary>
///  由于FVM_2D::readContinueFile()函数包含对指针表初始化等操作，很复杂，
/// 目前先不忙完成该类
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