#ifndef SYSTEM_INFO_H
#define SYSTEM_INFO_H
#include <string>

class SystemInfo {
public:
	static std::string getCurrentDateTime();
	// 去除非法字符":"等，适合用作文件名
	static std::string getCurrentDateTime_suitableForFileName();
};

#endif