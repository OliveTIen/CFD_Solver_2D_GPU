#ifndef SYSTEM_INFO_H
#define SYSTEM_INFO_H
#include <string>

class SystemInfo {
public:
	static std::string getCurrentDateTime();
	// ȥ���Ƿ��ַ�":"�ȣ��ʺ������ļ���
	static std::string getCurrentDateTime_suitableForFileName();
};

#endif