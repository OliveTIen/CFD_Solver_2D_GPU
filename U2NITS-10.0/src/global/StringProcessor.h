
#ifndef STRING_PROCESSOR_H
#define STRING_PROCESSOR_H

#include <vector>
#include <string>

class StringProcessor {
public:
	// �ַ�������->int����
	static std::vector<int> Words2Ints(std::vector<std::string> words);
	// ������顣����ת��Ϊ{}��ʾ����ʽ
	static std::string DoubleArray_2_String(double* U, int length);
	// ���ո���Ʊ���ָ�
	static std::vector<std::string> splitString(std::string tLine);
};

#endif // !STRING_PROCESSOR_H
