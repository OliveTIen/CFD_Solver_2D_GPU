
#ifndef STRING_PROCESSOR_H
#define STRING_PROCESSOR_H

#include <vector>
#include <string>

class StringProcessor {
public:
	// 字符串数组->int数组
	static std::vector<int> Words2Ints(std::vector<std::string> words);
	// 输出数组。数组转化为{}表示的形式
	static std::string DoubleArray_2_String(double* U, int length);
	// 按空格和制表符分割
	static std::vector<std::string> splitString(std::string tLine);
};

#endif // !STRING_PROCESSOR_H
