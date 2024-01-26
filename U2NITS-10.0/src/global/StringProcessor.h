
#ifndef STRING_PROCESSOR_H
#define STRING_PROCESSOR_H

#include <vector>
#include <string>

class StringProcessor {
public:
	//// 字符串数组转换
	// 字符串数组转int数组
	static std::vector<int> stringVector_2_intVector(std::vector<std::string> words);

	//// 字符串转换
	// 字符串转字符串数组。按空格和制表符分割
	static std::vector<std::string> splitString(std::string tLine);

	//// 字符串修改
	// 字符串中'='替换为" = " 源自CFD-bighw
	static std::string replaceCharInString(std::string word, char oldChar = '=', std::string newChar = " = ");
	// 删除flag及其后面的字符
	static std::string deleteCommentByChar(std::string word, char flag = '%');

	//// 字符串生成
	// 时间转字符串。输入：秒，输出：时:分:秒
	static std::string timeFormat(int dt);
	// double数组转字符串。数组转化为{}表示的形式
	static std::string doubleArray_2_string(double* U, int length);

};

#endif // !STRING_PROCESSOR_H
