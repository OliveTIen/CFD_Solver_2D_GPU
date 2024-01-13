#include "StringProcessor.h"

std::vector<int> StringProcessor::Words2Ints(std::vector<std::string> words) {
	std::vector<int> ints;
	for (int i = 0; i < words.size(); i++) {
		ints.push_back(std::stoi(words[i]));
	}
	return ints;
}

std::string StringProcessor::DoubleArray_2_String(double* U, int length) {
	std::string str;
	str += "{";
	for (int i = 0; i < length; i++) {
		str += std::to_string(U[i]) + ",";
	}
	str[str.length() - 1] = '}';//将最后的","替换成"}"
	return str;
}
