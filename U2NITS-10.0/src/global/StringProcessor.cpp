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
	str[str.length() - 1] = '}';//������","�滻��"}"
	return str;
}

std::vector<std::string> StringProcessor::splitString(std::string tLine) {
	tLine = tLine + ' ';//Ŀ���Ƿ�ֹ���һ������û�б�push_back��tWords
	std::vector<std::string> tWords;
	std::string tWord;
	bool isCurrentMeaningful = 0;
	bool isLastMeaningful = 0;
	const int length = (int)tLine.size();
	//std::cout << "tLine= " << tLine << std::endl;
	//std::cout << "size= " << length << std::endl; �ɹ�����tLine����
	for (int j = 0; j < length; j++) {
		if (tLine[j] != ' ' && tLine[j] != '\t' && tLine[j] != ',')
			isCurrentMeaningful = 1;
		else
			isCurrentMeaningful = 0;
		if (isCurrentMeaningful)
			tWord.push_back(tLine[j]);
		else if (isLastMeaningful) {
			tWords.push_back(tWord);
			tWord.clear();
		}
		isLastMeaningful = isCurrentMeaningful;
	}
	return tWords;

}
