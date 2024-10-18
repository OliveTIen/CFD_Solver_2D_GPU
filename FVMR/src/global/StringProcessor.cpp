#include "StringProcessor.h"

std::vector<int> StringProcessor::stringVector_2_intVector(std::vector<std::string> words) {
	std::vector<int> ints;
	for (int i = 0; i < words.size(); i++) {
		ints.push_back(std::stoi(words[i]));
	}
	return ints;
}

std::string StringProcessor::doubleArray_2_string(myfloat* U, int length) {
	std::string str;
	str += "{";
	for (int i = 0; i < length; i++) {
		str += std::to_string(U[i]) + ",";
	}
	str[str.length() - 1] = '}';//将最后的","替换成"}"
	return str;
}

std::vector<std::string> StringProcessor::splitString(std::string tLine) {
	tLine = tLine + ' ';//目的是防止最后一个单词没有被push_back进tWords
	std::vector<std::string> tWords;
	std::string tWord;
	bool isCurrentMeaningful = 0;
	bool isLastMeaningful = 0;
	const int length = (int)tLine.size();
	//std::cout << "tLine= " << tLine << std::endl;
	//std::cout << "size= " << length << std::endl; 成功返回tLine长度
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

std::string StringProcessor::replaceCharInString(std::string word, char oldChar, std::string newChar) {
	std::string newString;
	for (int i = 0; i < word.size(); i++) {
		if (word[i] == oldChar) {
			newString = newString + newChar;
		}
		else {
			newString = newString + word[i];
		}
	}
	return newString;
}

std::string StringProcessor::deleteCommentByChar(std::string word, char flag) {
	// 删除flag及其后面的字符
	int i_last_meaningful_char = -1;
	for (int i = 0; i < word.size(); i++) {
		// 遇到'%'停止
		if (word[i] == flag) {
			break;
		}
		// 记录位置
		i_last_meaningful_char = i;
	}
	// 若整行被删去，则返回空字符串
	if (i_last_meaningful_char == -1) {
		return std::string();
	}
	return word.substr(0, i_last_meaningful_char + 1);
}

std::string StringProcessor::timeFormat(int dt) {
	if (dt < 0) {
		return "0:00:00";
	}

	int h = dt / 3600;//专门用int的除法
	dt -= h * 3600;
	int m = dt / 60;
	dt -= m * 60;

	std::string str;
	str += std::to_string(h);
	char szBuffer[4]{};
	sprintf_s(szBuffer, _countof(szBuffer), "%02d", m);
	str += std::string(":") + szBuffer;
	sprintf_s(szBuffer, _countof(szBuffer), "%02d", dt);
	str += std::string(":") + szBuffer;
	return str;
}
