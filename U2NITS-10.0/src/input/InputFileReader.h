#ifndef INPUTFILEREADER_H
#define INPUTFILEREADER_H

#include <fstream>
#include <iostream>
#include <vector>
#include <string>

// 将文本文件转化为字符串矩阵
class InputFileReader {
private:
	std::vector<std::vector<std::string>> tWordsMatrix;

public:
	// 将文件转化为字符串矩阵(成员变量)。按照空格和制表符分割，并且会忽略空行和注释('#'和'//'后的内容)
	// 现存问题：1.如果一行文字超过300，会被截断 2.中文乱码(getline的问题)
	void FileToWordMatrix(std::string filepath_filename) {

		std::ifstream inputFile(filepath_filename);

		if (!inputFile) {
			std::cout << "Error: no inputFile, " << __FILE__ << ", " << __LINE__ << "\n";
			return;
		}

		const int bufferLength = 300;
		char buffer[bufferLength];
		std::string tLine;
		std::vector<std::string> tWords;

		while (1) {
			inputFile.getline(buffer, bufferLength);
			tLine = buffer;
			tLine = DeleteComment(tLine);// 删除注释
			tLine = ReplaceString(tLine, '=', " = ");// '='替换为" = "
			tWords = SplitString(tLine);// 分割字符串
			// 忽略空行
			if (tWords.size() != 0) {
				tWordsMatrix.push_back(SplitString(tLine));
			}
			if (inputFile.eof())break;
		}
		inputFile.close();

		return;
	};

	// 输出字符串矩阵
	void CoutWordMatrix() {
		for (int i = 0; i < tWordsMatrix.size(); i++) {
			for (int j = 0; j < tWordsMatrix[i].size(); j++) {
				std::cout << tWordsMatrix[i][j] << ' ';
			}
			std::cout << '\n';
		}
	}

	std::vector<std::vector<std::string>> GetWordMatrix() { return tWordsMatrix; }

	// 字符串中'='替换为“ = ”
	std::string ReplaceString(std::string word, char oldChar = '=', std::string newChar = " = ") {
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

	// 删除一行字符的注释部分
	std::string DeleteComment(std::string word) {

		int i_last_meaningful_char = -1;
		for (int i = 0; i < word.size(); i++) {
			// 遇到'#'停止
			if (word[i] == '#') {
				break;
			}
			// 遇到'//'停止
			if (word[i] == '/' && i != word.size() - 1) {
				if (word[i + 1] == '/') {
					break;
				}
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

	// 分割字符串，分隔符是空格和制表符。会忽略连续的分隔符
	std::vector<std::string> SplitString(std::string tLine) {
		tLine = tLine + ' ';//目的是防止最后一个单词没有被push_back进tWords
		std::vector<std::string> tWords;
		std::string tWord;
		bool isCurrentMeaningful = 0;
		bool isLastMeaningful = 0;
		const int length = tLine.size();

		for (int j = 0; j < length; j++) {
			if (tLine[j] != ' ' && tLine[j] != '\t')
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
};

#endif // !INPUTFILEREADER_H

