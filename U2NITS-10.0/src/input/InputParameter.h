#ifndef INPUTPARAMETER_H
#define INPUTPARAMETER_H

#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <map>

// 字符串矩阵翻译为键值对
class InputParameter {

	//// 单例模式 start
private:
	static InputParameter* m_instance;
	InputParameter();
public:
	static InputParameter* getInstance() {
		if (m_instance == nullptr) {
			m_instance = new InputParameter();
		}
		return m_instance;
	}
	//// 单例模式 end
private:
	std::map<std::string, std::string> mapParameter;

public:
	std::string outputFileName = "outputfile";
	const double PI = 3.14159265368979;
	int elementNumX = 30;
	int elementNumY = 30;
	double u_wall = 1.0;// upper boundary velocity
	double Re = 100.0;
	int maxIterationStep = 36000;
	int maxIterationStep_SOR = 20000;// max iteration step in SOR
	double epsilon = 1.0e-6;// epsilon of iteration
	double epsilon_SOR = 1.0e-6;

public:

	void InitializeParameter(std::vector<std::vector<std::string>> tWordsMatrix) {
		// 字符串矩阵转map
		for (int i = 0; i < tWordsMatrix.size(); i++) {
			if (tWordsMatrix[i].size() >= 3) {
				if (tWordsMatrix[i][1] == "=") {
					mapParameter.insert(std::pair<std::string, std::string>(tWordsMatrix[i][0], tWordsMatrix[i][2]));
				}
				else {
					std::cout << "Error: invalid input, " << __FILE__ << ", " << __LINE__ << "\n";
				}
			}
		}

		getIntIfExsits(elementNumX, "elementNumX");
		getIntIfExsits(elementNumY, "elementNumY");
		getDoubleIfExsits(Re, "Re");
		getStringIfExsits(outputFileName, "outputFileName");
	}

	void getIntIfExsits(int& value_to_be_changed, std::string key = "elementNumX") {
		// 若未找到，则不改变value的值
		std::map<std::string, std::string>::iterator iter = mapParameter.find(key);
		if (iter != mapParameter.end()) {
			std::cout << key << ": " << iter->second << std::endl;
			value_to_be_changed = std::stoi(iter->second);
		}
		else {
			std::cout << "Warning: cannot find int " << key << std::endl;
		}
	}

	void getDoubleIfExsits(double& value_to_be_changed, std::string key = "Re") {
		// 若未找到，则不改变value的值
		std::map<std::string, std::string>::iterator iter = mapParameter.find(key);
		if (iter != mapParameter.end()) {
			std::cout << key << ": " << iter->second << std::endl;
			value_to_be_changed = std::stod(iter->second);
		}
		else {
			std::cout << "Warning: cannot find double " << key << std::endl;
		}
	}

	void getStringIfExsits(std::string& value_to_be_changed, std::string key = "outputFileName") {
		// 若未找到，则不改变value的值
		std::map<std::string, std::string>::iterator iter = mapParameter.find(key);
		if (iter != mapParameter.end()) {
			std::cout << key << ": " << iter->second << std::endl;
			value_to_be_changed = iter->second;
		}
		else {
			std::cout << "Warning: cannot find string " << key << std::endl;
		}
	}
};

#endif
