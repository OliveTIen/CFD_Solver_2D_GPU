#pragma once
#include <string>
#include <vector>
#include <glm/glm.hpp>
#include "typeDefine.h"

namespace GUI {
	class UGridData;

	class UGridReader {
	public:

		static UGridData* getData();
		static void readData(std::string filePath, bool convertRectToTriangle);
		static void m_readSU2File(std::string filePath, bool convertRectToTriangle, int& maxNodeID, int& maxElementID);
		static void m_scaleNodeCoordiate(double factor);
		
	private:
		static std::vector<std::string> splitString(std::string tLine);
		static std::string replaceCharInString(std::string word, char oldChar, std::string newChar);
		static void logConsole(std::string content);

	};

	///////////////////////////////////////////

	class UGridData {
	public:
		class Node {
		public:
			int ID = -1;
			myfloat x = 0;
			myfloat y = 0;
		};

		class Element {
		public:
			int ID = -1;
			int nodes[4]{ -1,-1,-1,-1 };
		};

		std::vector<Node> nodes;
		std::vector<Element> elements;
		std::vector<std::pair<std::string, std::vector<glm::uvec2>>>boundaries;
		static UGridData* getInstance();

	private:
		UGridData() {};
		static UGridData* pData;
	};
}