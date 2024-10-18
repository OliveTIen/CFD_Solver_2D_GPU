#include <fstream>
#include "UGridReader.h"

using namespace GUI;

UGridData* UGridData::pData = nullptr;
UGridData* UGridData::getInstance() {
	if (pData == nullptr) {
		pData = new UGridData();
	}
	return pData;
}


std::vector<std::string> UGridReader::splitString(std::string tLine) {
	/*
	split string
	*/
	tLine = tLine + ' ';//目的是防止最后一个单词没有被push_back进tWords
	std::vector<std::string> tWords;
	std::string tWord;
	bool isCurrentMeaningful = 0;
	bool isLastMeaningful = 0;
	const int length = (int)tLine.size();
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

std::string UGridReader::replaceCharInString(std::string word, char oldChar, std::string newChar) {
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

/////////////////////////////////////////

UGridData* GUI::UGridReader::getData() {
	return UGridData::getInstance();
}

void UGridReader::readData(std::string filePath, bool convertRectToTriangle) {

	int maxnodeID = 1;
	int maxelementID = 1;
	double mesh_scale_factor = 1.0;

	m_readSU2File(filePath, convertRectToTriangle, maxnodeID, maxelementID);
	m_scaleNodeCoordiate(mesh_scale_factor);
}


void UGridReader::m_readSU2File(std::string filePath, bool convertRectToTriangle, int& maxnodeID, int& maxelementID) {
	/*
	Updated on August 13th, 2024
	
	Function Introduction:
	This functions reads SU2 mesh file, 
	including node indices, node coordinates, 
	element indices, and element-node relations
	
	*/
	// 状态。读取网格时，读取到某一位置，会有特定状态
	enum State {
		state_START,
		state_NDIME,
		state_NPOIN,
		state_NPOIN_data,
		state_NELEM,
		state_NELEM_data,
		state_NMARK,
		state_NMARK_data,
	};

	int offset = 0;// offset = 0 means index starts from zero
	logConsole("Mesh file: " + filePath + "\n");
	UGridData* pData = getData();
	auto& boundaries = pData->boundaries;
	std::ifstream infile(filePath);
	if (!infile) {
		logConsole("cannot open file " + filePath + " while reading su2 mesh\n");
		exit(-1);
		
	}
	State state = state_START;//1-Node, 2-Element
	const int bufferLength = 300;
	char buffer[bufferLength];
	std::string tLine;
	std::vector<std::string> tWords;
	while (infile.getline(buffer, bufferLength)) {

		tLine = buffer;
		tLine = replaceCharInString(tLine, '=', " = ");
		tWords = splitString(tLine);
		long tWordsSize = tWords.size();
		// update state machine
		if (tWordsSize == 0 || tWords[0] == "%") {
			continue;// skip to avoid illegal memory access (cannot read tWords[0] when tWord is empty)
		}
		if (tWords[0] == "NDIME") {
			state = state_NDIME;
		}
		else if (tWords[0] == "NPOIN") {
			state = state_NPOIN;
		}
		else if (tWords[0] == "NELEM") {
			state = state_NELEM;
		}
		else if (tWords[0] == "NMARK") {
			state = state_NMARK;
		}

		// operate according to states
		if (state == state_NDIME) {
			if (tWordsSize == 3 && tWords[2] == "2") {

			}
			else {
				logConsole("SU2 mesh dimension is not 2D\n");
				exit(-1);
			}
		}
		else if (state == state_NPOIN && tWords[0].substr(0, 1) != "N") {// meets "*NODE" in file
			UGridData::Node node;
			node.x = std::stod(tWords[0]);
			node.y = std::stod(tWords[1]);
			if (tWordsSize == 3) {// x y ID
				node.ID = (int)std::stod(tWords[2]) + offset;
			}
			else if (tWordsSize == 4) {// x y z ID
				node.ID = (int)std::stod(tWords[3]) + offset;
			}
			pData->nodes.push_back(node);
			maxnodeID = (std::max)(maxnodeID, node.ID);
		}
		else if (state == state_NELEM && tWords[0].substr(0, 1) != "N") {
			// triangle elements
			if (tWords[0] == "5") {
				UGridData::Element e;
				e.ID = pData->elements.size() + offset;
				e.nodes[0] = (int)std::stod(tWords[1]) + offset;
				e.nodes[1] = (int)std::stod(tWords[2]) + offset;
				e.nodes[2] = (int)std::stod(tWords[3]) + offset;
				pData->elements.push_back(e);
				maxelementID = (std::max)(maxelementID, e.ID);
			}
			// quads
			else if (tWords[0] == "9") {
				if (convertRectToTriangle) {
					UGridData::Element e;
					e.ID = pData->elements.size() + offset;
					e.nodes[0] = (int)std::stod(tWords[1]) + offset;
					e.nodes[1] = (int)std::stod(tWords[2]) + offset;
					e.nodes[2] = (int)std::stod(tWords[4]) + offset;
					pData->elements.push_back(e);
					maxelementID = (std::max)(maxelementID, e.ID);

					e.ID = pData->elements.size() + offset;
					e.nodes[0] = (int)std::stod(tWords[4]) + offset;
					e.nodes[1] = (int)std::stod(tWords[2]) + offset;
					e.nodes[2] = (int)std::stod(tWords[3]) + offset;
					pData->elements.push_back(e);
					maxelementID = (std::max)(maxelementID, e.ID);
				}
				else {
					UGridData::Element e;
					e.ID = pData->elements.size() + offset;
					e.nodes[0] = (int)std::stod(tWords[1]) + offset;
					e.nodes[1] = (int)std::stod(tWords[2]) + offset;
					e.nodes[2] = (int)std::stod(tWords[3]) + offset;
					e.nodes[3] = (int)std::stod(tWords[4]) + offset;
					pData->elements.push_back(e);
					maxelementID = (std::max)(maxelementID, e.ID);
				}
			}
		}
		else if (state == state_NMARK) {//meets "NMARK" in file
			if (tWords[0] == "NMARK") {

			}
			else if (tWords[0] == "MARKER_TAG") {
				// boundary is made up of name and vertices
				// add a boundary and set its name with tWords[2]
				auto boundary = std::pair<std::string, std::vector<glm::uvec2>>(tWords[2], 0);
				boundaries.push_back(boundary);

			}
			else if (tWords[0] == "MARKER_ELEMS") {
				// do nothing
			}
			else {
				// add an edge to the current boundary
				int vertex0 = std::stoi(tWords[1]) + offset;
				int vertex1 = std::stoi(tWords[2]) + offset;
				auto& boundary = boundaries[boundaries.size() - 1];
				boundary.second.push_back(glm::uvec2(vertex0, vertex1));
			}
		}
	}
	infile.close();

}

void UGridReader::m_scaleNodeCoordiate(double factor) {
	if (factor < 0.0) {
		logConsole("Mesh scale factor should be positive. Function doesn't take effect\n");
		return;
	}
	if (abs(factor - 1.0) < 0.01) {
		//logConsole("Mesh scale factor is too close to 1.0. Function doesn't take effect\n");
		return;
	}
	for (auto& node : getData()->nodes) {
		node.x *= factor;
		node.y *= factor;
	}
}

void GUI::UGridReader::logConsole(std::string content) {
	printf(content.c_str());
}
