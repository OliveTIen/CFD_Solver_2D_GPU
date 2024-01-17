#include "GmeshMeshReader.h"
#include <iostream>
#include <fstream>
#include "../global/FilePathManager.h"
#include "../GlobalPara.h"
#include "../output/LogWriter.h"
#include "../BoundaryManager_2D.h"
#include "../global/StringProcessor.h"
#include "../FVM_2D.h"

int GmeshMeshReader::readMeshFile(std::string suffix) {
	FVM_2D* pFVM2D = FVM_2D::pFVM2D;
	std::cout << "read mesh file\n";
	std::string dir = FilePathManager::getInstance()->getExePath_withSlash() + "input\\";
	std::ifstream infile(dir + GlobalPara::basic::filename + suffix);
	if (!infile) {
		LogWriter::writeLogAndCout("Warning: fail to read mesh(*.inp). (FVM_2D::readMeshFile) \n");
		return -1;
	}
	int state = 0;//1-Node, 2-Element
	int substate = 0;//1-Line
	int maxnodeID = 1;
	int maxedgeID = 1;
	int maxelementID = 1;
	const int bufferLength = 300;
	char buffer[bufferLength];
	std::string tLine;
	std::vector<std::string> tWords;
	std::vector<std::vector<std::string>> tWordsMatrix;
	std::vector<VirtualEdge_2D> vBoundaryEdges;//临时变量。存储边界元节点ID
	vBoundaryEdges.push_back(VirtualEdge_2D());//填充1个null元素，以保证行号表示ID
	std::vector<int>vInts;//临时变量
	std::string setName;//临时变量
	while (1) {
		//get words and set state
		infile.getline(buffer, bufferLength);
		if (infile.eof())break;
		tLine = buffer;
		tWords = StringProcessor::splitString(tLine);
		if (tWords.size() == 0)continue;//对于空行，强迫进入下一次循环，防止读取tWords[0]出现内存错误
		else {
			if (tWords[0] == "*NODE")state = 1;
			if (tWords[0] == "*ELEMENT") {
				state = 2;
				if (tWords[2].substr(0, 7) == "ELSET=L") substate = 1;//ELSET=Line[x]
				else if (tWords[2].substr(0, 7) == "ELSET=S") substate = 2;
			}
			if (tWords[0] == "*ELSET") {
				state = 3;
			}
			if (tWords[0] == "*NSET") state = 4;//目前没用到。因此gmsh导出inp时可以不导出节点
		}

		//translate//注意慎用stoi
		if (state == 1 && tWords[0].substr(0, 1) != "*") {//"*NODE" 读取节点ID和坐标
			Node_2D node;
			node.ID = (int)std::stod(tWords[0]);
			node.x = std::stod(tWords[1]);
			node.y = std::stod(tWords[2]);
			pFVM2D->nodes.push_back(node);
			maxnodeID = (std::max)(maxnodeID, node.ID);
		}
		else if (state == 2 && substate == 1 && tWords[0].substr(0, 1) != "*") {//边界元 "*ELEMENT""ELSET=Line"
			VirtualEdge_2D ve;//成员：nodes[2]
			ve.nodes[0] = std::stoi(tWords[1]);//1, 1, 6; EdgeID,node0,node1
			ve.nodes[1] = std::stoi(tWords[2]);
			vBoundaryEdges.push_back(ve);
		}
		else if (state == 2 && substate == 2 && tWords[0].substr(0, 1) != "*") {//面单元 "*ELEMENT""ELSET=Surface"
			Element_T3 e;
			e.ID = (int)std::stod(tWords[0]);
			e.nodes[0] = (int)std::stod(tWords[1]);
			e.nodes[1] = (int)std::stod(tWords[2]);
			e.nodes[2] = (int)std::stod(tWords[3]);
			pFVM2D->elements.push_back(e);
			maxelementID = (std::max)(maxelementID, e.ID);
		}
		else if (state == 3) {//"*ELSET"
			if (tWords[0].substr(0, 1) == "*") {
				//初始化vBoundarySets
				std::vector<int> start_end = pFVM2D->boundaryManager.splitInts(vInts);
				int nSets = (int)start_end.size() / 2;
				for (int is = 0; is < nSets; is++) {//暗含nSets!=0，说明setName已经被初始化
					VirtualBoundarySet_2D vb;
					vb.name = setName;
					vb.ID = (int)pFVM2D->boundaryManager.vBoundarySets.size() + 1;
					vb.startID = start_end[is * 2 + 0];//0,2,4,... 0=is*2+0
					vb.endID = start_end[is * 2 + 1];
					pFVM2D->boundaryManager.vBoundarySets.push_back(vb);
				}
				//重置
				vInts.clear();
				setName = tWords[1].substr(6);//"wall""outlet""inlet"...
			}
			else {
				for (int iword = 0; iword < tWords.size(); iword++) {
					vInts.push_back(std::stoi(tWords[iword]));
				}
			}
		}
	}
	infile.close();

	std::cout << "Generate and initiate elements..." << std::endl;

	pFVM2D->iniPNodeTable(maxnodeID);
	pFVM2D->iniEdges();
	pFVM2D->iniPEdgeTable();
	pFVM2D->iniPElementTable(maxelementID);
	pFVM2D->iniElement_xy_pEdges();
	pFVM2D->iniNode_neighborElements();
	pFVM2D->iniEdges_lengths();
	pFVM2D->boundaryManager.iniBoundarySetPEdges_in_readMeshFile(pFVM2D, vBoundaryEdges);//初始化vBoundarySet的pEdges


	return 0;

}
