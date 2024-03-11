#include "SU2MeshReader.h"
#include "../FVM_2D.h"
#include "../global/FilePathManager.h"
#include "../global/StringProcessor.h"
#include "../output/LogWriter.h"
//#include "../global/VectorProcessor.h"


int SU2MeshReader::readFile(std::string filePath, bool convertRectToTriangle) {
	/*
	需要初始化以下内容
	1. nodes。包括坐标、ID、neighboringElements
	2. elements。包括ID、单元类型、节点ID。此外，为减少后续计算，需要提前计算单元中心坐标
	3. edges。包括ID、节点ID、所属边集合ID、左右单元指针。此外，也需要提前计算length和refLength

	其他的潜在问题
	1. 有的su2文件会用%写注释。目前没遇到过

	SU2文件读取规范 https://zhuanlan.zhihu.com/p/641146110
	*/

	std::cout << "Read \".su2\" mesh file\n";
	FVM_2D* pFVM2D = FVM_2D::getInstance();
	std::ifstream infile(filePath);
	if (!infile) {
		return ERROR_READ_FILE;
	}
	State state = state_START;//1-Node, 2-Element
	int maxnodeID = 1;
	int maxelementID = 1;
	const int bufferLength = 300;
	char buffer[bufferLength];
	std::string tLine;
	std::vector<std::string> tWords;
	std::vector<VirtualEdge_2D> vBoundaryEdges;//临时变量。存储边界元节点ID
	vBoundaryEdges.push_back(VirtualEdge_2D());//填充1个null元素，以保证行号表示ID
	std::vector<SimpleBoundary>tmp_boudaries;
	while (infile.getline(buffer, bufferLength)) {

		tLine = buffer;
		tLine = StringProcessor::replaceCharInString(tLine, '=', " = ");
		tWords = StringProcessor::splitString(tLine);
		long tWordsSize = tWords.size();
		// 更新状态
		if (tWordsSize == 0) {
			continue;//对于空行，强迫进入下一次循环，防止读取tWords[0]出现内存错误
		}
		if (tWords[0] == "NDIME") {
			state = state_NDIME;
		}
		else if (tWords[0] == "NPOIN") {
			state = state_NPOIN;
		}
		else if (tWords[0] == "NELEM") {
			state = state_NELEM;
			std::cout << "Read SU2 elements\n";
		}
		else if (tWords[0] == "NMARK") {
			state = state_NMARK;
		}
			
		// 根据状态进行特定操作
		if (state == state_NDIME) {
			if (tWordsSize == 3 && tWords[2] == "2") {

			}
			else {
				LogWriter::writeLogAndCout("Error: SU2 mesh dimension is not 2D\n", LogWriter::Error);
				exit(-1);
			}
		}
		else if (state == state_NPOIN && tWords[0].substr(0, 1) != "N") {//"*NODE" 读取节点ID和坐标
			Node_2D node;
			node.x = std::stod(tWords[0]);
			node.y = std::stod(tWords[1]);
			if (tWordsSize == 3) {// x y ID
				node.ID = (int)std::stod(tWords[2]) + 1;
			}
			else if (tWordsSize == 4) {// x y z ID
				node.ID = (int)std::stod(tWords[3]) + 1;// 跟inp保持一致，从1开始
			}
			pFVM2D->nodes.push_back(node);
			maxnodeID = (std::max)(maxnodeID, node.ID);
		}
		else if (state == state_NELEM && tWords[0].substr(0, 1) != "N") {
			// 三角元
			if (tWords[0] == "5") {
				Element_2D e;
				e.ID = pFVM2D->elements.size() + 1;
				e.nodes[0] = (int)std::stod(tWords[1]) + 1;
				e.nodes[1] = (int)std::stod(tWords[2]) + 1;
				e.nodes[2] = (int)std::stod(tWords[3]) + 1;
				pFVM2D->elements.push_back(e);
				maxelementID = (std::max)(maxelementID, e.ID);
			}
			// 四边元
			else if (tWords[0] == "9") {
				if (convertRectToTriangle) {
					Element_2D e;
					e.ID = pFVM2D->elements.size() + 1;
					e.nodes[0] = (int)std::stod(tWords[1]) + 1;
					e.nodes[1] = (int)std::stod(tWords[2]) + 1;
					e.nodes[2] = (int)std::stod(tWords[4]) + 1;
					pFVM2D->elements.push_back(e);
					maxelementID = (std::max)(maxelementID, e.ID);

					e.ID = pFVM2D->elements.size() + 1;
					e.nodes[0] = (int)std::stod(tWords[4]) + 1;
					e.nodes[1] = (int)std::stod(tWords[2]) + 1;
					e.nodes[2] = (int)std::stod(tWords[3]) + 1;
					pFVM2D->elements.push_back(e);
					maxelementID = (std::max)(maxelementID, e.ID);
				}
				else {
					Element_2D e;
					e.ID = pFVM2D->elements.size() + 1;
					e.nodes[0] = (int)std::stod(tWords[1]) + 1;
					e.nodes[1] = (int)std::stod(tWords[2]) + 1;
					e.nodes[2] = (int)std::stod(tWords[3]) + 1;
					e.nodes[3] = (int)std::stod(tWords[4]) + 1;
					pFVM2D->elements.push_back(e);
					maxelementID = (std::max)(maxelementID, e.ID);
				}
			}
		}
		else if (state == state_NMARK) {//"NMARK"
			if (tWords[0] == "NMARK") {

			}
			else if (tWords[0] == "MARKER_TAG") {
				// 添加一个boundary，初始化名称
				tmp_boudaries.push_back(SimpleBoundary(tWords[2]));

			}
			else if(tWords[0] == "MARKER_ELEMS"){

			}
			else {
				// 当前boundary添加一个边
				SimpleBoundary& tmp_currentBoundary = tmp_boudaries[tmp_boudaries.size() - 1];
				tmp_currentBoundary.edges.push_back(SimpleEdge(
					std::stoi(tWords[1]) + 1, std::stoi(tWords[2]) + 1
				));
			}
		}
	}
	infile.close();


	pFVM2D->iniEdges();
	std::cout << "Initiate pNodeTable, pElementTable, pEdgeTable" << std::endl;
	pFVM2D->iniPNodeTable(maxnodeID);
	pFVM2D->iniPElementTable(maxelementID);
	pFVM2D->iniPEdgeTable();
	pFVM2D->iniNode_neighborElements();//前置条件 elements, elements.nodes, pNodeTable

	std::cout << "Calculate Element xy, edge length" << std::endl;
	pFVM2D->iniElement_xy_pEdges();
	pFVM2D->iniEdges_lengths();

	std::cout << "Initialize boundary condition\n";
	// 初始化boundaryManager.boundaries的pEdges、type
	// 初始化edge的setID
	// 初始化周期边界关系，由boundaryManager.periodPairs维护
	// 前置条件：有boundaries向量，有edges向量且edges已初始化nodeIDs，boundaries有name
	{
		std::vector<VirtualBoundarySet_2D>& real_boundaries = pFVM2D->boundaryManager.boundaries;
		real_boundaries.resize(tmp_boudaries.size());
		for (auto i = 0; i < tmp_boudaries.size(); i++) {
			real_boundaries[i].name = tmp_boudaries[i].name;
			real_boundaries[i].type = BoundaryManager_2D::getBoundaryTypeByName(real_boundaries[i].name);
			real_boundaries[i].ID = i + 1;
			for (auto& edge : tmp_boudaries[i].edges) {
				Edge_2D* pEdge = pFVM2D->getEdgeByNodeIDs(edge.nodeIDs[0], edge.nodeIDs[1]);
				if (pEdge == nullptr) {
					throw "Error: pEdge == nullptr, SU2MeshReader::readFile()";
				}
				else {
					pEdge->setID = real_boundaries[i].ID;
				}
				real_boundaries[i].pEdges.push_back(pEdge);
			}

			//! 注意该段代码存在于多处。修改时可能需要修改多处
			// 初始化periodPairs，periodPairs存储周期边界的配对信息，用来检查周期边界完整性
			// periodPairs是boundaryManager的成员，每个pair存储int bType, int setID_0, int setID_1
			const int& bType = real_boundaries[i].type;
			const int& bID = real_boundaries[i].ID;
			auto& periodPairs = pFVM2D->boundaryManager.periodPairs;
			if (_BC_periodic_0 <= bType && bType <= _BC_periodic_9) {// 两个不等式需拆开，用&&连接
				//检查periodPairs是否已经录入该bType，若是，则使用已有的；若否，则新建一个PeriodPair tmp，存入periodPairs
				int index_pairs = -1;
				for (int i = 0; i < periodPairs.size(); i++) {
					if (periodPairs[i].bType == bType)index_pairs = i;
				}
				// 存在，则直接设置
				if (index_pairs != -1) {
					periodPairs[index_pairs].setID_1 = bID;
				}
				// 不存在，则新建
				else {
					BoundaryManager_2D::PeriodPair tmp;
					tmp.bType = bType;
					tmp.setID_0 = bID;
					periodPairs.push_back(tmp);
				}
			}
		}
	}
	//pFVM2D->boundaryManager.iniBoundaryEdgeSetID_and_iniBoundaryType(pFVM2D);
	// 检查周期边界正确性
	pFVM2D->boundaryManager.checkPeriodPairs();


    return 0;
}
