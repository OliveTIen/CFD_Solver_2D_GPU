#include "SU2MeshReader.h"
#include "../FVM_2D.h"
#include "../global/FilePathManager.h"
#include "../global/StringProcessor.h"
#include "../output/LogWriter.h"
#include "../global/CExit.h"
#include "../boundary_condition/BoundaryManager.h"

void SU2MeshReader::readFile_2(std::string filePath, bool convertRectToTriangle) {

	int maxnodeID = 1;
	int maxelementID = 1;
	std::vector<SimpleBoundary>tmp_boudaries;

	readMesh(filePath, convertRectToTriangle, maxnodeID, maxelementID, tmp_boudaries);
	process(maxnodeID, maxelementID, tmp_boudaries);
}

void SU2MeshReader::readMesh(std::string filePath, bool convertRectToTriangle, int& maxnodeID, int& maxelementID, std::vector<SimpleBoundary>& tmp_boudaries) {
	/*
	��ȡSU2�����ļ���
	�����ڵ��š��ڵ�����
	��Ԫ��š�ÿ����Ԫ�Ľڵ���
	*/
	LogWriter::logAndPrint("Mesh file: " + filePath + "\n");
	FVM_2D* pFVM2D = FVM_2D::getInstance();
	std::ifstream infile(filePath);
	if (!infile) {
		std::stringstream ss;
		ss << "cannot open file " << filePath << ". @SU2MeshReader::readMesh\n";
		LogWriter::logAndPrintError(ss.str());
		exit(-1);
		
	}
	State state = state_START;//1-Node, 2-Element
	const int bufferLength = 300;
	char buffer[bufferLength];
	std::string tLine;
	std::vector<std::string> tWords;
	std::vector<VirtualEdge_2D> vBoundaryEdges;//��ʱ�������洢�߽�Ԫ�ڵ�ID
	vBoundaryEdges.push_back(VirtualEdge_2D());//���1��nullԪ�أ��Ա�֤�кű�ʾID
	while (infile.getline(buffer, bufferLength)) {

		tLine = buffer;
		tLine = StringProcessor::replaceCharInString(tLine, '=', " = ");
		tWords = StringProcessor::splitString(tLine);
		long tWordsSize = tWords.size();
		// ����״̬
		if (tWordsSize == 0 || tWords[0] == "%") {
			continue;//���ڿ��У�ǿ�Ƚ�����һ��ѭ������ֹ��ȡtWords[0]�����ڴ����
		}
		if (tWords[0] == "NDIME") {
			state = state_NDIME;
		}
		else if (tWords[0] == "NPOIN") {
			state = state_NPOIN;
		}
		else if (tWords[0] == "NELEM") {
			state = state_NELEM;
			//std::cout << "Read SU2 elements\n";
		}
		else if (tWords[0] == "NMARK") {
			state = state_NMARK;
		}

		// ����״̬�����ض�����
		if (state == state_NDIME) {
			if (tWordsSize == 3 && tWords[2] == "2") {

			}
			else {
				LogWriter::logAndPrintError("SU2 mesh dimension is not 2D\n");
				CExit::saveAndExit(-1);
			}
		}
		else if (state == state_NPOIN && tWords[0].substr(0, 1) != "N") {//"*NODE" ��ȡ�ڵ�ID������
			Node_2D node;
			node.x = std::stod(tWords[0]);
			node.y = std::stod(tWords[1]);
			if (tWordsSize == 3) {// x y ID
				node.ID = (int)std::stod(tWords[2]) + 1;
			}
			else if (tWordsSize == 4) {// x y z ID
				node.ID = (int)std::stod(tWords[3]) + 1;// ��inp����һ�£���1��ʼ
			}
			pFVM2D->nodes.push_back(node);
			maxnodeID = (std::max)(maxnodeID, node.ID);
		}
		else if (state == state_NELEM && tWords[0].substr(0, 1) != "N") {
			// ����Ԫ
			if (tWords[0] == "5") {
				Element_2D e;
				e.ID = pFVM2D->elements.size() + 1;
				e.nodes[0] = (int)std::stod(tWords[1]) + 1;
				e.nodes[1] = (int)std::stod(tWords[2]) + 1;
				e.nodes[2] = (int)std::stod(tWords[3]) + 1;
				pFVM2D->elements.push_back(e);
				maxelementID = (std::max)(maxelementID, e.ID);
			}
			// �ı�Ԫ
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
				// ���һ��boundary����ʼ������
				tmp_boudaries.push_back(SimpleBoundary(tWords[2]));

			}
			else if (tWords[0] == "MARKER_ELEMS") {

			}
			else {
				// ��ǰboundary���һ����
				SimpleBoundary& tmp_currentBoundary = tmp_boudaries[tmp_boudaries.size() - 1];
				tmp_currentBoundary.edges.push_back(SimpleEdge(
					std::stoi(tWords[1]) + 1, std::stoi(tWords[2]) + 1
				));
			}
		}
	}
	infile.close();

}

void SU2MeshReader::process(int maxNodeID, int maxElementID, std::vector<SimpleBoundary>& tmp_boudaries) {
	FVM_2D* pFVM2D = FVM_2D::getInstance();

	std::stringstream ss;
	ss << "Initiate pNodeTable, pElementTable, pEdgeTable.\b";
	ss << pFVM2D->elements.size() << " elements, " << pFVM2D->nodes.size() << " nodes.\n";
	LogWriter::logAndPrint(ss.str());
	pFVM2D->iniPNodeTable(maxNodeID);
	pFVM2D->iniPElementTable(maxElementID);
	pFVM2D->iniEdges();
	pFVM2D->iniPEdgeTable();
	pFVM2D->iniNode_neighborElements();//ǰ������ elements, elements.nodes, pNodeTable

	LogWriter::logAndPrint("Calculate Element xy, edge length\n");
	pFVM2D->iniElement_xy_pEdges();
	pFVM2D->iniEdges_lengths();

	LogWriter::logAndPrint("Initialize boundary condition\n");
	// ��ʼ��boundaryManager.boundaries��pEdges��type
	// ��ʼ��edge��setID
	// ��ʼ�����ڱ߽��ϵ����boundaryManager.periodPairsά��
	// ǰ����������boundaries��������edges������edges�ѳ�ʼ��nodeIDs��boundaries��name
	/*
	std::vector<SimpleBoundary>& tmp_boudaries����ȡsu2 meshʱ������SimpleBoundary���б߽����ƣ�����"periodic_2"���ͱߵ����жԣ�����{51,52}��
	real_boundaries����pFVM2D->boundaryManager.boundaries�����ã��ڴ˴�����ʼ��
	*/
	{
		std::vector<VirtualBoundarySet_2D>& real_boundaries = pFVM2D->boundaryManager.boundaries;
		real_boundaries.resize(tmp_boudaries.size());
		for (auto i = 0; i < tmp_boudaries.size(); i++) {
			real_boundaries[i].name = tmp_boudaries[i].name;
			real_boundaries[i].type = U2NITS::BoundaryManager::boundaryNameToType(real_boundaries[i].name);
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

			//! ע��öδ�������ڶദ���޸�ʱ������Ҫ�޸Ķദ
			// ��ʼ��periodPairs��periodPairs�洢���ڱ߽�������Ϣ������������ڱ߽�������
			// periodPairs��boundaryManager�ĳ�Ա��ÿ��pair�洢int bType, int setID_0, int setID_1
			const int& bType = real_boundaries[i].type;
			const int& bID = real_boundaries[i].ID;
			auto& periodPairs = pFVM2D->boundaryManager.periodPairs;
			if (_BC_periodic_0 <= bType && bType <= _BC_periodic_9) {// ��������ʽ��𿪣���&&����
				//���periodPairs�Ƿ��Ѿ�¼���bType�����ǣ���ʹ�����еģ��������½�һ��PeriodPair tmp������periodPairs
				int index_pairs = -1;
				for (int i = 0; i < periodPairs.size(); i++) {
					if (periodPairs[i].bType == bType)index_pairs = i;
				}
				// ���ڣ���ֱ������
				if (index_pairs != -1) {
					periodPairs[index_pairs].setID_1 = bID;
				}
				// �����ڣ����½�
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
	// ������ڱ߽���ȷ��
	pFVM2D->boundaryManager.checkPeriodPairs();

	//LogWriter::log("elements.size="+std::to_string())

}
