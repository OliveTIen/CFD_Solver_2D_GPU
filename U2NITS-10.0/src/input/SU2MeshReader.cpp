#include "SU2MeshReader.h"
#include "../FVM_2D.h"
#include "../global/FilePathManager.h"
#include "../global/StringProcessor.h"
//#include "../global/VectorProcessor.h"


int SU2MeshReader::readFile(std::string filePath, bool convertRectToTriangle) {
	/*
	��Ҫ��ʼ����������
	1. nodes���������ꡢID��neighboringElements
	2. elements������ID����Ԫ���͡��ڵ�ID�����⣬Ϊ���ٺ������㣬��Ҫ��ǰ���㵥Ԫ��������
	3. edges������ID���ڵ�ID�������߼���ID�����ҵ�Ԫָ�롣���⣬Ҳ��Ҫ��ǰ����length��refLength

	������Ǳ������
	1. �е�su2�ļ�����%дע�͡�Ŀǰû������

	SU2�ļ���ȡ�淶 https://zhuanlan.zhihu.com/p/641146110
	*/

	std::cout << "Read \".su2\" mesh file\n";
	FVM_2D* pFVM2D = FVM_2D::pFVM2D;
	std::ifstream infile(filePath);
	if (!infile) {
		return ERROR_READ_FILE;
	}
	int state = -1;//1-Node, 2-Element
	int maxnodeID = 1;
	int maxelementID = 1;
	const int bufferLength = 300;
	char buffer[bufferLength];
	std::string tLine;
	std::vector<std::string> tWords;
	std::vector<VirtualEdge_2D> vBoundaryEdges;//��ʱ�������洢�߽�Ԫ�ڵ�ID
	vBoundaryEdges.push_back(VirtualEdge_2D());//���1��nullԪ�أ��Ա�֤�кű�ʾID
	std::vector<SimpleBoundary>tmp_boudaries;
	while (infile.getline(buffer, bufferLength)) {
		//get words and set state
		//infile.getline(buffer, bufferLength);
		//if (infile.eof())break;
		tLine = buffer;
		tLine = StringProcessor::replaceCharInString(tLine, '=', " = ");
		tWords = StringProcessor::splitString(tLine);
		if (tWords.size() == 0)continue;//���ڿ��У�ǿ�Ƚ�����һ��ѭ������ֹ��ȡtWords[0]�����ڴ����
		else {
			// ά��
			if (tWords[0] == "NDIME")state = 0;
			// ��
			if (tWords[0] == "NPOIN")state = 1;
			// ��Ԫ
			if (tWords[0] == "NELEM") {
				state = 2;
			}
			// �߽�
			if (tWords[0] == "NMARK") {
				state = 3;
			}
			
		}

		//translate
		if (state == 1 && tWords[0].substr(0, 1) != "N") {//"*NODE" ��ȡ�ڵ�ID������
			Node_2D node;
			node.x = std::stod(tWords[0]);
			node.y = std::stod(tWords[1]);
			node.ID = (int)std::stod(tWords[3]) + 1;// ��inp����һ�£���1��ʼ
			pFVM2D->nodes.push_back(node);
			maxnodeID = (std::max)(maxnodeID, node.ID);
		}
		else if (state == 2 && tWords[0].substr(0, 1) != "N") {
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
		else if (state == 3) {//"NMARK"
			if (tWords[0] == "NMARK") {

			}
			else if (tWords[0] == "MARKER_TAG") {
				// ����һ��boundary����ʼ������
				tmp_boudaries.push_back(SimpleBoundary(tWords[2]));

			}
			else if(tWords[0] == "MARKER_ELEMS"){

			}
			else {
				// ��ǰboundary����һ����
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
	pFVM2D->iniNode_neighborElements();//ǰ������ elements, elements.nodes, pNodeTable

	std::cout << "Calculate Element xy, edge length" << std::endl;
	pFVM2D->iniElement_xy_pEdges();
	pFVM2D->iniEdges_lengths();

	std::cout << "Initialize boundary condition\n";
	// ��ʼ��boundaryManager.boundaries��pEdges��type
	// ��ʼ��edge��setID
	// ��ʼ�����ڱ߽��ϵ����boundaryManager.periodPairsά��
	// ǰ����������boundaries��������edges������edges�ѳ�ʼ��nodeIDs��boundaries��name
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
	// ������ڱ߽���ȷ��
	pFVM2D->boundaryManager.checkPeriodPairs();


    return 0;
}