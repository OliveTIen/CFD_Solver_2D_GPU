#include "ContinueFileReader.h"
#include <fstream>
#include <iostream>
#include "../legacy/FVM_2D.h"
#include "../global/FilePathManager.h"
#include "../global/GlobalPara.h"
#include "../global/StringProcessor.h"
#include "../output/LogWriter.h"
#include "../drivers/RoutineController.h"
#include "SU2MeshReader.h"

static void jointVectors(std::vector<int>& accepter, const std::vector<int>& giver) {
	// ƴ������
	for (int i = 0; i < giver.size(); i++) {
		accepter.push_back(giver[i]);
	}
}

int ContinueFileReader::readContinueFile_1() {
	/*
	Ѱ�ҵ�ǰĿ¼���Ƿ�����Ϊ��pause_�����ļ������У����ȡ�����ޣ��򷵻�-1
	���뷵��һ��ֵ���������-1�����Դ�ͷ��ʼ
	*/

	FVM_2D* pFVM2D = FVM_2D::getInstance();

	std::string m_path = FilePathManager::getInstance()->getOutputDirectory();
	std::vector<std::string> files = FilePathManager::getInstance()->getFiles(m_path);//filesΪ�ַ����������洢��·���������ļ���
	int index_maxNstep = -1;
	//��Ѱ�Ƿ���pause_filename[xxx].dat�ļ������У�����xxx����
	int maxNstep = 0;
	int tmp_index = 0;//����ָʾ"[", "]"
	std::string tmp_str;
	for (int i = 0; i < files.size(); i++) {
		//��Ѱ��pause_filename��ͷ���ļ���
		std::string str_match = "pause_" + GlobalPara::basic::filename;
		if (files[i].substr(0, int(str_match.size())) == str_match) {
			//�޳�"pause_filename["
			int pre_length = int(str_match.size() + 1);//"pause_filename["�ĳ���
			std::string str = files[i].substr(pre_length);
			//�޳�"].dat"
			int post_length = 5;//"].dat"�ĳ���
			for (int i = 0; i < post_length; i++) {
				str.pop_back();
			}
			//��xxx����num
			int num = std::stoi(str);
			if (num > maxNstep) {
				index_maxNstep = i;
				maxNstep = num;
			}
		}
	}
	//���ޣ��򷵻�-1
	if (index_maxNstep == -1) {
		return -1;
	}
	std::ifstream infile(m_path + files[index_maxNstep]);
	if (!infile) {
		return -1;
	}

	LogWriter::logAndPrint("continue from: " + files[index_maxNstep] + "\n");

	int state = -1;//1-Node, 2-Element
	int maxnodeID = 1;
	int maxedgeID = 1;
	int maxelementID = 1;
	int iline_set = 0;//��ʱ���������ڼ�¼��ȡsetʱ�����˶�����
	std::vector<std::vector<int>> edges_of_all_sets;//��ʱ�������洢��set��edgeID
	const int bufferLength = 600;//!!!С�ĳ��Ȳ���
	char buffer[bufferLength];
	std::string tLine;
	std::vector<std::string> tWords;
	std::vector<std::vector<std::string>> tWordsMatrix;
	std::vector<int> edge_ID_long;
	int nNullRow = 0;
	bool loop = 1;
	while (loop) {
		//get words and set state
		infile.getline(buffer, bufferLength);
		if (infile.eof())loop = 0;
		tLine = buffer;
		tWords = StringProcessor::splitString(tLine);
		if (nNullRow >= 10)loop = 0;
		if (tWords.size() == 0) {
			nNullRow++;
			continue;
		}//���У�ǿ�ȿ�ʼ��һ��ѭ������ֹtWords[0]�ڴ����
		else {
			if (tWords[0] == "t")state = 0;
			if (tWords[0] == "nodes:")state = 1;
			if (tWords[0] == "edges:")state = 2;
			if (tWords[0] == "elements:")state = 3;
			if (tWords[0] == "boundaries:")state = 4;
		}

		//translate
		if (state == 0 && tWords[0].substr(0, 1) != "t") {//t_solve, istep_solve, CFL
			GlobalPara::time::t_previous = std::stod(tWords[0]);
			GlobalPara::time::istep_previous = (int)std::stod(tWords[1]);
			if (tWords.size() == 3) {
				GlobalPara::time::CFL = std::stod(tWords[2]);// ��ȡ�ϴα����CFL
			}
		}
		else if (state == 1 && tWords[0].substr(0, 1) != "n") {//nodes: ID, x, y
			Node_2D node;
			node.ID = (int)std::stod(tWords[0]);
			node.x = std::stod(tWords[1]);
			node.y = std::stod(tWords[2]);
			pFVM2D->nodes.push_back(node);
			maxnodeID = (std::max)(maxnodeID, node.ID);
		}
		else if (state == 3 && tWords[0].substr(0, 1) != "e") {//elements: ID, nodes, U
			Element_2D ele;
			ele.ID = (int)std::stod(tWords[0]);
			ele.nodes[0] = (int)std::stod(tWords[1]);
			ele.nodes[1] = (int)std::stod(tWords[2]);
			ele.nodes[2] = (int)std::stod(tWords[3]);
			ele.U[0] = std::stod(tWords[4]);
			ele.U[1] = std::stod(tWords[5]);
			ele.U[2] = std::stod(tWords[6]);
			ele.U[3] = std::stod(tWords[7]);
			pFVM2D->elements.push_back(ele);
			maxelementID = (std::max)(maxelementID, ele.ID);
		}
		else if (state == 4 && tWords[0].substr(0, 1) != "b") {//boundaries: ID, name, edgeIDs
			iline_set++;
			if (tWords.size() >= 2 && !isdigit(tWords[1][0])) {
				/*
				"1 inf"
				1.vBoundarySets������Ŀ������Ŀ�Ĳ��ֳ�ʼ��
				2.edges_of_all_sets������Ŀ
				*/
				//set��ID,name��ʼ��
				VirtualBoundarySet_2D vb;
				vb.ID = (int)std::stod(tWords[0]);
				vb.name = tWords[1];
				pFVM2D->boundaryManager.boundaries.push_back(vb);
				//edges_of_all_sets������Ŀ
				edges_of_all_sets.push_back(std::vector<int>());
			}
			else {
				// 1.edges_of_all_sets�ĳ�ʼ��
				// �������һ��BoundarySet
				std::vector<int>& edges_of_current_set = edges_of_all_sets[edges_of_all_sets.size() - 1];
				std::vector<int> intVector = StringProcessor::stringVector_2_intVector(tWords);
				jointVectors(edges_of_current_set, intVector);
			}
		}
	}

	infile.close();
	process(maxnodeID, maxelementID, edges_of_all_sets);

	return 0;

}

int ContinueFileReader::readContinueFile_2_unused_addUxUy() {
	/*
	Ѱ�ҵ�ǰĿ¼���Ƿ�����Ϊ��pause_�����ļ������У����ȡ�����ޣ��򷵻�-1
	������־��20240405�����Ux Uy
	*/

	FVM_2D* pFVM2D = FVM_2D::getInstance();

	std::string m_path = FilePathManager::getInstance()->getOutputDirectory();
	std::vector<std::string> files = FilePathManager::getInstance()->getFiles(m_path);//filesΪ�ַ����������洢��·���������ļ���
	int index_maxNstep = -1;
	//��Ѱ�Ƿ���pause_filename[xxx].dat�ļ������У�����xxx����
	int maxNstep = 0;
	int tmp_index = 0;//����ָʾ"[", "]"
	std::string tmp_str;
	for (int i = 0; i < files.size(); i++) {
		//��Ѱ��pause_filename��ͷ���ļ���
		std::string str_match = "pause_" + GlobalPara::basic::filename;
		if (files[i].substr(0, int(str_match.size())) == str_match) {
			//�޳�"pause_filename["
			int pre_length = int(str_match.size() + 1);//"pause_filename["�ĳ���
			std::string str = files[i].substr(pre_length);
			//�޳�"].dat"
			int post_length = 5;//"].dat"�ĳ���
			for (int i = 0; i < post_length; i++) {
				str.pop_back();
			}
			//��xxx����num
			int num = std::stoi(str);
			if (num > maxNstep) {
				index_maxNstep = i;
				maxNstep = num;
			}
		}
	}
	//���ޣ��򷵻�-1
	if (index_maxNstep == -1)return -1;
	std::ifstream infile(m_path + files[index_maxNstep]);
	if (!infile) {
		return -1;
	}

	std::cout << "Continue from: " << files[index_maxNstep] << std::endl;

	int state = -1;//1-Node, 2-Element
	int maxnodeID = 1;
	int maxedgeID = 1;
	int maxelementID = 1;
	int iline_set = 0;//��ʱ���������ڼ�¼��ȡsetʱ�����˶�����
	std::vector<std::vector<int>> edges_of_all_sets;//��ʱ�������洢��set��edgeID
	const int bufferLength = 600;//!!!С�ĳ��Ȳ���
	char buffer[bufferLength];
	std::string tLine;
	std::vector<std::string> tWords;
	std::vector<std::vector<std::string>> tWordsMatrix;
	std::vector<int> edge_ID_long;
	int nNullRow = 0;
	bool loop = 1;
	while (loop) {
		//get words and set state
		infile.getline(buffer, bufferLength);
		if (infile.eof())loop = 0;
		tLine = buffer;
		tWords = StringProcessor::splitString(tLine);
		if (nNullRow >= 10)loop = 0;
		if (tWords.size() == 0) {
			nNullRow++;
			continue;
		}//���У�ǿ�ȿ�ʼ��һ��ѭ������ֹtWords[0]�ڴ����
		else {
			if (tWords[0] == "t")state = 0;
			if (tWords[0] == "nodes:")state = 1;
			if (tWords[0] == "edges:")state = 2;
			if (tWords[0] == "elements:")state = 3;
			if (tWords[0] == "boundaries:")state = 4;
		}

		//translate
		if (state == 0 && tWords[0].substr(0, 1) != "t") {//t_solve, istep_solve
			GlobalPara::time::t_previous = std::stod(tWords[0]);
			GlobalPara::time::istep_previous = (int)std::stod(tWords[1]);
		}
		else if (state == 1 && tWords[0].substr(0, 1) != "n") {//nodes: ID, x, y
			Node_2D node;
			node.ID = (int)std::stod(tWords[0]);
			node.x = std::stod(tWords[1]);
			node.y = std::stod(tWords[2]);
			pFVM2D->nodes.push_back(node);
			maxnodeID = (std::max)(maxnodeID, node.ID);
		}
		else if (state == 3 && tWords[0].substr(0, 1) != "e") {//elements: ID, nodes, U, Ux, Uy
			Element_2D ele;
			ele.ID = (int)std::stod(tWords[0]);
			ele.nodes[0] = (int)std::stod(tWords[1]);
			ele.nodes[1] = (int)std::stod(tWords[2]);
			ele.nodes[2] = (int)std::stod(tWords[3]);
			ele.U[0] = std::stod(tWords[4]);
			ele.U[1] = std::stod(tWords[5]);
			ele.U[2] = std::stod(tWords[6]);
			ele.U[3] = std::stod(tWords[7]);
			ele.Ux[0] = std::stod(tWords[8]);
			ele.Ux[1] = std::stod(tWords[9]);
			ele.Ux[2] = std::stod(tWords[10]);
			ele.Ux[3] = std::stod(tWords[11]);
			ele.Uy[0] = std::stod(tWords[12]);
			ele.Uy[1] = std::stod(tWords[13]);
			ele.Uy[2] = std::stod(tWords[14]);
			ele.Uy[3] = std::stod(tWords[15]);
			pFVM2D->elements.push_back(ele);
			maxelementID = (std::max)(maxelementID, ele.ID);
		}
		else if (state == 4 && tWords[0].substr(0, 1) != "b") {//boundaries: ID, name, edgeIDs
			iline_set++;
			if (!isdigit(tWords[1][0])) {
				/*
				"1 inf"
				1.vBoundarySets������Ŀ������Ŀ�Ĳ��ֳ�ʼ��
				2.edges_of_all_sets������Ŀ
				*/
				//set��ID,name��ʼ��
				VirtualBoundarySet_2D vb;
				vb.ID = (int)std::stod(tWords[0]);
				vb.name = tWords[1];
				pFVM2D->boundaryManager.boundaries.push_back(vb);
				//edges_of_all_sets������Ŀ
				edges_of_all_sets.push_back(std::vector<int>());
			}
			else {
				// 1.edges_of_all_sets�ĳ�ʼ��
				// �������һ��BoundarySet
				std::vector<int>& edges_of_current_set = edges_of_all_sets[edges_of_all_sets.size() - 1];
				std::vector<int> intVector = StringProcessor::stringVector_2_intVector(tWords);
				jointVectors(edges_of_current_set, intVector);
			}
		}
	}

	infile.close();


	process(maxnodeID, maxelementID, edges_of_all_sets);

	return 0;
}

std::string tryFindPauseFilePath() {
	/*
	
	*/
	FVM_2D* pFVM2D = FVM_2D::getInstance();
	std::string output_dir = FilePathManager::getInstance()->getOutputDirectory();
	std::vector<std::string> files = FilePathManager::getInstance()->getFiles(output_dir);//filesΪ�ַ����������洢��·���������ļ���
	int index_maxNstep = -1;
	//��Ѱ�Ƿ���pause_filename[xxx].dat�ļ������У�����xxx����
	int maxNstep = 0;
	int tmp_index = 0;//����ָʾ"[", "]"
	std::string tmp_str;
	for (int i = 0; i < files.size(); i++) {
		//��Ѱ��pause_filename��ͷ���ļ���
		std::string str_match = "pause_" + GlobalPara::basic::filename;
		if (files[i].substr(0, int(str_match.size())) == str_match) {
			//�޳�"pause_filename["
			int pre_length = int(str_match.size() + 1);//"pause_filename["�ĳ���
			std::string str = files[i].substr(pre_length);
			//�޳�"].dat"
			int post_length = 5;//"].dat"�ĳ���
			for (int i = 0; i < post_length; i++) {
				str.pop_back();
			}
			//��xxx����num
			int num = std::stoi(str);
			if (num > maxNstep) {
				index_maxNstep = i;
				maxNstep = num;
			}
		}
	}
	//���ޣ��򷵻ؿ�
	if (index_maxNstep == -1) {
		return "";
	}
	else {
		return output_dir + files[index_maxNstep];
	}
}

int read_1_1(int& maxnodeID, int& maxelementID, std::vector<SU2MeshReader::SimpleBoundary>& tmp_boudaries) {
	FVM_2D* pFVM2D = FVM_2D::getInstance();
	std::string pause_file_path = tryFindPauseFilePath();
	if (pause_file_path == "") {
		return -1;
	}
	std::ifstream infile(pause_file_path);
	if (!infile) {
		return -1;
	}
	std::cout << "Continue from: " << pause_file_path << std::endl;

	int state = -1;//1-Node, 2-Element
	const int bufferLength = 600;//!!!С�ĳ��Ȳ���
	char buffer[bufferLength];
	std::string tLine;
	std::vector<std::string> tWords;
	std::vector<std::vector<std::string>> tWordsMatrix;
	std::vector<int> edge_ID_long;
	int nNullRow = 0;
	bool loop = 1;
	while (loop) {
		//get words and set state
		infile.getline(buffer, bufferLength);
		if (infile.eof())loop = 0;
		tLine = buffer;
		tWords = StringProcessor::splitString(tLine);
		if (nNullRow >= 10)loop = 0;
		if (tWords.size() == 0) {
			nNullRow++;
			continue;
		}//���У�ǿ�ȿ�ʼ��һ��ѭ������ֹtWords[0]�ڴ����
		else {
			if (tWords[0] == "t")state = 0;
			if (tWords[0] == "nodes:")state = 1;
			if (tWords[0] == "edges:")state = 2;
			if (tWords[0] == "elements:")state = 3;
			if (tWords[0] == "boundaries:")state = 4;
		}

		//translate
		if (state == 0 && tWords[0].substr(0, 1) != "t") {//t_solve, istep_solve, CFL
			GlobalPara::time::t_previous = std::stod(tWords[0]);
			GlobalPara::time::istep_previous = (int)std::stod(tWords[1]);
			if (tWords.size() == 3) {
				GlobalPara::time::CFL = std::stod(tWords[2]);// ��ȡ�ϴα����CFL
			}
		}
		else if (state == 1 && tWords[0].substr(0, 1) != "n") {//nodes: ID, x, y
			Node_2D node;
			node.ID = (int)std::stod(tWords[0]);
			node.x = std::stod(tWords[1]);
			node.y = std::stod(tWords[2]);
			pFVM2D->nodes.push_back(node);
			maxnodeID = (std::max)(maxnodeID, node.ID);
		}
		else if (state == 3 && tWords[0].substr(0, 1) != "e") {//elements: ID, nodes, U
			Element_2D ele;
			ele.ID = (int)std::stod(tWords[0]);
			ele.nodes[0] = (int)std::stod(tWords[1]);
			ele.nodes[1] = (int)std::stod(tWords[2]);
			ele.nodes[2] = (int)std::stod(tWords[3]);
			ele.U[0] = std::stod(tWords[4]);
			ele.U[1] = std::stod(tWords[5]);
			ele.U[2] = std::stod(tWords[6]);
			ele.U[3] = std::stod(tWords[7]);
			pFVM2D->elements.push_back(ele);
			maxelementID = (std::max)(maxelementID, ele.ID);
		}
		else if (state == 4 && tWords[0].substr(0, 1) != "b") {//boundaries: ID, name, edgeIDs
			if (tWords.size() >= 2 && !isdigit(tWords[1][0])) {
				// ��������2���ʣ���2�����ʲ������֣�˵����boundary name��Ҫ��ȡname
				tmp_boudaries.push_back(SU2MeshReader::SimpleBoundary(tWords[2]));
			}
			else {
				// ����˵����edge��ţ���˵�ǰboundary���һ����
				// SU2�Ǵ�0��ʼ��Ҫ+1��continue file�Ǵ�1��ʼ������+1
				SU2MeshReader::SimpleBoundary& tmp_currentBoundary = tmp_boudaries[tmp_boudaries.size() - 1];
				tmp_currentBoundary.edges.push_back(SU2MeshReader::SimpleEdge(
					std::stoi(tWords[1]) + 0, std::stoi(tWords[2]) + 0
				));
			}
		}
	}

	infile.close();
	return 0;
}

void process_1_1(int maxNodeID, int maxElementID, std::vector<SU2MeshReader::SimpleBoundary>& tmp_boudaries) {
	LogWriter::logAndPrintError("unimplemented.@process_1_1\n");
	exit(-1);

}

int ContinueFileReader::readContinueFile_1_1() {
	/*
	Ѱ�ҵ�ǰĿ¼���Ƿ�����Ϊ��pause_�����ļ������У����ȡ�����ޣ��򷵻�-1
	���뷵��һ��ֵ���������-1�����Դ�ͷ��ʼ
	*/
	int maxnodeID = 1;
	int maxelementID = 1;
	std::vector<SU2MeshReader::SimpleBoundary>tmp_boudaries;
	read_1_1(maxnodeID, maxelementID, tmp_boudaries);
	process_1_1(maxnodeID, maxelementID, tmp_boudaries);

	return 0;


}

void ContinueFileReader::process(int maxnodeID, int maxelementID, std::vector<std::vector<int>>& edges_of_all_sets) {

	FVM_2D* pFVM2D = FVM_2D::getInstance();
	LogWriter::logAndPrint("Initiate pNodeTable, pElementTable, pEdgeTable\n");
	pFVM2D->iniPNodeTable(maxnodeID);// ��Ҫ�Ƚ���nodes���飬֪��ÿ��node��ID��������nodes
	pFVM2D->iniPElementTable(maxelementID);// ��Ҫ�Ƚ���elements���飬֪��ÿ��element��ID��������elements
	pFVM2D->iniEdges();// ��װedges
	pFVM2D->iniPEdgeTable();// ����pEdges���顣��Ҫ�Ƚ���edges���飬֪��ÿ��edge��ID��������edges����Ҫ����iniEdges()����
	pFVM2D->iniNode_neighborElements();
	LogWriter::logAndPrint("Calculate Element xy, edge length\n");
	pFVM2D->iniElement_xy_pEdges();
	pFVM2D->iniEdges_lengths();

	LogWriter::logAndPrint("Initialize boundary condition\n");
	// ��ʼ��boundaryManager.boundaries��pEdges��������ʼ��pEdge��ָ��edge��setID����Ϣ���ò��ֹ����������溯��
	// ǰ����������Ҫ��ʼ��pEdgeTable
	if (pFVM2D->pEdgeTable.size() == 0) {
		LogWriter::logAndPrint("Error: uninitialized pEdgeTable. (BoundaryManager_2D::iniBoundarySetPEdges)\n");
		throw "uninitialized pEdgeTable";
	}
	else {
		//��ʼ��boundaryManager.boundaries��ÿ��set��pEdges
		for (int iset = 0; iset < edges_of_all_sets.size(); iset++) {
			std::vector<int>& edgeIDs = edges_of_all_sets[iset];//��iset��set��edgeIDs
			for (int iw = 0; iw < edgeIDs.size(); iw++) {//��iset��set�ĵ�iw��edge
				//��f->pEdgeTable�У�����edgeID��ѯpEdge����ɳ�ʼ��
				int edge_ID = edgeIDs[iw];
				Edge_2D* pE = pFVM2D->pEdgeTable[edge_ID];
				pFVM2D->boundaryManager.boundaries[iset].pEdges.push_back(pE);
			}
		}
	}

	// ����edges�е�setID������boundaryManager.boundaries��type
	// ǰ����������edges��������boundaries��������boundaries��name��pEdges��ID
	pFVM2D->boundaryManager.iniBoundaryEdgeSetID_and_iniBoundaryType(pFVM2D);

}
