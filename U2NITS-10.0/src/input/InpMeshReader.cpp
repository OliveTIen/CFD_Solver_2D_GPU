#include "InpMeshReader.h"
#include <iostream>
#include <fstream>
#include "../global/FilePathManager.h"
#include "../global/GlobalPara.h"
#include "../output/LogWriter.h"
#include "../BoundaryManager_2D.h"
#include "../global/StringProcessor.h"
#include "../FVM_2D.h"

int InpMeshReader::readGmeshFile(std::string filepath) {
	/*
	读取网格文件，要求初始化以下内容：
	1. nodes。包括坐标、ID、neighboringElements
	2. elements。包括ID、单元类型、节点ID。此外，为减少后续计算，需要提前计算单元中心坐标
	3. edges。包括ID、节点ID、所属边集合ID、左右单元指针。此外，也需要提前计算length和refLength

	用boundarySet维护一个pEdgeTable而不是直接在edge中存boundaryType，是为了处理周期边界。
	在处理周期边界时，要保证边与边的一一对应性
	
	*/


	FVM_2D* pFVM2D = FVM_2D::getInstance();
	std::cout << "read mesh file\n";
	std::ifstream infile(filepath);
	if (!infile) {
		LogWriter::writeLogAndCout("Warning: fail to read mesh(*.inp). (InpMeshReader::readGmeshFile) \n");
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
	std::vector<VirtualEdge_2D> vBoundaryEdges;//临时变量。存储边界元节点ID。数据结构为nx2的矩阵
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
		if (state == 1 && tWords[0].substr(0, 1) != "*") {
			//"*NODE" 
			// 4个数，第1个是节点ID，后3个是坐标

			// 存储点 存进nodes
			Node_2D node;
			node.ID = (int)std::stod(tWords[0]);
			node.x = std::stod(tWords[1]);
			node.y = std::stod(tWords[2]);
			pFVM2D->nodes.push_back(node);
			maxnodeID = (std::max)(maxnodeID, node.ID);
		}
		else if (state == 2 && substate == 1 && tWords[0].substr(0, 1) != "*") {
			// *ELEMENT, type=T3D2, ELSET=Line1
			// 3个数，第1个是线元ID，后2个是节点编号
			// 名称(Line1)不重要，之后都会统一成inf之类的边界名称

			// 存储边界线元 存进vBoundaryEdges
			VirtualEdge_2D ve;//成员：nodes[2]
			ve.nodes[0] = std::stoi(tWords[1]);//1, 1, 6; EdgeID,node0,node1
			ve.nodes[1] = std::stoi(tWords[2]);
			vBoundaryEdges.push_back(ve);
		}
		else if (state == 2 && substate == 2 && tWords[0].substr(0, 1) != "*") {
			// *ELEMENT, type=CPS3, ELSET=Surface1
			// 4个数，第1个是三角元ID，后3个是节点编号

			// 存储三角元 存进elements
			Element_2D e;
			e.ID = (int)std::stod(tWords[0]);
			e.nodes[0] = (int)std::stod(tWords[1]);
			e.nodes[1] = (int)std::stod(tWords[2]);
			e.nodes[2] = (int)std::stod(tWords[3]);
			pFVM2D->elements.push_back(e);
			maxelementID = (std::max)(maxelementID, e.ID);
		}
		else if (state == 3) {
			// *ELSET,ELSET=inf
			// 每行最多10个数，都是单元ID
			if (tWords[0].substr(0, 1) == "*") {
				// 根据vInts，确定有多少段连续数列，即有多少边界
				// 下面函数用于压缩存储分段连续的数列，例如"1-5,11-14,21-26"->"1,5,11,14,21,26"，存进start_end
				std::vector<int> start_end = pFVM2D->boundaryManager.compressSeveralSequences(vInts);
				int nSets = (int)start_end.size() / 2;
				for (int is = 0; is < nSets; is++) {//暗含nSets!=0，说明setName已经被初始化
					VirtualBoundarySet_2D vb;// 一条具有名称的边界，由一串连续的Edge构成的集合
					vb.name = setName;// 边界名称
					vb.ID = (int)pFVM2D->boundaryManager.boundaries.size() + 1;
					vb.startID = start_end[is * 2 + 0];//0,2,4,... 0=is*2+0 // 首EdgeID
					vb.endID = start_end[is * 2 + 1];// 尾EdgeID
					pFVM2D->boundaryManager.boundaries.push_back(vb);
				}
				//重置
				vInts.clear();
				setName = tWords[1].substr(6);//"wall""outlet""inlet"...
			}
			else {
				for (int iword = 0; iword < tWords.size(); iword++) {
					/* 
					vInts用于存储"*ELSET,ELSET=inf"下面所有数字。一般是连续的，如果
					不连续，表示是2条边界
					*/
					vInts.push_back(std::stoi(tWords[iword]));
				}
			}
		}
	}
	infile.close();

	std::cout << "Generate and initiate elements..." << std::endl;
	/*
	已完成：
	建立nodes数组，初始化node的ID、x、y
	建立elements数组，初始化element的ID、nodeIDs
	*/

	pFVM2D->iniPNodeTable(maxnodeID);// 需要先建立nodes数组，知道每个node的ID。已完成
	pFVM2D->iniPElementTable(maxelementID);// 需要先建立elements数组，知道每个element的ID。已完成
	pFVM2D->iniEdges();// 建立edges数组。需要知道每个element的nodeIDs。已完成

	pFVM2D->iniPEdgeTable();// 建立pEdges数组。需要先建立edges数组，知道每个edge的ID。依赖于iniEdges
	
	pFVM2D->iniElement_xy_pEdges();// 需要
	pFVM2D->iniNode_neighborElements();
	pFVM2D->iniEdges_lengths();
	

	// 初始化boundaryManager.boundaries的pEdges
	// 前置条件：有boundaries向量，有edges向量且edges已初始化nodeIDs
	{
		/*
		函数目的：初始化vBoundarySets的每个set的pEdges
		vBoundaryEdges包含未经提取的边界edge信息
		boundaries[iset]的startID、endID指示了应当提取vBoundaryEdges的哪一段

		vBoundarySets的每个元素为VirtualBoundarySet_2D类型
		VirtualBoundarySet_2D类型：定义了一条边界，成员变量有：边界ID、name、startID、endID、pEdges，
		分别存储边界ID、边界名称、起始点和终止点的ID、各edge单元指针

		vBoundaryEdges的每个元素为VirtualEdge_2D类型
		VirtualEdge_2D类型：定义了一个edge单元，成员变量为nodes，存储了两个节点的ID
		*/
		std::vector<VirtualBoundarySet_2D>& boundaries = pFVM2D->boundaryManager.boundaries;
		int istart = 1; int iend = 10;
		for (int iset = 0; iset < boundaries.size(); iset++) {
			//指示应当提取vBoundaryEdges的哪一段
			istart = boundaries[iset].startID;
			iend = boundaries[iset].endID;
			//提取这一段，存入vBoundarySets[iset].pEdges
			for (int ie = istart; ie <= iend; ie++) {
				int n0 = vBoundaryEdges[ie].nodes[0];
				int n1 = vBoundaryEdges[ie].nodes[1];
				Edge_2D* pEdge = pFVM2D->getEdgeByNodeIDs(n0, n1);
				//pEdge->set = boundaries[iset].ID;
				boundaries[iset].pEdges.push_back(pEdge);
			}
		}
	}

	// 设置edges中的setID，设置boundaryManager.boundaries的type
	// 前置条件：有edges向量，有boundaries向量，且boundaries有name
	{
		int ret = pFVM2D->boundaryManager.iniBoundaryEdgeSetID_and_iniBoundaryType(pFVM2D);
		if (ret != 0)return ret;
	}

	return 0;

}

int InpMeshReader::readSU2File(std::string suffix) {
	std::cout << "read mesh file" << "\n";
	int error_code = SU2::readFile(suffix);
	// 错误处理
	switch (error_code) {
	case SU2::E_invalid_first_word:
		std::cout << "invalid first word.\n";
		break;
	case SU2::E_open_file_error:
		LogWriter::writeLogAndCout("Warning: fail to read mesh(*.su2). (InpMeshReader::readSU2File) \n");
		break;
	case SU2::E_meet_number_in_the_start:
		std::cout << "meet number in the start. please check su2 file.\n";
	case SU2::E_invalid_word_num:
		std::cout << "E_invalid_word_num.\n";
	}
	return 0;
}

int InpMeshReader::SU2::readFile(std::string suffix) {
	/*
设置状态为START

S_START:
读取一行，
如果是末尾，则设置状态为END
如果首字符是字母，则读取第一个单词，设置状态为相应状态

NPOIN:
读取一行，
如果是double double int int，则存入点坐标
如果首字符是字母，则读取状态


*/
	FVM_2D* pFVM2D = FVM_2D::getInstance();
	std::string dir = FilePathManager::getInstance()->getInputDirectory();
	std::ifstream infile(dir + GlobalPara::basic::filename + suffix);
	if (!infile) {
		return E_open_file_error;
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

	SU2::State currentState = SU2::State::S_START;
	std::string message;// 上一个状态需要传递的消息
	bool loopContinue = true;
	while (loopContinue) {
		// 根据状态进行相应处理
		switch (currentState) {
			//// 开始
		case S_START:
		{
			// 试图读取一行，存入buffer
			infile.getline(buffer, bufferLength);
			// 若到文件末尾，设置状态为END
			if (infile.eof()) {
				currentState = S_END;
				break;
			}
			tLine = buffer;
			// 若为空行，无任何处理
			if (tLine.empty()) {
				break;
			}
			// 若以字母开头，则根据第一个单词设置当前状态
			if (isupper(tLine[0])) {
				tLine = StringProcessor::replaceCharInString(tLine, '=', " = ");// 等号前后加空格
				tWords = StringProcessor::splitString(tLine);// 分词
				std::string firstWord = tWords[0];
				if (firstWord == "NDIME") currentState = S_NDIME;
				else {
					infile.close();
					return E_invalid_first_word;
				}
				break;
			}
			else {
				// 一开始就遇见数字，不正常
				return E_meet_number_in_the_start;
			}
		}

		case S_NDIME:
		{
			// 处理当前行
			if (tWords.size() != 3) {
				return E_invalid_word_num;
			}
			if (std::stoi(tWords[2]) != 2) {
				std::cout << "Invalid dimension.\n";
				return E_invalid_word;
			}

			// 读取一行，设置状态
			infile.getline(buffer, bufferLength);
			if (infile.eof()) {
				currentState = S_END;
				break;
			}
			tLine = buffer;
			if (tLine.empty()) {
				break;
			}
			if (isupper(tLine[0])) {
				tLine = StringProcessor::replaceCharInString(tLine, '=', " = ");// 等号前后加空格
				tWords = StringProcessor::splitString(tLine);// 分词
				std::string firstWord = tWords[0];
				if (firstWord == "NPOIN")currentState = S_NPOIN;
				else {
					infile.close();
					return E_invalid_first_word;
				}
				break;
			}
		}
			break;

		case S_NPOIN:
		{
			// 处理当前行
			if (tWords.size() == 4) {
				Node_2D node;
				node.ID = (int)std::stod(tWords[3]);// 此处与inp不同
				node.x = std::stod(tWords[0]);
				node.y = std::stod(tWords[1]);
				pFVM2D->nodes.push_back(node);
				maxnodeID = (std::max)(maxnodeID, node.ID);

			}

			// 读取下一行，设置状态
			infile.getline(buffer, bufferLength);
			if (infile.eof()) {
				currentState = S_END;
				break;
			}
			tLine = buffer;
			if (tLine.empty()) {
				break;
			}
			if (isupper(tLine[0])) {
				tLine = StringProcessor::replaceCharInString(tLine, '=', " = ");// 等号前后加空格
				tWords = StringProcessor::splitString(tLine);// 分词
				std::string firstWord = tWords[0];
				if (firstWord == "NELEM")currentState = S_NELEM;
				else {
					infile.close();
					return E_invalid_first_word;
				}
				break;
			}
		}
			break;

		case S_NELEM:
		{

		}

		case S_END:
			loopContinue = false;
			break;
		}
		
	}

	infile.close();
	return E_no_error;
}
