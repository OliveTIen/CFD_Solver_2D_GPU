#include "FVM_2D.h"
#include "../output/ConsolePrinter.h"
//#include "output/ResidualCalculator.h"
#include "../output/LogWriter.h"
#include "../global/StringProcessor.h"
#include "../global/SystemInfo.h"
#include "../global/FilePathManager.h"
#include "../output/HistWriter.h"
#include "../input/InpMeshReader.h"
#include "../input/SU2MeshReader.h"
#include "../solvers/GPUSolver2.h"
#include "../space/FluxGPU.h"
#include "../time/CalculateDt.h"
#include "../output/FieldWriter.h"

FVM_2D* FVM_2D::pFVM2D = nullptr;

FVM_2D* FVM_2D::getInstance() {
	if (pFVM2D == nullptr) {
		pFVM2D = new FVM_2D();
	}
	return pFVM2D;
}

void FVM_2D::iniPNodeTable(int maxnodeID) {
	if (nodes.size() == 0) {
		LogWriter::logError("null nodes exception, @iniPNodeTable\n");
		exit(-1);
	}
	pNodeTable.resize(maxnodeID + 1);
	for (int in = 0; in < nodes.size(); in++) {
		pNodeTable[nodes[in].ID] = &(nodes[in]);
	}
}

void FVM_2D::iniEdges() {
	// 确保单元、节点信息已读取。这是组装edge的前提条件
	if (elements.size() == 0) {
		LogWriter::logError("null elements exception, @FVM_2D::iniEdges\n");
		exit(-1);
	}
	if (nodes.size() == 0) {
		LogWriter::logError("null nodes exception, @FVM_2D::iniEdges\n");
		exit(-1);
	}
	if (pNodeTable.size() == 0) {
		LogWriter::logError("null pNodeTable exception, @FVM_2D::iniEdges\n");
		exit(-1);
	}
	// 叉乘计算单元面积，确保节点顺序是逆时针。
	for (int ie = 0; ie < elements.size(); ie++) {
		// 叉乘计算单元面积，若面积为负，应交换某两个节点顺序，使得节点顺序是逆时针
		Element_2D& element_i = elements[ie];
		double xn[3]{}, yn[3]{};
		for (int i = 0; i < 3; i++) {
			xn[i] = this->getNodeByID(element_i.nodes[i])->x;
			yn[i] = this->getNodeByID(element_i.nodes[i])->y;
		}
		double area = 0.5 * (xn[0] * (yn[1] - yn[2]) + xn[1] * (yn[2] - yn[0]) + xn[2] * (yn[0] - yn[1]));
		if (area < 0) {
			element_i.area = -area;
			int n1 = element_i.nodes[1];
			int n2 = element_i.nodes[2];
			element_i.nodes[1] = n2;
			element_i.nodes[2] = n1;
		}
		else {
			element_i.area = area;
		}
	}
	// 组装edge
	for (int ie = 0; ie < elements.size(); ie++) {

		int n0, n1, n2;
		n0 = elements[ie].nodes[0];
		n1 = elements[ie].nodes[1];
		n2 = elements[ie].nodes[2];
		iniEdges_registerSingle(n0, n1, &(elements[ie]));
		iniEdges_registerSingle(n1, n2, &(elements[ie]));
		iniEdges_registerSingle(n2, n0, &(elements[ie]));
	}
}

void FVM_2D::iniPEdgeTable() {
	if (edges.size() == 0) {
		LogWriter::logError("null edges exception, @iniPEdgeTable\n");
		exit(-1);
	}
	pEdgeTable.resize(edges.size() + 1);
	for (int i = 0; i < edges.size(); i++) {
		pEdgeTable[edges[i].ID] = &(edges[i]);
	}
}

void FVM_2D::iniNode_neighborElements() {
	// 前置条件：elements, elements.nodes, pNodeTable
	for (int ie = 0; ie < elements.size(); ie++) {
		for (int jn = 0; jn < 3; jn++) {
			Node_2D* pN = getNodeByID(elements[ie].nodes[jn]);
			pN->neighborElements.push_back(&(elements[ie]));
		}
	}
}

Edge_2D* FVM_2D::getEdgeByNodeIDs(int n0, int n1) {
	// 在edges中根据nodeIDs查询edge
	//n0, n1表示节点ID
	if(edges.size()==0)return nullptr;
	for (int i = 0; i < edges.size(); i++) {
		if (edges[i].nodes[0] == n0 && edges[i].nodes[1] == n1)return &(edges[i]);
		if (edges[i].nodes[0] == n1 && edges[i].nodes[1] == n0)return &(edges[i]);
	}
	return nullptr;
}

void FVM_2D::iniEdges_registerSingle(int n0, int n1, Element_2D* pE) {
	// 依赖于elements, edges,
	Edge_2D* pEdge = getEdgeByNodeIDs(n0, n1);
	if (pEdge == nullptr) {
		Edge_2D tmp_edge;
		tmp_edge.ID = (int)edges.size() + 1;
		tmp_edge.nodes[0] = n0;
		tmp_edge.nodes[1] = n1;
		tmp_edge.pElement_L = pE;
		edges.push_back(tmp_edge);
	}
	else {
		/*
		隐含假设：一个edge如果已经属于某个element A，则A是该edge的左element
		一个edge最多同时属于2个element
		*/ 
		pEdge->pElement_R = pE;
	}
}

void FVM_2D::iniPElementTable(int maxelementID) {
	if (elements.size() == 0) {
		LogWriter::logError("null elements exception, @iniPElementTable\n");
		exit(-1);
	}
	pElementTable.resize(maxelementID + 1);
	for (int ie = 0; ie < elements.size(); ie++) {
		pElementTable[elements[ie].ID] = &(elements[ie]);
	}
}

void FVM_2D::iniElement_xy_pEdges() {
	/*
	初始化单元的x,y,pEdges，便于以后查找
	前置条件：
	  需要有elements数组
	  需要已知element的nodeIDs
	  node的坐标已初始化
	  有edges数组，且edge的nodeIDs已初始化
	隐患：
	  三角形边界，pEdges只有3个。对于四边形，pEdges的node需要重新确定
	不足：
	  该方式时间复杂度为O(n^2)，首先要遍历单元，然后getEdgeByNodeIDs要遍历所有edge
	  究竟该如何并行？
	*/
	if (elements.size() == 0) {
		LogWriter::logError("null elements exception, @iniElement_xy_pEdges\n");
		exit(-1);
	}
	for (int ie = 0; ie < elements.size(); ie++) {
		//已经calxy，以后不必担心。但是必须要先读取node再读取element

		Element_2D* pElement = &(elements[ie]);
		//初始化单元中心坐标
		double sum_x = 0;
		double sum_y = 0;
		for (int i = 0; i < 3; i++) {
			sum_x += this->getNodeByID(pElement->nodes[i])->x;
			sum_y += this->getNodeByID(pElement->nodes[i])->y;
		}
		pElement->x = sum_x / 3.0;
		pElement->y = sum_y / 3.0;


		pElement->pEdges[0] = this->getEdgeByNodeIDs(pElement->nodes[0], pElement->nodes[1]);
		pElement->pEdges[1] = this->getEdgeByNodeIDs(pElement->nodes[1], pElement->nodes[2]);
		pElement->pEdges[2] = this->getEdgeByNodeIDs(pElement->nodes[2], pElement->nodes[0]);


	}
	hasInitElementXY = true;//已经初始化单元中心坐标。
}

void FVM_2D::iniEdges_lengths() {
	if (hasInitElementXY == false) {
		LogWriter::logAndPrintError("hasInitElementXY == false, in FVM_2D::iniEdges_lengths()\n");
	}
	for (int i = 0; i < edges.size(); i++) {
		edges[i].length = (float)edges[i].getLength();
		edges[i].refLength = (float)edges[i].getRefLength();
	}
	hasInitEdgeLengths = true;
}

Node_2D* FVM_2D::getNodeByID(int ID) {
	// 防止越界
	if (ID < 0 || ID >= pNodeTable.size()) {
		std::string error_msg = "ID out of range in FVM_2D::getNodeByID(), ID=" + std::to_string(ID) + "\n";
		LogWriter::logAndPrintError(error_msg);
		throw error_msg;
	}
	return pNodeTable[ID];
}

