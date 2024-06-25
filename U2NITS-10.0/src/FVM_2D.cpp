#include "FVM_2D.h"
#include "output/ConsolePrinter.h"
//#include "output/ResidualCalculator.h"
#include "output/LogWriter.h"
#include "global/StringProcessor.h"
#include "global/SystemInfo.h"
#include "global/FilePathManager.h"
#include "output/HistWriter.h"
#include "input/InpMeshReader.h"
#include "input/SU2MeshReader.h"
#include "solvers/GPUSolver2.h"
#include "space/FluxGPU.h"
#include "time/CalculateDt.h"
#include "output/FieldWriter.h"

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
	// ȷ����Ԫ���ڵ���Ϣ�Ѷ�ȡ��������װedge��ǰ������
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
	// ��˼��㵥Ԫ�����ȷ���ڵ�˳������ʱ�롣
	for (int ie = 0; ie < elements.size(); ie++) {
		// ��˼��㵥Ԫ����������Ϊ����Ӧ����ĳ�����ڵ�˳��ʹ�ýڵ�˳������ʱ��
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
	// ��װedge
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
	// ǰ��������elements, elements.nodes, pNodeTable
	for (int ie = 0; ie < elements.size(); ie++) {
		for (int jn = 0; jn < 3; jn++) {
			Node_2D* pN = getNodeByID(elements[ie].nodes[jn]);
			pN->neighborElements.push_back(&(elements[ie]));
		}
	}
}

Edge_2D* FVM_2D::getEdgeByNodeIDs(int n0, int n1) {
	// ��edges�и���nodeIDs��ѯedge
	//n0, n1��ʾ�ڵ�ID
	if(edges.size()==0)return nullptr;
	for (int i = 0; i < edges.size(); i++) {
		if (edges[i].nodes[0] == n0 && edges[i].nodes[1] == n1)return &(edges[i]);
		if (edges[i].nodes[0] == n1 && edges[i].nodes[1] == n0)return &(edges[i]);
	}
	return nullptr;
}

void FVM_2D::iniEdges_registerSingle(int n0, int n1, Element_2D* pE) {
	// ������elements, edges,
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
		�������裺һ��edge����Ѿ�����ĳ��element A����A�Ǹ�edge����element
		һ��edge���ͬʱ����2��element
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
	��ʼ����Ԫ��x,y,pEdges�������Ժ����
	ǰ��������
	  ��Ҫ��elements����
	  ��Ҫ��֪element��nodeIDs
	  node�������ѳ�ʼ��
	  ��edges���飬��edge��nodeIDs�ѳ�ʼ��
	������
	  �����α߽磬pEdgesֻ��3���������ı��Σ�pEdges��node��Ҫ����ȷ��
	���㣺
	  �÷�ʽʱ�临�Ӷ�ΪO(n^2)������Ҫ������Ԫ��Ȼ��getEdgeByNodeIDsҪ��������edge
	  ��������β��У�
	*/
	if (elements.size() == 0) {
		LogWriter::logError("null elements exception, @iniElement_xy_pEdges\n");
		exit(-1);
	}
	for (int ie = 0; ie < elements.size(); ie++) {
		//�Ѿ�calxy���Ժ󲻱ص��ġ����Ǳ���Ҫ�ȶ�ȡnode�ٶ�ȡelement

		Element_2D* pElement = &(elements[ie]);
		//��ʼ����Ԫ��������
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
	hasInitElementXY = true;//�Ѿ���ʼ����Ԫ�������ꡣ
}

void FVM_2D::iniElement_xy_pEdges_parallel() {

	
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
	// ��ֹԽ��
	if (ID < 0 || ID >= pNodeTable.size()) {
		std::string error_msg = "ID out of range in FVM_2D::getNodeByID(), ID=" + std::to_string(ID) + "\n";
		LogWriter::logAndPrintError(error_msg);
		throw error_msg;
	}
	return pNodeTable[ID];
}

bool FVM_2D::isStable(std::vector<Element_2D> old) {
	double sum = 0;
	double res;
	for (int ie = 0; ie < elements.size(); ie++) {
		//sum += abs(elements[ie].calculateu() - old[ie].calculateu());
		res = elements[ie].U[1] - old[ie].U[1];
		sum += abs(res);
	}
	//std::cout << sum << std::endl;
	if (sum <= 1e-12)return 1;
	else return 0;
}

bool FVM_2D::isNan() {
	// ���ÿ��Ԫ�ص��غ�������������쳣ֵ
	bool is_nan = 0;
	for (int ie = 0; ie < elements.size(); ie++) {
		
		std::string str;
		if (isnan(elements[ie].U[0])) {// ������õ���ϵͳ��isnan����
			is_nan = 1;
			const double& x = elements[ie].x;
			const double& y = elements[ie].y;
			str = 
				"rho ==\"NaN\", in element (x=" + std::to_string(x) + ", y=" + std::to_string(y) 
				+ ", U[0,1,2,3]=" + std::to_string(elements[ie].U[0]) + ", " + std::to_string(elements[ie].U[1])
				+ ", " + std::to_string(elements[ie].U[2]) + ", " + std::to_string(elements[ie].U[3]) + "\n";
			//break;
		}
		else if(elements[ie].U[0] < 0) {
			const double& x = elements[ie].x;
			const double& y = elements[ie].y;
			str =
				"rho < 0, in element (x=" + std::to_string(x) + ", y=" + std::to_string(y)
				+ ", U[0,1,2,3]=" + std::to_string(elements[ie].U[0]) + ", " + std::to_string(elements[ie].U[1])
				+ ", " + std::to_string(elements[ie].U[2]) + ", " + std::to_string(elements[ie].U[3]) + "\n";
		}
		LogWriter::logAndPrint(str,LogWriter::Warning);
	}
	return is_nan;
}
