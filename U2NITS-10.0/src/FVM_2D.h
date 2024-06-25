#ifndef FVM2D_H
#define FVM2D_H
#include "Node_2D.h"
#include "head.h"
#include "Math.h"
#include "Element_2D.h"
#include "Edge_2D.h"
#include "BoundaryManager_2D.h"
#include "Solver_2D.h"
class FVM_2D {
private:
	FVM_2D() {};

public:
	static FVM_2D* pFVM2D;
	BoundaryManager_2D boundaryManager;
	std::vector<Node_2D> nodes;
	std::vector<Element_2D> elements;
	std::vector<Edge_2D> edges;
	std::vector<Node_2D*> pNodeTable;// ��ID��Ϊ�кţ�����ָ��
	std::vector<Element_2D*> pElementTable;
	std::vector<Edge_2D*> pEdgeTable;// pEdgeTable[edges[i].ID] = &(edges[i]);

	bool hasInitElementXY = false;// �Ѿ���ʼ����Ԫ��������
	bool hasInitEdgeLengths = false;// �Ѿ���ʼ��edges��length, refLength

public:
	//static FVM_2D* getInstance();
	static FVM_2D* getInstance();
	
	// ���
	//void solve_CPU2(std::string suffix_out, std::string suffix_info);
	// ���Ƿ�ֵ
	bool isNan();
	// �Ƿ�ﵽ��̬
	bool isStable(std::vector<Element_2D> old);


	// ��ʼ��pNodeTable������nodes��ID��Ϊ�кţ�����ָ��
	void iniPNodeTable(int maxnodeID);
	// ��ʼ��edges��ID, nodes, pElement_L, pElement_R�������㵥Ԫ����������Ԫarea�����ڶ�ȡ��Ԫ�ͽڵ�����̵���
	void iniEdges();
	// �������������á�����edges����ӻ�����edge��Ϣ����ӣ�ID,nodes,Element_L�����ƣ�Element_R
	void iniEdges_registerSingle(int n0, int n1, Element_2D* pE);
	// ��ʼ��pEdgeTable������edges��ID��Ϊ�кţ�����ָ��
	void iniPEdgeTable();
	// ��ʼ��pElementTable
	void iniPElementTable(int maxelementID);
	// ��ʼ��elements������xy��pEdges
	void iniElement_xy_pEdges();
	// ��ʼ��nodes���ھӵ�Ԫ
	void iniNode_neighborElements();
	// ��ʼ��edges��length��refLength������edge����ֹ�ظ�����
	void iniEdges_lengths();

	// ��Ѱedges�����ݽڵ�ID�ж�edge�Ƿ���ڣ���������return nullptr��������n0, n1�Ⱥ�
	Edge_2D* getEdgeByNodeIDs(int n0, int n1);
	// ���ݽڵ�ID����Nodeָ��
	Node_2D* getNodeByID(int ID);

private:
	void iniElement_xy_pEdges_parallel();

};

#endif