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
	std::vector<Node_2D*> pNodeTable;// 将ID作为行号，存入指针
	std::vector<Element_2D*> pElementTable;
	std::vector<Edge_2D*> pEdgeTable;// pEdgeTable[edges[i].ID] = &(edges[i]);

	bool hasInitElementXY = false;// 已经初始化单元中心坐标
	bool hasInitEdgeLengths = false;// 已经初始化edges的length, refLength

public:
	//static FVM_2D* getInstance();
	static FVM_2D* getInstance();
	
	// 求解
	//void solve_CPU2(std::string suffix_out, std::string suffix_info);
	// 检查非法值
	bool isNan();
	// 是否达到稳态
	bool isStable(std::vector<Element_2D> old);


	// 初始化pNodeTable。根据nodes将ID作为行号，存入指针
	void iniPNodeTable(int maxnodeID);
	// 初始化edges的ID, nodes, pElement_L, pElement_R，并计算单元面积，存进单元area。须在读取单元和节点后立刻调用
	void iniEdges();
	// 被上述函数调用。操作edges，添加或完善edge信息。添加：ID,nodes,Element_L，完善：Element_R
	void iniEdges_registerSingle(int n0, int n1, Element_2D* pE);
	// 初始化pEdgeTable。根据edges将ID作为行号，存入指针
	void iniPEdgeTable();
	// 初始化pElementTable
	void iniPElementTable(int maxelementID);
	// 初始化elements的坐标xy、pEdges
	void iniElement_xy_pEdges();
	// 初始化nodes的邻居单元
	void iniNode_neighborElements();
	// 初始化edges的length、refLength，存入edge，防止重复计算
	void iniEdges_lengths();

	// 搜寻edges，根据节点ID判断edge是否存在，不存在则return nullptr。不区分n0, n1先后
	Edge_2D* getEdgeByNodeIDs(int n0, int n1);
	// 根据节点ID返回Node指针
	Node_2D* getNodeByID(int ID);

private:
	void iniElement_xy_pEdges_parallel();

};

#endif