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
	// 计算dt。根据CFL数计算各单元局部dt，取最小值。如果t+dt超出设定T，则限制dt
	double caldt(double t, double T);
	// 检查非法值
	bool isNan();
	// 更新自动保存文件
	void updateAutoSaveFile(double t, int istep, int& iCurrentAutosaveFile);
	// 输出流场到文件 旧函数，不推荐使用，仅用于solve_CPU()
	// solve_GPU中，用FieldWriter实现该功能
	void writeTecplotFile(std::string f_name, double t_current,
		std::vector<double>& rho_nodes,
		std::vector<double>& u_nodes,
		std::vector<double>& v_nodes,
		std::vector<double>& p_nodes);
	// 计算节点函数值，输出文件时用到
	void calculateNodeValue(std::vector<double>& rho_nodes,
		std::vector<double>& u_nodes,
		std::vector<double>& v_nodes,
		std::vector<double>& p_nodes);
	// GPU版本
	//void calculateNodeValue_GPU(void* pSolver);
	// 续算文件
	void writeContinueFile(std::string f_name, double t, int istep);
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

	//// 等熵涡 Isentropic vortex
	//void isentropicVortex(double x, double y, double xc, double yc, double chi, double& deltau, double& deltav, double& deltaT);
	//void isentropicVortex_2(double xc, double yc, double chi, const double* ruvp0);
	//// 等熵涡误差文件头
	//void writeFileHeader_isentropicVortex();
private:
	void iniElement_xy_pEdges_parallel();

};

#endif