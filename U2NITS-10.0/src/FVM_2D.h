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
public:
	static FVM_2D* pFVM2D;
	BoundaryManager_2D boundaryManager;
	Solver_2D solver;
	Math_2D math;
	std::vector<Node_2D> nodes;
	std::vector<Element_T3> elements;
	std::vector<Edge_2D> edges;
	std::vector<Node_2D*> pNodeTable;//将ID作为行号，存入指针
	std::vector<Element_T3*> pElement_T3Table;
	std::vector<Edge_2D*> pEdgeTable;
	std::vector<double> u_nodes;//节点u, size=nodes.size()。输出流场时，用作缓存。
	std::vector<double> v_nodes;
	std::vector<double> rho_nodes;
	std::vector<double> p_nodes;

	double dt = 1e10;
	double t_solve = 0;//续算时用到
	int istep_solve = 0;
	int iCurrentAutosaveFile = 0;//自动保存文件的指标，构成循环队列

public:
	FVM_2D();
	//主函数
	void run();
	//读取output文件夹中的续算文件(pause_*.dat)(绝对路径)，初始化网格变量，更新t_solve和istep_solve
	int readContinueFile();
	//初始化流场。续算不要使用该函数，否则会重新开始
	void setInitialCondition();
	//日志记录边界参数
	void logBoundaryCondition();
	//求解
	void solve_(std::string suffix_out, std::string suffix_info);
	//计算dt
	void caldt();
	//检查非法值
	bool isNan();
	//更新自动保存文件
	void updateAutoSaveFile(double t, int istep);
	//输出流场到文件
	void writeTecplotFile(std::string f_name, double t_current);
	//计算节点函数值，输出文件时用到
	void calculateNodeValue();
	//续算文件
	void writeContinueFile(std::string f_name, double t, int istep);
	//是否达到稳态
	bool isStable(std::vector<Element_T3> old);


	//初始化pNodeTable。根据nodes将ID作为行号，存入指针
	void iniPNodeTable(int maxnodeID);
	//初始化edges所有项。根据elements，注册edges，并对edges初始化
	void iniEdges();
	//初始化pEdgeTable。根据edges将ID作为行号，存入指针
	void iniPEdgeTable();
	//初始化节点剩余项。初始化邻居单元
	void iniNode_neighborElements();
	//搜寻edges，根据节点ID判断edge是否存在，不存在则return nullptr。不区分n0, n1先后
	Edge_2D* isEdgeExisted(int n0, int n1);
	//操作edges，添加或完善edge信息。添加：ID,nodes,Element_L，完善：Element_R
	void iniEdges_registerSingle(int n0, int n1, Element_T3* pE);
	//初始化pElementTable
	void iniPElementTable(int maxelementID);
	//初始化单元剩余项。计算element的xy、pEdges
	void iniElement_xy_pEdges();
	//根据节点ID返回Node指针
	Node_2D* getNodeByID(int ID);

	//等熵涡 Isentropic vortex
	void isentropicVortex(double x, double y, double xc, double yc, double chi, double& deltau, double& deltav, double& deltaT);
	void isentropicVortex_2(double xc, double yc, double chi, const double* ruvp0);
	//等熵涡误差文件头
	void writeFileHeader_isentropicVortex();
};

#endif