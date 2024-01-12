#ifndef BOUNDARYMANAGER_2D
#define BOUNDARYMANAGER_2D
#include "Node_2D.h"
#include "Element_2D.h"
#include "Edge_2D.h"

class FVM_2D;

//虚拟边单元。在最初读取inp边界元时作为缓存
class VirtualEdge_2D {
public:
	int nodes[2]{};
};

class VirtualBoundarySet_2D {
public:
	int ID = -1;//从1开始。push_back之前自动赋ID，等于数组指标+1。
	std::string name;
	int type = -1;//边界类型 参见head.h _BC_xxx。由name翻译而来
	
	int startID;//临时变量，仅初始化时有意义
	int endID;//临时变量，仅初始化时有意义

	std::vector<Edge_2D*> pEdges;

public:
	//查询pEdge是否在pEdges中，若是，则返回指标(从0开始)，否则返回-1
	int get_pEdge_index(Edge_2D* pEdge);
};

class BoundaryManager_2D {
public:
	//
	std::vector<VirtualBoundarySet_2D> vBoundarySets;

	//类中定义结构体，便于封装：一对周期边界
	struct PeriodPair {
		int bType = -1;
		int setID_0 = -1;//从1开始
		int setID_1 = -1;
	};
	std::vector<PeriodPair> periodPairs;//用来检查周期边界完整性


public:
	//初始化vBoundarySet的pEdges。须在f->edges初始化后使用(仅限于readMeshFile中使用)
	void iniBoundarySetPEdges_in_readMeshFile(FVM_2D* f, std::vector<VirtualEdge_2D>& vBoundaryEdges);
	//初始化vBoundarySet的pEdges。须在f->pEdgeTable初始化后使用(仅限于readContinueFile中使用)
	void iniBoundarySetPEdges_in_readContinueFile(FVM_2D* f, std::vector<std::vector<int>>& set_edge_ID);
	//将ints按照连续性分割成若干数组。输出：1首，1尾，2首，2尾，...
	std::vector<int> splitInts(const std::vector<int>& ints);
	//将第2个数组加到第1个数组末尾
	void attachToVector(std::vector<int>& v1, const std::vector<int>& v2);
	//
	VirtualBoundarySet_2D* findSetByID(int ID);
	//根据BoundarySet的name得到type；给边界edge打上setID标签
	int iniBoundaryEdgeSetID_and_iniBoundaryType(FVM_2D* f);
	//初始化远场边界的ruvp
	void ini_infBoundary_ruvp();
	//[debug]远场边界单元守恒量强制赋值 0-Ma,AOA 101-[debug]设置特定值
	void setBoundaryElementU(int tag = 0);

	//根据setID得到对应的边界指针
	VirtualBoundarySet_2D* getBoundarySetByID(const int setID);

	//周期边界
	//找到与之配对的set。仅限于周期边界
	VirtualBoundarySet_2D* getPairByID_periodicBoundary(const int setID);
	//找到某edge的虚拟pElement_R。仅限于周期边界
	Element_T3* get_pElement_R_periodic(Edge_2D* pEdge);
	//找到某edge对应的edge。仅限于周期边界
	Edge_2D* get_pairEdge_periodic(Edge_2D* pEdge);

	//找到最大和最小坐标
};

#endif