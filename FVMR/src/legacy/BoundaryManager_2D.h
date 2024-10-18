#ifndef BOUNDARYMANAGER_2D
#define BOUNDARYMANAGER_2D
#include "Node_2D.h"
#include "Element_2D.h"
#include "Edge_2D.h"

/*
该类的目的是通过pEdge快速查找pEdge所属的边界类型。
查询方式是先根据pEdge查询setID，再根据setID查询set的类型
每个pEdge都有一个setID，即所属边集合的ID
所有边的集合由

明确几个概念：
Node 节点，包含坐标信息
Element 单元，包含节点ID
Edge 单元的边。三维情况为面
Boundary 边界，包含一系列pEdge。相当于边集合EdgeSet

*/

class FVM_2D;

//虚拟边单元。在最初读取inp边界元时作为缓存
class VirtualEdge_2D {
public:
	int nodes[2]{};
};

// 虚拟边界集合，是对一条边的抽象表示
class VirtualBoundarySet_2D {
public:
	int ID = -1;     // 从1开始。push_back之前自动赋ID，等于数组指标+1。
	std::string name;// 需要保留，因为writeContinueFile中需要输出
	int type = -1;   // 边界类型 参见head.h _BC_xxx。由name翻译而来
	
	int startID;     // 临时变量，仅初始化时有意义
	int endID;       // 临时变量，仅初始化时有意义

	std::vector<Edge_2D*> pEdges;

public:
	//查询pEdge是否在pEdges中，若是，则返回指标(从0开始)，否则返回-1
	int getEdgeIndex(Edge_2D* pEdge);
};

// 边界管理器
class BoundaryManager_2D {
public:
	struct PeriodPair {
	public:
		int bType = -1;
		int setID_0 = -1;//从1开始
		int setID_1 = -1;
	public:
		
	};

	// 边界集合数组
	std::vector<VirtualBoundarySet_2D> boundaries;
	// 周期边界数组
	std::vector<PeriodPair> periodPairs;//用来检查周期边界完整性


public:
	//将ints按照连续性分割成若干数组。输出：1首，1尾，2首，2尾，...
	static std::vector<int> compressSeveralSequences(const std::vector<int>& ints);

	VirtualBoundarySet_2D* findSetByID(int ID);
	//根据BoundarySet的name得到type；给边界edge打上setID标签
	void iniBoundaryEdgeSetID_and_iniBoundaryType(FVM_2D* f);
	//初始化远场边界的ruvp
	void ini_infBoundary_ruvp();
	//[debug]远场边界单元守恒量强制赋值 0-Ma,AOA 101-[debug]设置特定值
	void setBoundaryElementU(int tag = 0);

	//根据setID得到对应的边界指针
	VirtualBoundarySet_2D* getBoundarySetByID(const int setID);

	//周期边界
	//检查周期边界完整性
	void checkPeriodPairs();
	//找到与之配对的set。仅限于周期边界
	VirtualBoundarySet_2D* getPairByID_periodicBoundary(const int setID);
	//找到某edge的虚拟pElement_R。仅限于周期边界
	Element_2D* get_pElement_R_periodic(Edge_2D* pEdge);
	//找到某edge对应的edge。仅限于周期边界
	Edge_2D* get_pairEdge_periodic(Edge_2D* pEdge);

	//找到最大和最小坐标
};

#endif