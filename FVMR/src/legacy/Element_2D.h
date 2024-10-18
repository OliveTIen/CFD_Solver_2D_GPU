#ifndef ELEMENT_2D_H
#define ELEMENT_2D_H
#include "Node_2D.h"

class FVM_2D;
class Edge_2D;

class Element_2D {
	//三角形单元用普通数值积分，不用高斯积分
	//四边形单元
public:
	//结构信息
	int ID = -1;
	int node_num = 3;//单元类型
	int nodes[4]{-1,-1,-1,-1};//node ID 之所以用ID而不是指针，是因为读取文件是直接读的ElementID-NodeID-NodeID-NodeID
	Edge_2D* pEdges[4]{};//已经在读取文件时初始化，放心用
	//数据信息
	myfloat x = 0;//单元中心坐标。使用前务必calxy!
	myfloat y = 0;//单元中心坐标
	myfloat area = 0.0;
	myfloat U[4] = { 1,0,0,1.1e5 / 0.4 };//守恒量ρ,ρu,ρv,ρE
	myfloat Ux[4] = { 0,0,0,0 };//
	myfloat Uy[4] = { 0,0,0,0 };//
	myfloat Flux[4]{};//4个守恒量的数值通量。每次加减前要清零


	int GPUID = -1;

public:
	//计算面积 叉乘法计算三角形面积
	myfloat calArea(FVM_2D*);
	//寻找相邻单元，返回指针数组，数组大小为3，元素按edge顺序排列，可能为nullptr
	std::vector<Element_2D*> findNeighbor();
	//寻找相邻单元，返回指针数组，剔除空指针
	std::vector<Element_2D*> findNeighbor_withoutNullptr();
	
};

#endif