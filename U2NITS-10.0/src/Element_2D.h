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
	//附加信息
	myfloat deltaeig;//Roe格式计算通量时用到。每一轮需清零

	//static Eigen::MatrixXi si_ti;//顶点的参数坐标
	//static Eigen::MatrixXd GaussPointMatrix;//顶点的参数坐标
	int GPUID = -1;

public:
	//计算面积 叉乘法计算三角形面积
	myfloat calArea(FVM_2D*);
	//寻找相邻单元，返回指针数组，数组大小为3，元素按edge顺序排列，可能为nullptr
	std::vector<Element_2D*> findNeighbor();
	//寻找相邻单元，返回指针数组，剔除空指针
	std::vector<Element_2D*> findNeighbor_withoutNullptr();
	//计算到最近邻居的距离
	myfloat calDistanceFromNearestNeighbor(FVM_2D* f);
	//[未使用]计算Ux、Uy，并进行限制
	void updateSlope_Barth(FVM_2D* f);
	//梯度限制器。修正Ux、Uy，防止过大
	void restructor_in_updateSlope_Barth(FVM_2D* f);
	//(重构)根据单元分布函数，计算某点U。使用线性分布前，须保证UxUy是最新的
	void get_U(myfloat xpoint, myfloat ypoint, myfloat* _U);
	Eigen::Vector4d get_U(myfloat xpoint, myfloat ypoint);
	//计算任意点的非守恒量U2。线性插值
									//未完成
	//守恒量转速度
	std::vector<myfloat> U2uv(const Eigen::Vector4d& U_);
	//[old]生成边界元
	void generateElementEdge(FVM_2D* f);
	//[old]子函数，用于注册某条边
	void generateElementEdge_registerSingle(FVM_2D* f, int ID_0, int ID_1, int iEdge);
	myfloat calLambda(const myfloat gamma);
	//计算Λ
	myfloat calLambdaFlux(FVM_2D* f);
};

#endif