#ifndef ELEMENT_2D_H
#define ELEMENT_2D_H
#include "Node_2D.h"

class FVM_2D;
class Edge_2D;

class Data_2D {
public:
	double U[4] = { 1,0,0,1e5 / 0.4 };//守恒量ρ,ρu,ρv,ρE
	double Ux[4] = { 0,0,0,0 };// 梯度
	double Uy[4] = { 0,0,0,0 };// 梯度
	double Flux[4]{};//4个守恒量的数值通量
};


class Element_T3 {
	//三角形单元 三角形单元用普通数值积分，不用高斯积分
public:
	//结构信息
	int ID = -1;
	int set = -1;//未使用
	int nodes[3];//node ID 之所以用ID而不是指针，是因为读取文件是直接读的ElementID-NodeID-NodeID-NodeID
	Edge_2D* pEdges[3];//已经在读取文件时初始化，放心用
	//数据信息
	double x = 0;//单元中心坐标。使用前务必calxy!
	double y = 0;//单元中心坐标
	double U[4] = { 1,0,0,1.1e5 / 0.4 };//守恒量ρ,ρu,ρv,ρE
	double Ux[4] = { 0,0,0,0 };//
	double Uy[4] = { 0,0,0,0 };//
	double Flux[4]{};//4个守恒量的数值通量。每次加减前要清零
	//附加信息
	double deltaeig;//Roe格式计算通量时用到。每一轮需清零

	//static Eigen::MatrixXi si_ti;//顶点的参数坐标
	//static Eigen::MatrixXd GaussPointMatrix;//顶点的参数坐标

public:
	//计算面积
	double calArea(FVM_2D*);
	//初始化单元中心坐标。对于静网格，只需要初始化一次即可
	void inixy(FVM_2D*);
	//寻找相邻单元，返回指针数组
	std::vector<Element_T3*> findNeighbor();
	std::vector<Element_T3*> findNeighbor_withoutNullptr();
	//计算到最近邻居的距离
	double calDistanceFromNearestNeighbor(FVM_2D* f);
	//计算Ux、Uy，并进行限制
	void updateSlope_Barth(FVM_2D* f);
	//梯度限制器。修正Ux、Uy，防止过大
	void restructor_in_updateSlope_Barth(FVM_2D* f);
	//(重构)根据单元分布函数，计算某点U。使用线性分布前，须保证UxUy是最新的
	void get_U(double xpoint, double ypoint, double* _U);
	Eigen::Vector4d get_U(double xpoint, double ypoint);
	//计算任意点的非守恒量U2。线性插值
									//未完成
	//守恒量转速度
	std::vector<double> U2uv(const Eigen::Vector4d& U_);
	//[old]生成边界元
	void generateElementEdge(FVM_2D* f);
	//[old]子函数，用于注册某条边
	void generateElementEdge_registerSingle(FVM_2D* f, int ID_0, int ID_1, int iEdge);
	void iniPEdges(FVM_2D* f);
	double calLambda(const double gamma);
	//计算Λ
	double calLambdaFlux(FVM_2D* f);
};

#endif