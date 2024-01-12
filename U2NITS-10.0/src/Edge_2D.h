#ifndef EDGE_2D_H
#define EDGE_2D_H
#include "head.h"

class FVM_2D;
class Element_T3;

//单元相邻的边，一维时为点
class Edge_2D {
public:
	int ID;
	int nodes[2];//节点ID
	int setID = -1;//边界ID(数组指标+1，从1开始)。通过setID可以得到setType
	Element_T3* pElement_L = nullptr;
	Element_T3* pElement_R = nullptr;

public:
	void getxy(FVM_2D* f, double& x, double& y);
	double getx();//用内联会报错：未定义类型FVM_2D
	double gety();
	double getLength(FVM_2D* f);
	double getLength();
	void U_2_F_lambda(const Eigen::Vector4d U, Eigen::Vector4d& F, double& lambda);

	//计算法向，从Element_L到Element_R
	std::vector<double> getDirectionN(FVM_2D* f);
	std::vector<double> getDirectionN();
	void getDirectionN(double& nx, double& ny);
	//计算边朝向，从node0到node1
	std::vector<double> getDirectionT(FVM_2D* f);
	//旋转矩阵。flag=-1表示逆矩阵
	Eigen::Matrix4d calT(FVM_2D* f, double flag=1);
	//获取边中点坐标，代入左/右单元分布函数，获取U
	Eigen::Vector4d get_UL();
	Eigen::Vector4d get_UR();
	void get_ULUR(Eigen::Vector4d& U_L, Eigen::Vector4d&U_R);
	void get_ULUR(double* U_L, double* U_R);
};

#endif