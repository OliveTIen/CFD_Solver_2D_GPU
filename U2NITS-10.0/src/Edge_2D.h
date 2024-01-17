#ifndef EDGE_2D_H
#define EDGE_2D_H
#include "head.h"

class FVM_2D;
class Element_T3;

//��Ԫ���ڵıߣ�һάʱΪ��
class Edge_2D {
public:
	int ID;
	int nodes[2];//�ڵ�ID
	int setID = -1;//�߽�ID(����ָ��+1����1��ʼ)��ͨ��setID���Եõ�setType
	Element_T3* pElement_L = nullptr;
	Element_T3* pElement_R = nullptr;
	float length{};// �߳��ȡ���ʼ����ã�����������㡣Ϊ��С�ڴ�ռ�ã���float
	float refLength{};// ���൥Ԫ���ľ��룬�����ݶȼ��㡣

public:
	void getxy(FVM_2D* f, double& x, double& y);
	double getx();//�������ᱨ��δ��������FVM_2D
	double gety();
	//double getLength(FVM_2D* f);
	double getLength();
	// ����������Ԫ���ĵľ���
	float getRefLength();
	void U_2_F_lambda(const Eigen::Vector4d U, Eigen::Vector4d& F, double& lambda);

	//���㷨�򣬴�Element_L��Element_R
	std::vector<double> getDirectionN();
	void getDirectionN(double& nx, double& ny);
	//����߳��򣬴�node0��node1
	std::vector<double> getDirectionT();
	//��ת����flag=-1��ʾ�����
	Eigen::Matrix4d calT(FVM_2D* f, double flag=1);
	//��ȡ���е����꣬������/�ҵ�Ԫ�ֲ���������ȡU
	Eigen::Vector4d get_UL();
	Eigen::Vector4d get_UR();
	void get_ULUR(Eigen::Vector4d& U_L, Eigen::Vector4d&U_R);
	void get_ULUR(double* U_L, double* U_R);
};

#endif