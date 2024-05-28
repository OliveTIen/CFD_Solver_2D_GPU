#ifndef EDGE_2D_H
#define EDGE_2D_H
#include "head.h"

class FVM_2D;
class Element_2D;

//��Ԫ���ڵıߣ�һάʱΪ��
class Edge_2D {
public:
	int ID;
	int nodes[2];//�ڵ�ID
	int setID = -1;//�߽�ID(����ָ��+1����1��ʼ)��ͨ��setID���Եõ�setType
	Element_2D* pElement_L = nullptr;
	Element_2D* pElement_R = nullptr;
	myfloat length{};// �߳��ȡ���ʼ����ã�����������㡣Ϊ��С�ڴ�ռ�ã���float
	myfloat refLength{};// ���൥Ԫ���ľ��룬�����ݶȼ��㡣
	int GPUID = -1;

public:
	void getxy(FVM_2D* f, myfloat& x, myfloat& y);
	myfloat getx();//�������ᱨ��δ��������FVM_2D
	myfloat gety();
	//myfloat getLength(FVM_2D* f);
	myfloat getLength();
	// ����������Ԫ���ĵľ���
	myfloat getRefLength();
	//void U_2_F_lambda(const Eigen::Vector4d U, Eigen::Vector4d& F, myfloat& lambda);

	//���㷨�򣬴�Element_L��Element_R
	std::vector<myfloat> getDirectionN();
	void getDirectionN(myfloat& nx, myfloat& ny);
	//����߳��򣬴�node0��node1
	std::vector<myfloat> getDirectionT();
	//��ת����flag=-1��ʾ�����
	Eigen::Matrix4d calT(FVM_2D* f, myfloat flag=1);
	//��ȡ���е����꣬������/�ҵ�Ԫ�ֲ���������ȡU
	//Eigen::Vector4d get_UL();
	//Eigen::Vector4d get_UR();
	//void get_ULUR(Eigen::Vector4d& U_L, Eigen::Vector4d&U_R);
	void get_ULUR(myfloat* U_L, myfloat* U_R);
};

#endif