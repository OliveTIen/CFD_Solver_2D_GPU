#ifndef ELEMENT_2D_H
#define ELEMENT_2D_H
#include "Node_2D.h"

class FVM_2D;
class Edge_2D;

class Element_2D {
	//�����ε�Ԫ����ͨ��ֵ���֣����ø�˹����
	//�ı��ε�Ԫ
public:
	//�ṹ��Ϣ
	int ID = -1;
	int node_num = 3;//��Ԫ����
	int nodes[4]{-1,-1,-1,-1};//node ID ֮������ID������ָ�룬����Ϊ��ȡ�ļ���ֱ�Ӷ���ElementID-NodeID-NodeID-NodeID
	Edge_2D* pEdges[4]{};//�Ѿ��ڶ�ȡ�ļ�ʱ��ʼ����������
	//������Ϣ
	myfloat x = 0;//��Ԫ�������ꡣʹ��ǰ���calxy!
	myfloat y = 0;//��Ԫ��������
	myfloat area = 0.0;
	myfloat U[4] = { 1,0,0,1.1e5 / 0.4 };//�غ�����,��u,��v,��E
	myfloat Ux[4] = { 0,0,0,0 };//
	myfloat Uy[4] = { 0,0,0,0 };//
	myfloat Flux[4]{};//4���غ�������ֵͨ����ÿ�μӼ�ǰҪ����
	//������Ϣ
	myfloat deltaeig;//Roe��ʽ����ͨ��ʱ�õ���ÿһ��������

	//static Eigen::MatrixXi si_ti;//����Ĳ�������
	//static Eigen::MatrixXd GaussPointMatrix;//����Ĳ�������
	int GPUID = -1;

public:
	//������� ��˷��������������
	myfloat calArea(FVM_2D*);
	//Ѱ�����ڵ�Ԫ������ָ�����飬�����СΪ3��Ԫ�ذ�edge˳�����У�����Ϊnullptr
	std::vector<Element_2D*> findNeighbor();
	//Ѱ�����ڵ�Ԫ������ָ�����飬�޳���ָ��
	std::vector<Element_2D*> findNeighbor_withoutNullptr();
	//���㵽����ھӵľ���
	myfloat calDistanceFromNearestNeighbor(FVM_2D* f);
	//[δʹ��]����Ux��Uy������������
	void updateSlope_Barth(FVM_2D* f);
	//�ݶ�������������Ux��Uy����ֹ����
	void restructor_in_updateSlope_Barth(FVM_2D* f);
	//(�ع�)���ݵ�Ԫ�ֲ�����������ĳ��U��ʹ�����Էֲ�ǰ���뱣֤UxUy�����µ�
	void get_U(myfloat xpoint, myfloat ypoint, myfloat* _U);
	Eigen::Vector4d get_U(myfloat xpoint, myfloat ypoint);
	//���������ķ��غ���U2�����Բ�ֵ
									//δ���
	//�غ���ת�ٶ�
	std::vector<myfloat> U2uv(const Eigen::Vector4d& U_);
	//[old]���ɱ߽�Ԫ
	void generateElementEdge(FVM_2D* f);
	//[old]�Ӻ���������ע��ĳ����
	void generateElementEdge_registerSingle(FVM_2D* f, int ID_0, int ID_1, int iEdge);
	myfloat calLambda(const myfloat gamma);
	//���㦫
	myfloat calLambdaFlux(FVM_2D* f);
};

#endif