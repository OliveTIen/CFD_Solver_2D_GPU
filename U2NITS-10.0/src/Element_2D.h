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
	double x = 0;//��Ԫ�������ꡣʹ��ǰ���calxy!
	double y = 0;//��Ԫ��������
	double U[4] = { 1,0,0,1.1e5 / 0.4 };//�غ�����,��u,��v,��E
	double Ux[4] = { 0,0,0,0 };//
	double Uy[4] = { 0,0,0,0 };//
	double Flux[4]{};//4���غ�������ֵͨ����ÿ�μӼ�ǰҪ����
	//������Ϣ
	double deltaeig;//Roe��ʽ����ͨ��ʱ�õ���ÿһ��������

	//static Eigen::MatrixXi si_ti;//����Ĳ�������
	//static Eigen::MatrixXd GaussPointMatrix;//����Ĳ�������
	int GPUindex = -1;

public:
	//�������
	double calArea(FVM_2D*);
	//Ѱ�����ڵ�Ԫ������ָ������
	std::vector<Element_2D*> findNeighbor();
	std::vector<Element_2D*> findNeighbor_withoutNullptr();
	//���㵽����ھӵľ���
	double calDistanceFromNearestNeighbor(FVM_2D* f);
	//[δʹ��]����Ux��Uy������������
	void updateSlope_Barth(FVM_2D* f);
	//�ݶ�������������Ux��Uy����ֹ����
	void restructor_in_updateSlope_Barth(FVM_2D* f);
	//(�ع�)���ݵ�Ԫ�ֲ�����������ĳ��U��ʹ�����Էֲ�ǰ���뱣֤UxUy�����µ�
	void get_U(double xpoint, double ypoint, double* _U);
	Eigen::Vector4d get_U(double xpoint, double ypoint);
	//���������ķ��غ���U2�����Բ�ֵ
									//δ���
	//�غ���ת�ٶ�
	std::vector<double> U2uv(const Eigen::Vector4d& U_);
	//[old]���ɱ߽�Ԫ
	void generateElementEdge(FVM_2D* f);
	//[old]�Ӻ���������ע��ĳ����
	void generateElementEdge_registerSingle(FVM_2D* f, int ID_0, int ID_1, int iEdge);
	double calLambda(const double gamma);
	//���㦫
	double calLambdaFlux(FVM_2D* f);
};

#endif