#ifndef FVM2D_H
#define FVM2D_H
#include "Node_2D.h"
#include "head.h"
#include "Math.h"
#include "Element_2D.h"
#include "Edge_2D.h"
#include "BoundaryManager_2D.h"
#include "Solver_2D.h"
class FVM_2D {
public:
	static FVM_2D* pFVM2D;
	BoundaryManager_2D boundaryManager;
	Solver_2D solver;
	Math_2D math;
	std::vector<Node_2D> nodes;
	std::vector<Element_T3> elements;
	std::vector<Edge_2D> edges;
	std::vector<Node_2D*> pNodeTable;//��ID��Ϊ�кţ�����ָ��
	std::vector<Element_T3*> pElement_T3Table;
	std::vector<Edge_2D*> pEdgeTable;
	std::vector<double> u_nodes;//�ڵ�u, size=nodes.size()���������ʱ���������档
	std::vector<double> v_nodes;
	std::vector<double> rho_nodes;
	std::vector<double> p_nodes;

	double dt = 1e10;
	double t_solve = 0;//����ʱ�õ�
	int istep_solve = 0;
	int iCurrentAutosaveFile = 0;//�Զ������ļ���ָ�꣬����ѭ������

public:
	FVM_2D();
	//������
	void run();
	//��ȡoutput�ļ����е������ļ�(pause_*.dat)(����·��)����ʼ���������������t_solve��istep_solve
	int readContinueFile();
	//��ʼ�����������㲻Ҫʹ�øú�������������¿�ʼ
	void setInitialCondition();
	//��־��¼�߽����
	void logBoundaryCondition();
	//���
	void solve_(std::string suffix_out, std::string suffix_info);
	//����dt
	void caldt();
	//���Ƿ�ֵ
	bool isNan();
	//�����Զ������ļ�
	void updateAutoSaveFile(double t, int istep);
	//����������ļ�
	void writeTecplotFile(std::string f_name, double t_current);
	//����ڵ㺯��ֵ������ļ�ʱ�õ�
	void calculateNodeValue();
	//�����ļ�
	void writeContinueFile(std::string f_name, double t, int istep);
	//�Ƿ�ﵽ��̬
	bool isStable(std::vector<Element_T3> old);


	//��ʼ��pNodeTable������nodes��ID��Ϊ�кţ�����ָ��
	void iniPNodeTable(int maxnodeID);
	//��ʼ��edges���������elements��ע��edges������edges��ʼ��
	void iniEdges();
	//��ʼ��pEdgeTable������edges��ID��Ϊ�кţ�����ָ��
	void iniPEdgeTable();
	//��ʼ���ڵ�ʣ�����ʼ���ھӵ�Ԫ
	void iniNode_neighborElements();
	//��Ѱedges�����ݽڵ�ID�ж�edge�Ƿ���ڣ���������return nullptr��������n0, n1�Ⱥ�
	Edge_2D* isEdgeExisted(int n0, int n1);
	//����edges����ӻ�����edge��Ϣ����ӣ�ID,nodes,Element_L�����ƣ�Element_R
	void iniEdges_registerSingle(int n0, int n1, Element_T3* pE);
	//��ʼ��pElementTable
	void iniPElementTable(int maxelementID);
	//��ʼ����Ԫʣ�������element��xy��pEdges
	void iniElement_xy_pEdges();
	//���ݽڵ�ID����Nodeָ��
	Node_2D* getNodeByID(int ID);

	//������ Isentropic vortex
	void isentropicVortex(double x, double y, double xc, double yc, double chi, double& deltau, double& deltav, double& deltaT);
	void isentropicVortex_2(double xc, double yc, double chi, const double* ruvp0);
	//����������ļ�ͷ
	void writeFileHeader_isentropicVortex();
};

#endif