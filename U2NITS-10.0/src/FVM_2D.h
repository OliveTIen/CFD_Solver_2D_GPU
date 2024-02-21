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
	//Solver_2D solver;
	Math_2D math;
	std::vector<Node_2D> nodes;
	std::vector<Element_2D> elements;
	std::vector<Edge_2D> edges;
	std::vector<Node_2D*> pNodeTable;//��ID��Ϊ�кţ�����ָ��
	std::vector<Element_2D*> pElementTable;
	std::vector<Edge_2D*> pEdgeTable;
	std::vector<double> u_nodes;//�ڵ�u, size=nodes.size()���������ʱ���������档
	std::vector<double> v_nodes;
	std::vector<double> rho_nodes;
	std::vector<double> p_nodes;

	double dt = 1e10;
	double t_previous = 0;//����ʱ����ʼʱ��
	int istep_previous = 0;//����ʱ����ʼ��
	int iCurrentAutosaveFile = 0;//�Զ������ļ���ָ�꣬����ѭ������
	bool hasInitElementXY = false;//�Ѿ���ʼ����Ԫ��������
	bool hasInitEdgeLengths = false;//�Ѿ���ʼ��edges��length, refLength

public:
	void run_GPU();
	void solve_GPU();


	FVM_2D();
	//������
	void run();
	//��ȡoutput�ļ����е������ļ�(pause_*.dat)(����·��)����ʼ���������������t_solve��istep_solve
	int readContinueFile();
	//��ʼ�����������㲻Ҫʹ�øú�������������¿�ʼ
	void setInitialCondition();
	//���
	void solve_CPU(std::string suffix_out, std::string suffix_info);
	void solve_CPU2(std::string suffix_out, std::string suffix_info);
	//����dt
	void caldt();
	//���Ƿ�ֵ
	bool isNan();
	//�����Զ������ļ�
	void updateAutoSaveFile(double t, int istep);
	//����������ļ� �ɺ��������Ƽ�ʹ�ã�������solve_CPU()
	// solve_GPU�У���FieldWriterʵ�ָù���
	void writeTecplotFile(std::string f_name, double t_current);
	//����ڵ㺯��ֵ������ļ�ʱ�õ�
	void calculateNodeValue();
	//�����ļ�
	void writeContinueFile(std::string f_name, double t, int istep);
	//�Ƿ�ﵽ��̬
	bool isStable(std::vector<Element_2D> old);


	//��ʼ��pNodeTable������nodes��ID��Ϊ�кţ�����ָ��
	void iniPNodeTable(int maxnodeID);
	//��ʼ��edges��ID, nodes, pElement_L, pElement_R������elements��ע��edges������edges��ʼ��
	void iniEdges();
	//�������������á�����edges����ӻ�����edge��Ϣ����ӣ�ID,nodes,Element_L�����ƣ�Element_R
	void iniEdges_registerSingle(int n0, int n1, Element_2D* pE);
	//��ʼ��pEdgeTable������edges��ID��Ϊ�кţ�����ָ��
	void iniPEdgeTable();
	//��ʼ��pElementTable
	void iniPElementTable(int maxelementID);
	//��ʼ��elements������xy��pEdges
	void iniElement_xy_pEdges();
	//��ʼ��nodes���ھӵ�Ԫ
	void iniNode_neighborElements();
	//��ʼ��edges��length��refLength������edge����ֹ�ظ�����
	void iniEdges_lengths();

	//��Ѱedges�����ݽڵ�ID�ж�edge�Ƿ���ڣ���������return nullptr��������n0, n1�Ⱥ�
	Edge_2D* getEdgeByNodeIDs(int n0, int n1);
	//���ݽڵ�ID����Nodeָ��
	Node_2D* getNodeByID(int ID);

	//������ Isentropic vortex
	void isentropicVortex(double x, double y, double xc, double yc, double chi, double& deltau, double& deltav, double& deltaT);
	void isentropicVortex_2(double xc, double yc, double chi, const double* ruvp0);
	//����������ļ�ͷ
	void writeFileHeader_isentropicVortex();
};

#endif