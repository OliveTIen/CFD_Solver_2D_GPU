#ifndef BOUNDARYMANAGER_2D
#define BOUNDARYMANAGER_2D
#include "Node_2D.h"
#include "Element_2D.h"
#include "Edge_2D.h"

class FVM_2D;

//����ߵ�Ԫ���������ȡinp�߽�Ԫʱ��Ϊ����
class VirtualEdge_2D {
public:
	int nodes[2]{};
};

class VirtualBoundarySet_2D {
public:
	int ID = -1;//��1��ʼ��push_back֮ǰ�Զ���ID����������ָ��+1��
	std::string name;
	int type = -1;//�߽����� �μ�head.h _BC_xxx����name�������
	
	int startID;//��ʱ����������ʼ��ʱ������
	int endID;//��ʱ����������ʼ��ʱ������

	std::vector<Edge_2D*> pEdges;

public:
	//��ѯpEdge�Ƿ���pEdges�У����ǣ��򷵻�ָ��(��0��ʼ)�����򷵻�-1
	int get_pEdge_index(Edge_2D* pEdge);
};

class BoundaryManager_2D {
public:
	//
	std::vector<VirtualBoundarySet_2D> vBoundarySets;

	//���ж���ṹ�壬���ڷ�װ��һ�����ڱ߽�
	struct PeriodPair {
		int bType = -1;
		int setID_0 = -1;//��1��ʼ
		int setID_1 = -1;
	};
	std::vector<PeriodPair> periodPairs;//����������ڱ߽�������


public:
	//��ʼ��vBoundarySet��pEdges������f->edges��ʼ����ʹ��(������readMeshFile��ʹ��)
	void iniBoundarySetPEdges_in_readMeshFile(FVM_2D* f, std::vector<VirtualEdge_2D>& vBoundaryEdges);
	//��ʼ��vBoundarySet��pEdges������f->pEdgeTable��ʼ����ʹ��(������readContinueFile��ʹ��)
	void iniBoundarySetPEdges_in_readContinueFile(FVM_2D* f, std::vector<std::vector<int>>& set_edge_ID);
	//��ints���������Էָ���������顣�����1�ף�1β��2�ף�2β��...
	std::vector<int> splitInts(const std::vector<int>& ints);
	//����2������ӵ���1������ĩβ
	void attachToVector(std::vector<int>& v1, const std::vector<int>& v2);
	//
	VirtualBoundarySet_2D* findSetByID(int ID);
	//����BoundarySet��name�õ�type�����߽�edge����setID��ǩ
	int iniBoundaryEdgeSetID_and_iniBoundaryType(FVM_2D* f);
	//��ʼ��Զ���߽��ruvp
	void ini_infBoundary_ruvp();
	//[debug]Զ���߽絥Ԫ�غ���ǿ�Ƹ�ֵ 0-Ma,AOA 101-[debug]�����ض�ֵ
	void setBoundaryElementU(int tag = 0);

	//����setID�õ���Ӧ�ı߽�ָ��
	VirtualBoundarySet_2D* getBoundarySetByID(const int setID);

	//���ڱ߽�
	//�ҵ���֮��Ե�set�����������ڱ߽�
	VirtualBoundarySet_2D* getPairByID_periodicBoundary(const int setID);
	//�ҵ�ĳedge������pElement_R�����������ڱ߽�
	Element_T3* get_pElement_R_periodic(Edge_2D* pEdge);
	//�ҵ�ĳedge��Ӧ��edge�����������ڱ߽�
	Edge_2D* get_pairEdge_periodic(Edge_2D* pEdge);

	//�ҵ�������С����
};

#endif