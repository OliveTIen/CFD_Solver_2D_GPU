#ifndef BOUNDARYMANAGER_2D
#define BOUNDARYMANAGER_2D
#include "Node_2D.h"
#include "Element_2D.h"
#include "Edge_2D.h"

/*
�����Ŀ����ͨ��pEdge���ٲ���pEdge�����ı߽����͡�
��ѯ��ʽ���ȸ���pEdge��ѯsetID���ٸ���setID��ѯset������
ÿ��pEdge����һ��setID���������߼��ϵ�ID
���бߵļ�����

��ȷ�������
Node �ڵ㣬����������Ϣ
Element ��Ԫ�������ڵ�ID
Edge ��Ԫ�ıߡ���ά���Ϊ��
Boundary �߽磬����һϵ��pEdge���൱�ڱ߼���EdgeSet

*/

class FVM_2D;

//����ߵ�Ԫ���������ȡinp�߽�Ԫʱ��Ϊ����
class VirtualEdge_2D {
public:
	int nodes[2]{};
};

// ����߽缯�ϣ��Ƕ�һ���ߵĳ����ʾ
class VirtualBoundarySet_2D {
public:
	int ID = -1;     // ��1��ʼ��push_back֮ǰ�Զ���ID����������ָ��+1��
	std::string name;// ��Ҫ��������ΪwriteContinueFile����Ҫ���
	int type = -1;   // �߽����� �μ�head.h _BC_xxx����name�������
	
	int startID;     // ��ʱ����������ʼ��ʱ������
	int endID;       // ��ʱ����������ʼ��ʱ������

	std::vector<Edge_2D*> pEdges;

public:
	//��ѯpEdge�Ƿ���pEdges�У����ǣ��򷵻�ָ��(��0��ʼ)�����򷵻�-1
	int getEdgeIndex(Edge_2D* pEdge);
};

// �߽������
class BoundaryManager_2D {
public:
	struct PeriodPair {
	public:
		int bType = -1;
		int setID_0 = -1;//��1��ʼ
		int setID_1 = -1;
	public:
		
	};

	// �߽缯������
	std::vector<VirtualBoundarySet_2D> boundaries;
	// ���ڱ߽�����
	std::vector<PeriodPair> periodPairs;//����������ڱ߽�������


public:
	//��ints���������Էָ���������顣�����1�ף�1β��2�ף�2β��...
	static std::vector<int> compressSeveralSequences(const std::vector<int>& ints);
	// �߽���������תint
	static int getBoundaryTypeByName(std::string boundaryName);

	VirtualBoundarySet_2D* findSetByID(int ID);
	//����BoundarySet��name�õ�type�����߽�edge����setID��ǩ
	void iniBoundaryEdgeSetID_and_iniBoundaryType(FVM_2D* f);
	//��ʼ��Զ���߽��ruvp
	void ini_infBoundary_ruvp();
	//[debug]Զ���߽絥Ԫ�غ���ǿ�Ƹ�ֵ 0-Ma,AOA 101-[debug]�����ض�ֵ
	void setBoundaryElementU(int tag = 0);

	//����setID�õ���Ӧ�ı߽�ָ��
	VirtualBoundarySet_2D* getBoundarySetByID(const int setID);

	//���ڱ߽�
	//������ڱ߽�������
	void checkPeriodPairs();
	//�ҵ���֮��Ե�set�����������ڱ߽�
	VirtualBoundarySet_2D* getPairByID_periodicBoundary(const int setID);
	//�ҵ�ĳedge������pElement_R�����������ڱ߽�
	Element_2D* get_pElement_R_periodic(Edge_2D* pEdge);
	//�ҵ�ĳedge��Ӧ��edge�����������ڱ߽�
	Edge_2D* get_pairEdge_periodic(Edge_2D* pEdge);

	//�ҵ�������С����
};

#endif