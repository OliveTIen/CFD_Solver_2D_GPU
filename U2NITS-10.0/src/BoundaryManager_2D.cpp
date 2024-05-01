#include "FVM_2D.h"
#include "output/LogWriter.h"
#include "boundary_condition/BoundaryManager.h"


std::vector<int> BoundaryManager_2D::compressSeveralSequences(const std::vector<int>& ints) {
	/*
	�ú�������ѹ���洢���У���������ռ�ÿռ䵫����ʧ��Ϣ
	
	�������ints�Ƿֶ����������У�����Ǹ����и��ε����Ҷ˵�ֵ

	���磬"1-5,11-14,21-26"->"1,5,11,14,21,26"
	���룺1,2,3,4,5,11,12,13,14,21,22,23,24,25,26
	�����1,5,11,14,21,26

	ʵ��Ӧ���У����е�ÿ��Ԫ�ر�ʾһ���߽���Ԫ��ID��һ���������ж�Ӧһ���߽�

	*/

	std::vector<int> ret;
	if (ints.size() < 2) {//�գ�����ֻ��1��Ԫ��
		if (ints.size() == 1) {
			ret.push_back(ints[0]);
			ret.push_back(ints[0]);
		}
		return ret;//��-���ؿա���1��Ԫ��a-����a,a
	}
	ret.push_back(ints[0]);
	int last, current;
	last = ints[0];
	for (int i = 1; i < ints.size(); i++) {
		current = ints[i];
		if (current == last + 1) {
			if (i == (int)ints.size() - 1)ret.push_back(ints[i]);
		}
		else {//���
			ret.push_back(last);
			ret.push_back(current);
			if (i == (int)ints.size() - 1)
				ret.push_back(current);//��β����һ�������ظ�2��
		}
		last = current;
	}
	return ret;
}

VirtualBoundarySet_2D* BoundaryManager_2D::findSetByID(int ID) {
	if (ID < 1 || ID > boundaries.size()) {
		std::cout << "Error: illegal ID, in \"BoundaryManager_2D::findSetByID(int ID)\"\n";
		return nullptr;
	}
	return &(boundaries[ID - 1]);
}

void BoundaryManager_2D::iniBoundaryEdgeSetID_and_iniBoundaryType(FVM_2D* f) {
	//����Ŀ�ģ�����BoundarySet��name�õ�type�����߽�edge����setID��ǩ
	
	//���f->edges�Ƿ��ѳ�ʼ��
	if (f->edges.size() == 0) {
		LogWriter::logError("f->edges.size() == 0, @BoundaryManager_2D::iniBoundaryEdgeSetID_and_iniBoundaryType\n");
		exit(-1);
	}

	//���ǩ
	for (int is = 0; is < boundaries.size(); is++) {
		//�߽����͸���set��type���߽�ID����edge��setID
		int bType = U2NITS::BoundaryManager::boundaryNameToType(boundaries[is].name);
		boundaries[is].type = bType;
		for (int ie = 0; ie < boundaries[is].pEdges.size(); ie++) {
			boundaries[is].pEdges[ie]->setID = boundaries[is].ID;
		}

		//! ע��öδ�������ڶദ���޸�ʱ������Ҫ�޸Ķദ
		// ��ʼ��periodPairs��periodPairs�洢���ڱ߽�������Ϣ������������ڱ߽�������
		// periodPairs��boundaryManager�ĳ�Ա��ÿ��pair�洢int bType, int setID_0, int setID_1
		if (_BC_periodic_0 <= bType && bType <= _BC_periodic_9) {// ��������ʽ��𿪣���&&����
			//���periodPairs�Ƿ��Ѿ�¼���bType�����ǣ���ʹ�����еģ��������½�һ��PeriodPair tmp������periodPairs
			int index_pairs = -1;
			for (int i = 0; i < periodPairs.size(); i++) {
				if (periodPairs[i].bType == bType)index_pairs = i;
			}
			// ���ڣ���ֱ������
			if (index_pairs != -1) {
				periodPairs[index_pairs].setID_1 = boundaries[is].ID;
			}
			// �����ڣ����½�
			else {
				PeriodPair tmp;
				tmp.bType = bType;
				tmp.setID_0 = boundaries[is].ID;
				periodPairs.push_back(tmp);
			}
		}
	}

	//������ڱ߽����ȷ�Ժ�������
	checkPeriodPairs();
}

void BoundaryManager_2D::ini_infBoundary_ruvp() {
	//��ʼ��Զ���߽��ruvp
	//����Ϊdebug���趨�ض�ruvp
	//using namespace GlobalPara::boundaryCondition::_2D;
	//std::cout << "debug: BoundaryManager_2D::ini_infBoundary_ruvp\n";
	//inlet::ruvp[0] = 1.0;
	//inlet::ruvp[1] = 200;
	//inlet::ruvp[2] = 0;
	//inlet::ruvp[3] = 1.0e5;//101325
	//
	//outlet::ruvp[0] = 0.125;
	//outlet::ruvp[1] = 100;
	//outlet::ruvp[2] = 0;
	//outlet::ruvp[3] = 0.1e5;

	//��ʼ�ڲ�������FVM_2D::initField
}

void BoundaryManager_2D::setBoundaryElementU(int tag) {
	//[debug]Զ���߽絥Ԫ�غ���ǿ�Ƹ�ֵ
	using namespace GlobalPara::boundaryCondition::_2D;
	switch (tag) {
	case 0:
		for (int ie = 0; ie < FVM_2D::getInstance()->edges.size(); ie++) {
			Edge_2D& edge_tmp = FVM_2D::getInstance()->edges[ie];
			const int bType = boundaries[edge_tmp.setID - 1].type;
			if (bType == _BC_inf) {
				Math_2D::ruvp_2_U(inf::ruvp, edge_tmp.pElement_L->U, GlobalPara::constant::gamma);
			}
			else if (bType == _BC_inlet) {
				Math_2D::ruvp_2_U(inlet::ruvp, edge_tmp.pElement_L->U, GlobalPara::constant::gamma);
			}
			else if (bType == _BC_outlet) {
				Math_2D::ruvp_2_U(outlet::ruvp, edge_tmp.pElement_L->U, GlobalPara::constant::gamma);
			}
		}
		break;

	case 101:
		for (int ie = 0; ie < FVM_2D::getInstance()->edges.size(); ie++) {
			Edge_2D& edge_tmp = FVM_2D::getInstance()->edges[ie];
			const int bType = boundaries[edge_tmp.setID - 1].type;
			if (bType == _BC_wall_nonViscous) {
				//Math_2D::ruvp_2_U(inf::ruvp, edge_tmp.pElement_L->U, GlobalPara::constant::gamma);
				//edge_tmp.pElement_L->U[0] = ruvp[0];
				//edge_tmp.pElement_L->U[1] = ruvp[0] * ruvp[1];
				//edge_tmp.pElement_L->U[2] = ruvp[0] * ruvp[2];
			}
			else if (bType == _BC_inlet) {
				//Math_2D::ruvp_2_U(inlet::ruvp, edge_tmp.pElement_L->U, GlobalPara::constant::gamma);
				double* ruvp = inlet::ruvp;
				//edge_tmp.pElement_L->U[0] = ruvp[0];
				//edge_tmp.pElement_L->U[1] = ruvp[0] * ruvp[1];
				//edge_tmp.pElement_L->U[2] = ruvp[0] * ruvp[2];
				//edge_tmp.pElement_L->U[3] = Math_2D::get_rhoE(ruvp, GlobalPara::constant::gamma);
			}
			else if (bType == _BC_outlet) {
				//Math_2D::ruvp_2_U(outlet::ruvp, edge_tmp.pElement_L->U, GlobalPara::constant::gamma);
				double* ruvp = outlet::ruvp;
				//edge_tmp.pElement_L->U[0] = ruvp[0];
				//edge_tmp.pElement_L->U[1] = ruvp[0] * ruvp[1];
				//edge_tmp.pElement_L->U[2] = ruvp[0] * ruvp[2];

				//edge_tmp.pElement_L->U[3] = Math_2D::get_rhoE(ruvp, GlobalPara::constant::gamma);
			}
		}
		break;

	case 2:
		break;
	}

}

VirtualBoundarySet_2D* BoundaryManager_2D::getBoundarySetByID(const int setID) {
	if (setID <= 0 || setID > boundaries.size()) {//setID��ΧӦ��1��ʼ�����ΪvBoundarySets��size
		LogWriter::logAndPrintError("setID out of range. (getBoundarySetByID)\n");
		return nullptr;
	}
	return &(boundaries[setID - 1]);
}

void BoundaryManager_2D::checkPeriodPairs() {
	// ������ڱ߽����ȷ�Ժ�������
	// !�������������¼��裺
	// gmsh�����inp�ļ��У�edge����Ǵӵ�һ�������������ʱ�뻷��һȦ�����һ��periodPair�䷽���Ȼ���෴��
	bool b_loop = true;
	bool b_error_uninitialized = false;
	bool b_error_edge_num_not_match = false;
	bool b_error_not_translatable = false;
	for (int ip = 0; ip < periodPairs.size() && b_loop; ip++) {
		// ����Ӧ�߽��Ƿ���ڡ���������ڣ���������Ϊд���ڱ߽�ʱ��ÿ�����ڱ߽粻�ǳɶԳ��֣�����δ��ʼ��
		auto& pair = periodPairs[ip];
		if (pair.setID_0 == -1 || pair.setID_1 == -1) {
			b_error_uninitialized = true;
			b_loop = false;
			break;
		}

		const VirtualBoundarySet_2D& vbs0 = boundaries[pair.setID_0 - 1];
		const VirtualBoundarySet_2D& vbs1 = boundaries[pair.setID_1 - 1];

		// ����Ӧ�߽��Ƿ�������
		if (vbs0.pEdges.size() != vbs1.pEdges.size()) {
			b_error_edge_num_not_match = true;
			b_loop = false;
			break;
		}

		// ����Ӧ�߽��Ƿ�����ƽ�ƹ��� ���߽�1��ĳ��ƽ��(vx,vy)�õ�p0,�ж�p0��p1�Ƿ��غ�
		// ����ȡ��һ�Զ�Ӧ�㣬����λ��
		const int pEdge_size = int(vbs0.pEdges.size());
		double translate_x = vbs1.pEdges[pEdge_size - 1]->getx() - vbs0.pEdges[0]->getx();
		double translate_y = vbs1.pEdges[pEdge_size - 1]->gety() - vbs0.pEdges[0]->gety();
		for (int ie = 0; ie < pEdge_size && b_loop; ie++) {
			double x1_ref = vbs0.pEdges[ie]->getx() + translate_x;
			double y1_ref = vbs0.pEdges[ie]->gety() + translate_y;
			double x1 = vbs1.pEdges[pEdge_size - 1 - ie]->getx();
			double y1 = vbs1.pEdges[pEdge_size - 1 - ie]->gety();
			double residual_x = x1 - x1_ref;
			double residual_y = y1 - y1_ref;
			double residual_square = residual_x * residual_x + residual_y * residual_y;
			double translate_square = translate_x * translate_x + translate_y * translate_y;
			double residual_relative = residual_square / translate_square;// ��Բв�

			if (residual_relative >= 0.0001) {
				b_error_not_translatable = true;
				b_loop = false;
				break;
			}
		}
	}

	if (b_error_uninitialized) {
		std::stringstream ss;
		ss << "b_error_uninitialized, @BoundaryManager_2D::checkPeriodPairs().\n"
			<< "Maybe because perioic boundaries does not appear in pairs. Please write periodic boundary seperately.\n";
		LogWriter::logAndPrintError(ss.str());
		exit(-1);
	}
	if (b_error_edge_num_not_match) {
		std::stringstream ss;
		ss << "b_error_edge_num_not_match, @BoundaryManager_2D::checkPeriodPairs().\n"
			<< "Invalid periodic boundary, for nums of edges doesn't match.\n";
		LogWriter::logAndPrintError(ss.str());
		exit(-1);
	}
	if (b_error_not_translatable) {
		std::stringstream ss;
		ss << "b_error_not_translatable, @BoundaryManager_2D::checkPeriodPairs().\n"
			<< "Invalid periodic boundary, for not meeting translation conditions.\n";
		LogWriter::logAndPrintError(ss.str());
		exit(-1);
	}
}

VirtualBoundarySet_2D* BoundaryManager_2D::getPairByID_periodicBoundary(const int setID) {
	VirtualBoundarySet_2D* set_0 = getBoundarySetByID(setID);
	int setID_to_find = -1;//��Ѱ�ҵ�set
	const int bType = set_0->type;
	if (_BC_periodic_0 <= bType && bType <= _BC_periodic_9) {//�ж��Ƿ�Ϊ���ڱ߽�
		//�����ڱ߽磬��Ѱ�����
		for (int ip = 0; ip < periodPairs.size(); ip++) {
			if (periodPairs[ip].bType == bType) {//����ͬһ�����ڱ߽�
				if (periodPairs[ip].setID_0 == setID)setID_to_find = periodPairs[ip].setID_1;
				else setID_to_find = periodPairs[ip].setID_0;
				break;//����forѭ��
			}
		}
		//���ҵ���Ե�ID���򷵻�ָ�룬���򱨴�
		if (setID_to_find == -1) {
			LogWriter::logAndPrintError("pair not found. (getPairByID_periodicBoundary)\n");
			return nullptr;
		}
		else return getBoundarySetByID(setID_to_find);
	}
	else {
		LogWriter::logAndPrintError("not periodic type. (getPairByID_periodicBoundary)\n");
		return nullptr;
	}
}

Element_2D* BoundaryManager_2D::get_pElement_R_periodic(Edge_2D* pEdge_0) {
	//�������ܣ��ҵ�ĳedge������pElement_R�����������ڱ߽�

	//�ҵ�pEdge��Ӧ��boundarySet(pBoundary_0)���Լ����boundarySet(pBoundary_1)
	VirtualBoundarySet_2D* pBoundary_0 = FVM_2D::getInstance()->boundaryManager.getBoundarySetByID(pEdge_0->setID);//pEdge��Ӧ��boundarySet
	VirtualBoundarySet_2D* pBoundary_1 = FVM_2D::getInstance()->boundaryManager.getPairByID_periodicBoundary(pEdge_0->setID);
	//��ȡpEdge��pBoundary_0��pEdges�е����(��0��ʼ��ָ��)
	int index_0 = pBoundary_0->getEdgeIndex(pEdge_0);
	if (index_0 == -1) {
		LogWriter::logAndPrintError("edge not found. (BoundaryManager_2D::get_pElement_R_periodic)\n");
		return nullptr;
	}
	//����õ�����pBoundary_1�����
	const int eSize = int(pBoundary_0->pEdges.size());
	int index_1 = eSize - 1 - index_0;
	//�ҵ�pEdge_0��pBoundary_1�϶�Ӧ��pEdge_1
	Edge_2D* pEdge_1 = pBoundary_1->pEdges[index_1];
	//pEdge_1��pElement_L��ΪpEdge_0������pElement_R
	return pEdge_1->pElement_L;
}

Edge_2D* BoundaryManager_2D::get_pairEdge_periodic(Edge_2D* pEdge_0) {
	// !�ú���������������Ĭ��index_1 = eSize - 1 - index_0�����߽�ڵ�ID����ʱ�������
	// �������ܣ��ҵ�ĳedge��Ӧ��edge�����������ڱ߽�
	// 
	
	//�ҵ�pEdge��Ӧ��boundarySet(pBoundary_0)���Լ����boundarySet(pBoundary_1)
	VirtualBoundarySet_2D* pBoundary_0 = FVM_2D::getInstance()->boundaryManager.getBoundarySetByID(pEdge_0->setID);//pEdge��Ӧ��boundarySet
	VirtualBoundarySet_2D* pBoundary_1 = FVM_2D::getInstance()->boundaryManager.getPairByID_periodicBoundary(pEdge_0->setID);
	//��ȡpEdge��pBoundary_0��pEdges�е����(��0��ʼ��ָ��)
	int index_0 = pBoundary_0->getEdgeIndex(pEdge_0);
	if (index_0 == -1) {
		LogWriter::logAndPrintError("edge not found. (BoundaryManager_2D::get_pElement_R_periodic)\n");
		return nullptr;
	}
	//����õ�����pBoundary_1�����
	const int eSize = int(pBoundary_0->pEdges.size());
	int index_1 = eSize - 1 - index_0;
	//�ҵ�pEdge_0��pBoundary_1�϶�Ӧ��pEdge_1
	return pBoundary_1->pEdges[index_1];
}

int VirtualBoundarySet_2D::getEdgeIndex(Edge_2D* pEdge) {
	//��ѯpEdge�Ƿ���pEdges�У����ǣ��򷵻�ָ��(��0��ʼ)�����򷵻�-1
	for (int i = 0; i < pEdges.size(); i++) {
		if (pEdge == pEdges[i])return i;
	}
	return -1;
}
