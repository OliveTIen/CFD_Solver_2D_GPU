#include "FVM_2D.h"
#include "output/LogWriter.h"

void BoundaryManager_2D::iniBoundarySetPEdges_in_readMeshFile(FVM_2D* f, std::vector<VirtualEdge_2D>& vBoundaryEdges) {
	//����Ŀ�ģ���ʼ��vBoundarySets��ÿ��set��pEdges
	//vBoundaryEdges����δ����ȡ�ı߽�edge��Ϣ
	//vBoundarySets[iset]��startID��endIDָʾ��Ӧ����ȡvBoundaryEdges����һ��
	
	//vBoundarySets��ÿ��Ԫ��ΪVirtualBoundarySet_2D����
	//VirtualBoundarySet_2D���ͣ�������һ���߽磬��Ա�����У��߽�ID��name��startID��endID��pEdges��
	//�ֱ�洢�߽�ID���߽����ơ���ʼ�����ֹ���ID����edge��Ԫָ��

	//vBoundaryEdges��ÿ��Ԫ��ΪVirtualEdge_2D����
	//VirtualEdge_2D���ͣ�������һ��edge��Ԫ����Ա����Ϊnodes���洢�������ڵ��ID
	
	int istart = 1; int iend = 10;
	for (int iset = 0; iset < vBoundarySets.size(); iset++) {
		//ָʾӦ����ȡvBoundaryEdges����һ��
		istart = vBoundarySets[iset].startID;
		iend = vBoundarySets[iset].endID;
		//��ȡ��һ�Σ�����vBoundarySets[iset].pEdges
		for (int ie = istart; ie <= iend; ie++) {
			int n0 = vBoundaryEdges[ie].nodes[0];
			int n1 = vBoundaryEdges[ie].nodes[1];
			Edge_2D* pEdge = f->isEdgeExisted(n0, n1);
			//pEdge->set = vBoundarySets[iset].ID;
			vBoundarySets[iset].pEdges.push_back(pEdge);
		}
	}

}

void BoundaryManager_2D::iniBoundarySetPEdges_in_readContinueFile(FVM_2D* f, std::vector<std::vector<int>>& set_edge_ID) {
	//����Ŀ�ģ���ʼ��vBoundarySets��ÿ��set��pEdges��������ʼ��pEdge��ָ��edge��setID����Ϣ���ò��ֹ����������溯��
	//set_edge_IDΪ��ά���飬�������߽��edgeID
	//��f->pEdgeTable�У�����edgeID��ѯpEdge��һһ��Ӧ����

	//���ǰ������
	if (f->pEdgeTable.size() == 0) {
		LogWriter::writeLogAndCout("Error: uninitialized pEdgeTable. (BoundaryManager_2D::iniBoundarySetPEdges)\n");
		return;
	}
	//��ʼ��vBoundarySets��ÿ��set��pEdges���洢��edge��ָ��
	for (int iset = 0; iset < set_edge_ID.size(); iset++) {
		std::vector<int>& edgeIDs = set_edge_ID[iset];//��iset��set��edgeIDs
		for (int iw = 0; iw < edgeIDs.size(); iw++) {//��iset��set�ĵ�iw��edge
			//��f->pEdgeTable�У�����edgeID��ѯpEdge����ɳ�ʼ��
			int edge_ID = edgeIDs[iw];
			Edge_2D* pE = f->pEdgeTable[edge_ID];
			vBoundarySets[iset].pEdges.push_back(pE);
		}
	}
}

std::vector<int> BoundaryManager_2D::splitInts(const std::vector<int>& ints) {
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

void BoundaryManager_2D::attachToVector(std::vector<int>& v1, const std::vector<int>& v2) {
	for (int i = 0; i < v2.size(); i++) {
		v1.push_back(v2[i]);
	}
}

VirtualBoundarySet_2D* BoundaryManager_2D::findSetByID(int ID) {
	if (ID < 1 || ID > vBoundarySets.size()) {
		std::cout << "Error: illegal ID, in \"BoundaryManager_2D::findSetByID(int ID)\"\n";
		return nullptr;
	}
	return &(vBoundarySets[ID - 1]);
}

int BoundaryManager_2D::iniBoundaryEdgeSetID_and_iniBoundaryType(FVM_2D* f) {
	//����Ŀ�ģ�����BoundarySet��name�õ�type�����߽�edge����setID��ǩ
	
	//���f->edges�Ƿ��ѳ�ʼ��
	if (f->edges.size() == 0) {
		std::string str = "Warning: f->edges.size() == 0, in \"BoundaryManager_2D::iniEdgeBoundaryCondition(FVM_2D* f)\"\n";
		LogWriter::writeLog(str, 1);
		std::cout<< str;
		return -1;
	}

	//���ǩ
	for (int is = 0; is < vBoundarySets.size(); is++) {
		//�߽�����ת�߽�����
		int bType = -1;//�߽�����
		const std::string bName = vBoundarySets[is].name;//�߽�����
		if (bName == "wall" || bName == "obstacle") {
			if (GlobalPara::physicsModel::equation == _EQ_euler)bType = _BC_wall_nonViscous;
			else bType = _BC_wall_adiabat;
		}
		else if (bName == "inlet")bType = _BC_inlet;
		else if (bName == "outlet")bType = _BC_outlet;
		else if (bName == "inf")bType = _BC_inf;
		else if (bName == "symmetry")bType = _BC_symmetry;
		else if (bName.substr(0, 8) == "periodic") {//periodic_x
			int x = std::stoi(bName.substr(bName.size() - 1, 1));//ȡ���һ��(��size-1��)�ַ�
			bType = _BC_periodic_0 + x;//_BC_periodic_x = _BC_periodic_0 + x
		}

		//�߽����͸���set��type���߽�ID����edge��setID
		vBoundarySets[is].type = bType;
		for (int ie = 0; ie < vBoundarySets[is].pEdges.size(); ie++) {
			vBoundarySets[is].pEdges[ie]->setID = vBoundarySets[is].ID;
		}

		//��ʼ��periodPairs��periodPairs�洢���ڱ߽�������Ϣ������������ڱ߽�������
		if (_BC_periodic_0 <= bType && bType <= _BC_periodic_9) {//!ע�ⲻ����д��_BC_periodic_0 <= bType <= _BC_periodic_9
			//���periodPairs�Ƿ��Ѿ�¼���bType�����ǣ���ʹ�����еģ��������½�һ��PeriodPair tmp������periodPairs
			int check_p = -1;
			for (int ip = 0; ip < periodPairs.size(); ip++) {
				if (periodPairs[ip].bType == bType)check_p = ip;
			}
			if (check_p != -1) {
				periodPairs[check_p].setID_1 = vBoundarySets[is].ID;
			}
			else {
				PeriodPair tmp;
				tmp.bType = bType;
				tmp.setID_0 = vBoundarySets[is].ID;
				periodPairs.push_back(tmp);
			}
		}
	}

	//������ڱ߽����ȷ�Ժ�������
	//gmsh�����inp�ļ��У�edge����Ǵӵ�һ�������������ʱ�뻷��һȦ�����һ��periodPair�䷽���Ȼ���෴��
	int check_2 = 1;
	for (int ip = 0; ip < periodPairs.size() && check_2; ip++) {
		//ĳһ�Ա߽�
		const VirtualBoundarySet_2D& vbs0 = vBoundarySets[periodPairs[ip].setID_0 - 1];
		const VirtualBoundarySet_2D& vbs1 = vBoundarySets[periodPairs[ip].setID_1 - 1];
		//����Ӧ��������Ƿ�����ƽ�ƹ���
		if (vbs0.pEdges.size() != vbs1.pEdges.size()) {//����Ӧ����������ͬ
			check_2 *= 0;
		}
		else {
			const int eSize = int(vbs0.pEdges.size());
			//����ƽ������
			double vx, vy;//ƽ������v
			vx = vbs1.pEdges[eSize - 1]->getx() - vbs0.pEdges[0]->getx();
			vy = vbs1.pEdges[eSize - 1]->gety() - vbs0.pEdges[0]->gety();
			for (int ie = 0; ie < eSize && check_2; ie++) {
				//���߽�1��ĳ��ƽ��(vx,vy)�õ�p0,�ж�p0��p1�Ƿ��غ�
				double p0_x = vbs0.pEdges[ie]->getx() + vx;
				double p0_y = vbs0.pEdges[ie]->gety() + vy;
				double p1_x = vbs1.pEdges[eSize - 1 - ie]->getx();
				double p1_y = vbs1.pEdges[eSize - 1 - ie]->gety();
				double distance2 = (p1_x - p0_x) * (p1_x - p0_x) + (p1_y - p0_y) * (p1_y - p0_y);
				double vlength2 = vx * vx + vy * vy;
				if (distance2 >= 0.0001 * vlength2)check_2 *= 0;//|p1-p0|����0.01*|v|����Ϊ���غ�
			}
		}
		//�������򱨴�
		if (!check_2) {
			std::string str = "Error: invalid periodic boundary, for not meeting translation conditions. (BoundaryManager_2D::iniBoundaryEdge_SetType)\n";
			LogWriter::writeLogAndCout(str);
			return -1;
		}
	}
	return 0;
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
		for (int ie = 0; ie < FVM_2D::pFVM2D->edges.size(); ie++) {
			Edge_2D& edge_tmp = FVM_2D::pFVM2D->edges[ie];
			const int bType = vBoundarySets[edge_tmp.setID - 1].type;
			if (bType == _BC_inf) {
				Math_2D::ruvp_2_U(inf::ruvp, edge_tmp.pElement_L->U, Constant::gamma);
			}
			else if (bType == _BC_inlet) {
				Math_2D::ruvp_2_U(inlet::ruvp, edge_tmp.pElement_L->U, Constant::gamma);
			}
			else if (bType == _BC_outlet) {
				Math_2D::ruvp_2_U(outlet::ruvp, edge_tmp.pElement_L->U, Constant::gamma);
			}
		}
		break;

	case 101:
		for (int ie = 0; ie < FVM_2D::pFVM2D->edges.size(); ie++) {
			Edge_2D& edge_tmp = FVM_2D::pFVM2D->edges[ie];
			const int bType = vBoundarySets[edge_tmp.setID - 1].type;
			if (bType == _BC_wall_nonViscous) {
				//Math_2D::ruvp_2_U(inf::ruvp, edge_tmp.pElement_L->U, Constant::gamma);
				//edge_tmp.pElement_L->U[0] = ruvp[0];
				//edge_tmp.pElement_L->U[1] = ruvp[0] * ruvp[1];
				//edge_tmp.pElement_L->U[2] = ruvp[0] * ruvp[2];
			}
			else if (bType == _BC_inlet) {
				//Math_2D::ruvp_2_U(inlet::ruvp, edge_tmp.pElement_L->U, Constant::gamma);
				double* ruvp = inlet::ruvp;
				//edge_tmp.pElement_L->U[0] = ruvp[0];
				//edge_tmp.pElement_L->U[1] = ruvp[0] * ruvp[1];
				//edge_tmp.pElement_L->U[2] = ruvp[0] * ruvp[2];
				//edge_tmp.pElement_L->U[3] = Math_2D::get_rhoE(ruvp, Constant::gamma);
			}
			else if (bType == _BC_outlet) {
				//Math_2D::ruvp_2_U(outlet::ruvp, edge_tmp.pElement_L->U, Constant::gamma);
				double* ruvp = outlet::ruvp;
				//edge_tmp.pElement_L->U[0] = ruvp[0];
				//edge_tmp.pElement_L->U[1] = ruvp[0] * ruvp[1];
				//edge_tmp.pElement_L->U[2] = ruvp[0] * ruvp[2];

				//edge_tmp.pElement_L->U[3] = Math_2D::get_rhoE(ruvp, Constant::gamma);
			}
		}
		break;

	case 2:
		break;
	}

}

VirtualBoundarySet_2D* BoundaryManager_2D::getBoundarySetByID(const int setID) {
	if (setID <= 0 || setID > vBoundarySets.size()) {//setID��ΧӦ��1��ʼ�����ΪvBoundarySets��size
		LogWriter::writeLogAndCout("Error: setID out of range. (getBoundarySetByID)\n");
		return nullptr;
	}
	return &(vBoundarySets[setID - 1]);
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
			LogWriter::writeLogAndCout("Error: pair not found. (getPairByID_periodicBoundary)\n");
			return nullptr;
		}
		else return getBoundarySetByID(setID_to_find);
	}
	else {
		LogWriter::writeLogAndCout("Error: not periodic type. (getPairByID_periodicBoundary)\n");
		return nullptr;
	}
}

Element_T3* BoundaryManager_2D::get_pElement_R_periodic(Edge_2D* pEdge_0) {
	//�������ܣ��ҵ�ĳedge������pElement_R�����������ڱ߽�

	//�ҵ�pEdge��Ӧ��boundarySet(pBoundary_0)���Լ����boundarySet(pBoundary_1)
	VirtualBoundarySet_2D* pBoundary_0 = FVM_2D::pFVM2D->boundaryManager.getBoundarySetByID(pEdge_0->setID);//pEdge��Ӧ��boundarySet
	VirtualBoundarySet_2D* pBoundary_1 = FVM_2D::pFVM2D->boundaryManager.getPairByID_periodicBoundary(pEdge_0->setID);
	//��ȡpEdge��pBoundary_0��pEdges�е����(��0��ʼ��ָ��)
	int index_0 = pBoundary_0->get_pEdge_index(pEdge_0);
	if (index_0 == -1) {
		LogWriter::writeLogAndCout("Error: edge not found. (BoundaryManager_2D::get_pElement_R_periodic)\n");
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
	//�������ܣ��ҵ�ĳedge��Ӧ��edge�����������ڱ߽�
		//�ҵ�pEdge��Ӧ��boundarySet(pBoundary_0)���Լ����boundarySet(pBoundary_1)
	VirtualBoundarySet_2D* pBoundary_0 = FVM_2D::pFVM2D->boundaryManager.getBoundarySetByID(pEdge_0->setID);//pEdge��Ӧ��boundarySet
	VirtualBoundarySet_2D* pBoundary_1 = FVM_2D::pFVM2D->boundaryManager.getPairByID_periodicBoundary(pEdge_0->setID);
	//��ȡpEdge��pBoundary_0��pEdges�е����(��0��ʼ��ָ��)
	int index_0 = pBoundary_0->get_pEdge_index(pEdge_0);
	if (index_0 == -1) {
		LogWriter::writeLogAndCout("Error: edge not found. (BoundaryManager_2D::get_pElement_R_periodic)\n");
		return nullptr;
	}
	//����õ�����pBoundary_1�����
	const int eSize = int(pBoundary_0->pEdges.size());
	int index_1 = eSize - 1 - index_0;
	//�ҵ�pEdge_0��pBoundary_1�϶�Ӧ��pEdge_1
	return pBoundary_1->pEdges[index_1];
}

int VirtualBoundarySet_2D::get_pEdge_index(Edge_2D* pEdge) {
	//��ѯpEdge�Ƿ���pEdges�У����ǣ��򷵻�ָ��(��0��ʼ)�����򷵻�-1
	for (int i = 0; i < pEdges.size(); i++) {
		if (pEdge == pEdges[i])return i;
	}
	return -1;
}
