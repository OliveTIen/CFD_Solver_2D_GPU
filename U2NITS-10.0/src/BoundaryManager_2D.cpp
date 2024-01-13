#include "FVM_2D.h"
#include "output/LogWriter.h"

void BoundaryManager_2D::iniBoundarySetPEdges_in_readMeshFile(FVM_2D* f, std::vector<VirtualEdge_2D>& vBoundaryEdges) {
	//函数目的：初始化vBoundarySets的每个set的pEdges
	//vBoundaryEdges包含未经提取的边界edge信息
	//vBoundarySets[iset]的startID、endID指示了应当提取vBoundaryEdges的哪一段
	
	//vBoundarySets的每个元素为VirtualBoundarySet_2D类型
	//VirtualBoundarySet_2D类型：定义了一条边界，成员变量有：边界ID、name、startID、endID、pEdges，
	//分别存储边界ID、边界名称、起始点和终止点的ID、各edge单元指针

	//vBoundaryEdges的每个元素为VirtualEdge_2D类型
	//VirtualEdge_2D类型：定义了一个edge单元，成员变量为nodes，存储了两个节点的ID
	
	int istart = 1; int iend = 10;
	for (int iset = 0; iset < vBoundarySets.size(); iset++) {
		//指示应当提取vBoundaryEdges的哪一段
		istart = vBoundarySets[iset].startID;
		iend = vBoundarySets[iset].endID;
		//提取这一段，存入vBoundarySets[iset].pEdges
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
	//函数目的：初始化vBoundarySets的每个set的pEdges，但不初始化pEdge所指的edge的setID等信息，该部分工作留给后面函数
	//set_edge_ID为二维数组，包含各边界的edgeID
	//在f->pEdgeTable中，根据edgeID查询pEdge，一一对应即可

	//检查前提条件
	if (f->pEdgeTable.size() == 0) {
		LogWriter::writeLogAndCout("Error: uninitialized pEdgeTable. (BoundaryManager_2D::iniBoundarySetPEdges)\n");
		return;
	}
	//初始化vBoundarySets的每个set的pEdges，存储各edge的指针
	for (int iset = 0; iset < set_edge_ID.size(); iset++) {
		std::vector<int>& edgeIDs = set_edge_ID[iset];//第iset条set的edgeIDs
		for (int iw = 0; iw < edgeIDs.size(); iw++) {//第iset条set的第iw个edge
			//在f->pEdgeTable中，根据edgeID查询pEdge，完成初始化
			int edge_ID = edgeIDs[iw];
			Edge_2D* pE = f->pEdgeTable[edge_ID];
			vBoundarySets[iset].pEdges.push_back(pE);
		}
	}
}

std::vector<int> BoundaryManager_2D::splitInts(const std::vector<int>& ints) {
	std::vector<int> ret;
	if (ints.size() < 2) {//空，或者只有1个元素
		if (ints.size() == 1) {
			ret.push_back(ints[0]);
			ret.push_back(ints[0]);
		}
		return ret;//空-返回空。仅1个元素a-返回a,a
	}
	ret.push_back(ints[0]);
	int last, current;
	last = ints[0];
	for (int i = 1; i < ints.size(); i++) {
		current = ints[i];
		if (current == last + 1) {
			if (i == (int)ints.size() - 1)ret.push_back(ints[i]);
		}
		else {//间断
			ret.push_back(last);
			ret.push_back(current);
			if (i == (int)ints.size() - 1)
				ret.push_back(current);//结尾单独一个，则重复2次
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
	//函数目的：根据BoundarySet的name得到type；给边界edge打上setID标签
	
	//检查f->edges是否已初始化
	if (f->edges.size() == 0) {
		std::string str = "Warning: f->edges.size() == 0, in \"BoundaryManager_2D::iniEdgeBoundaryCondition(FVM_2D* f)\"\n";
		LogWriter::writeLog(str, 1);
		std::cout<< str;
		return -1;
	}

	//打标签
	for (int is = 0; is < vBoundarySets.size(); is++) {
		//边界名称转边界类型
		int bType = -1;//边界类型
		const std::string bName = vBoundarySets[is].name;//边界名称
		if (bName == "wall" || bName == "obstacle") {
			if (GlobalPara::physicsModel::equation == _EQ_euler)bType = _BC_wall_nonViscous;
			else bType = _BC_wall_adiabat;
		}
		else if (bName == "inlet")bType = _BC_inlet;
		else if (bName == "outlet")bType = _BC_outlet;
		else if (bName == "inf")bType = _BC_inf;
		else if (bName == "symmetry")bType = _BC_symmetry;
		else if (bName.substr(0, 8) == "periodic") {//periodic_x
			int x = std::stoi(bName.substr(bName.size() - 1, 1));//取最后一个(第size-1个)字符
			bType = _BC_periodic_0 + x;//_BC_periodic_x = _BC_periodic_0 + x
		}

		//边界类型赋给set的type，边界ID赋给edge的setID
		vBoundarySets[is].type = bType;
		for (int ie = 0; ie < vBoundarySets[is].pEdges.size(); ie++) {
			vBoundarySets[is].pEdges[ie]->setID = vBoundarySets[is].ID;
		}

		//初始化periodPairs，periodPairs存储周期边界的配对信息，用来检查周期边界完整性
		if (_BC_periodic_0 <= bType && bType <= _BC_periodic_9) {//!注意不能连写成_BC_periodic_0 <= bType <= _BC_periodic_9
			//检查periodPairs是否已经录入该bType，若是，则使用已有的；若否，则新建一个PeriodPair tmp，存入periodPairs
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

	//检查周期边界的正确性和完整性
	//gmsh输出的inp文件中，edge编号是从第一个顶点出发，逆时针环绕一圈。因此一对periodPair其方向必然是相反的
	int check_2 = 1;
	for (int ip = 0; ip < periodPairs.size() && check_2; ip++) {
		//某一对边界
		const VirtualBoundarySet_2D& vbs0 = vBoundarySets[periodPairs[ip].setID_0 - 1];
		const VirtualBoundarySet_2D& vbs1 = vBoundarySets[periodPairs[ip].setID_1 - 1];
		//检查对应点的坐标是否满足平移规律
		if (vbs0.pEdges.size() != vbs1.pEdges.size()) {//首先应满足点个数相同
			check_2 *= 0;
		}
		else {
			const int eSize = int(vbs0.pEdges.size());
			//计算平移向量
			double vx, vy;//平移向量v
			vx = vbs1.pEdges[eSize - 1]->getx() - vbs0.pEdges[0]->getx();
			vy = vbs1.pEdges[eSize - 1]->gety() - vbs0.pEdges[0]->gety();
			for (int ie = 0; ie < eSize && check_2; ie++) {
				//将边界1中某点平移(vx,vy)得到p0,判断p0和p1是否重合
				double p0_x = vbs0.pEdges[ie]->getx() + vx;
				double p0_y = vbs0.pEdges[ie]->gety() + vy;
				double p1_x = vbs1.pEdges[eSize - 1 - ie]->getx();
				double p1_y = vbs1.pEdges[eSize - 1 - ie]->gety();
				double distance2 = (p1_x - p0_x) * (p1_x - p0_x) + (p1_y - p0_y) * (p1_y - p0_y);
				double vlength2 = vx * vx + vy * vy;
				if (distance2 >= 0.0001 * vlength2)check_2 *= 0;//|p1-p0|大于0.01*|v|则认为不重合
			}
		}
		//不满足则报错
		if (!check_2) {
			std::string str = "Error: invalid periodic boundary, for not meeting translation conditions. (BoundaryManager_2D::iniBoundaryEdge_SetType)\n";
			LogWriter::writeLogAndCout(str);
			return -1;
		}
	}
	return 0;
}

void BoundaryManager_2D::ini_infBoundary_ruvp() {
	//初始化远场边界的ruvp
	//以下为debug，设定特定ruvp
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

	//初始内部条件见FVM_2D::initField
}

void BoundaryManager_2D::setBoundaryElementU(int tag) {
	//[debug]远场边界单元守恒量强制赋值
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
	if (setID <= 0 || setID > vBoundarySets.size()) {//setID范围应从1开始，最大为vBoundarySets的size
		LogWriter::writeLogAndCout("Error: setID out of range. (getBoundarySetByID)\n");
		return nullptr;
	}
	return &(vBoundarySets[setID - 1]);
}

VirtualBoundarySet_2D* BoundaryManager_2D::getPairByID_periodicBoundary(const int setID) {
	VirtualBoundarySet_2D* set_0 = getBoundarySetByID(setID);
	int setID_to_find = -1;//待寻找的set
	const int bType = set_0->type;
	if (_BC_periodic_0 <= bType && bType <= _BC_periodic_9) {//判断是否为周期边界
		//是周期边界，则寻找配对
		for (int ip = 0; ip < periodPairs.size(); ip++) {
			if (periodPairs[ip].bType == bType) {//若是同一对周期边界
				if (periodPairs[ip].setID_0 == setID)setID_to_find = periodPairs[ip].setID_1;
				else setID_to_find = periodPairs[ip].setID_0;
				break;//跳出for循环
			}
		}
		//若找到配对的ID，则返回指针，否则报错
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
	//函数功能：找到某edge的虚拟pElement_R。仅限于周期边界

	//找到pEdge对应的boundarySet(pBoundary_0)，以及配对boundarySet(pBoundary_1)
	VirtualBoundarySet_2D* pBoundary_0 = FVM_2D::pFVM2D->boundaryManager.getBoundarySetByID(pEdge_0->setID);//pEdge对应的boundarySet
	VirtualBoundarySet_2D* pBoundary_1 = FVM_2D::pFVM2D->boundaryManager.getPairByID_periodicBoundary(pEdge_0->setID);
	//获取pEdge在pBoundary_0的pEdges中的序号(从0开始的指标)
	int index_0 = pBoundary_0->get_pEdge_index(pEdge_0);
	if (index_0 == -1) {
		LogWriter::writeLogAndCout("Error: edge not found. (BoundaryManager_2D::get_pElement_R_periodic)\n");
		return nullptr;
	}
	//计算得到它在pBoundary_1的序号
	const int eSize = int(pBoundary_0->pEdges.size());
	int index_1 = eSize - 1 - index_0;
	//找到pEdge_0在pBoundary_1上对应的pEdge_1
	Edge_2D* pEdge_1 = pBoundary_1->pEdges[index_1];
	//pEdge_1的pElement_L即为pEdge_0的虚拟pElement_R
	return pEdge_1->pElement_L;
}

Edge_2D* BoundaryManager_2D::get_pairEdge_periodic(Edge_2D* pEdge_0) {
	//函数功能：找到某edge对应的edge。仅限于周期边界
		//找到pEdge对应的boundarySet(pBoundary_0)，以及配对boundarySet(pBoundary_1)
	VirtualBoundarySet_2D* pBoundary_0 = FVM_2D::pFVM2D->boundaryManager.getBoundarySetByID(pEdge_0->setID);//pEdge对应的boundarySet
	VirtualBoundarySet_2D* pBoundary_1 = FVM_2D::pFVM2D->boundaryManager.getPairByID_periodicBoundary(pEdge_0->setID);
	//获取pEdge在pBoundary_0的pEdges中的序号(从0开始的指标)
	int index_0 = pBoundary_0->get_pEdge_index(pEdge_0);
	if (index_0 == -1) {
		LogWriter::writeLogAndCout("Error: edge not found. (BoundaryManager_2D::get_pElement_R_periodic)\n");
		return nullptr;
	}
	//计算得到它在pBoundary_1的序号
	const int eSize = int(pBoundary_0->pEdges.size());
	int index_1 = eSize - 1 - index_0;
	//找到pEdge_0在pBoundary_1上对应的pEdge_1
	return pBoundary_1->pEdges[index_1];
}

int VirtualBoundarySet_2D::get_pEdge_index(Edge_2D* pEdge) {
	//查询pEdge是否在pEdges中，若是，则返回指标(从0开始)，否则返回-1
	for (int i = 0; i < pEdges.size(); i++) {
		if (pEdge == pEdges[i])return i;
	}
	return -1;
}
