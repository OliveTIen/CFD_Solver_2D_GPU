#include "FVM_2D.h"
#include "output/LogWriter.h"


std::vector<int> BoundaryManager_2D::compressSeveralSequences(const std::vector<int>& ints) {
	/*
	该函数用于压缩存储数列，用来减少占用空间但不损失信息
	
	输入参数ints是分段连续的数列，输出是该数列各段的左右端点值

	例如，"1-5,11-14,21-26"->"1,5,11,14,21,26"
	输入：1,2,3,4,5,11,12,13,14,21,22,23,24,25,26
	输出：1,5,11,14,21,26

	实际应用中，数列的每个元素表示一个边界线元的ID。一条连续数列对应一个边界

	*/

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

int BoundaryManager_2D::getBoundaryTypeByName(std::string boundaryName) {
	// 边界名称转边界类型
	int bType = -1;
	if (boundaryName == "wall" || boundaryName == "obstacle") {
		if (GlobalPara::physicsModel::equation == _EQ_euler)bType = _BC_wall_nonViscous;
		else bType = _BC_wall_adiabat;
	}
	else if (boundaryName == "inlet")bType = _BC_inlet;
	else if (boundaryName == "outlet")bType = _BC_outlet;
	else if (boundaryName == "inf")bType = _BC_inf;
	else if (boundaryName == "symmetry")bType = _BC_symmetry;
	else if (boundaryName.substr(0, 8) == "periodic") {//periodic_x
		int x = std::stoi(boundaryName.substr(boundaryName.size() - 1, 1));//取最后一个字符
		bType = _BC_periodic_0 + x;//_BC_periodic_x = _BC_periodic_0 + x
	}
	return bType;
}

VirtualBoundarySet_2D* BoundaryManager_2D::findSetByID(int ID) {
	if (ID < 1 || ID > boundaries.size()) {
		std::cout << "Error: illegal ID, in \"BoundaryManager_2D::findSetByID(int ID)\"\n";
		return nullptr;
	}
	return &(boundaries[ID - 1]);
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
	for (int is = 0; is < boundaries.size(); is++) {
		//边界类型赋给set的type，边界ID赋给edge的setID
		int bType = getBoundaryTypeByName(boundaries[is].name);
		boundaries[is].type = bType;
		for (int ie = 0; ie < boundaries[is].pEdges.size(); ie++) {
			boundaries[is].pEdges[ie]->setID = boundaries[is].ID;
		}

		//! 注意该段代码存在于多处。修改时可能需要修改多处
		// 初始化periodPairs，periodPairs存储周期边界的配对信息，用来检查周期边界完整性
		// periodPairs是boundaryManager的成员，每个pair存储int bType, int setID_0, int setID_1
		if (_BC_periodic_0 <= bType && bType <= _BC_periodic_9) {// 两个不等式需拆开，用&&连接
			//检查periodPairs是否已经录入该bType，若是，则使用已有的；若否，则新建一个PeriodPair tmp，存入periodPairs
			int index_pairs = -1;
			for (int i = 0; i < periodPairs.size(); i++) {
				if (periodPairs[i].bType == bType)index_pairs = i;
			}
			// 存在，则直接设置
			if (index_pairs != -1) {
				periodPairs[index_pairs].setID_1 = boundaries[is].ID;
			}
			// 不存在，则新建
			else {
				PeriodPair tmp;
				tmp.bType = bType;
				tmp.setID_0 = boundaries[is].ID;
				periodPairs.push_back(tmp);
			}
		}
	}

	//检查周期边界的正确性和完整性
	checkPeriodPairs();

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
			const int bType = boundaries[edge_tmp.setID - 1].type;
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
			const int bType = boundaries[edge_tmp.setID - 1].type;
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
	if (setID <= 0 || setID > boundaries.size()) {//setID范围应从1开始，最大为vBoundarySets的size
		LogWriter::writeLogAndCout("Error: setID out of range. (getBoundarySetByID)\n");
		return nullptr;
	}
	return &(boundaries[setID - 1]);
}

void BoundaryManager_2D::checkPeriodPairs() {
	// 检查周期边界的正确性和完整性
	// !隐患：存在以下假设：
	// gmsh输出的inp文件中，edge编号是从第一个顶点出发，逆时针环绕一圈。因此一对periodPair其方向必然是相反的
	int check_2 = 1;
	for (int ip = 0; ip < periodPairs.size() && check_2; ip++) {
		//某一对边界
		const VirtualBoundarySet_2D& vbs0 = boundaries[periodPairs[ip].setID_0 - 1];
		const VirtualBoundarySet_2D& vbs1 = boundaries[periodPairs[ip].setID_1 - 1];
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
			throw str;
		}
	}
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

Element_2D* BoundaryManager_2D::get_pElement_R_periodic(Edge_2D* pEdge_0) {
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
	// !该函数具有隐患，它默认index_1 = eSize - 1 - index_0，即边界节点ID是逆时针增大的
	// 函数功能：找到某edge对应的edge。仅限于周期边界
	// 
	
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
