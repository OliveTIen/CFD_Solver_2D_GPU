#include "OldDataConverter.h"

void U2NITS::OldDataConverter::Convert_FVM2D_to_HostData() {
	const int num_node = (int)m_pFVM2D->nodes.size();
	const int num_element = (int)m_pFVM2D->elements.size();
	const int num_edge = (int)m_pFVM2D->edges.size();
	const int num_boundary = (int)m_pFVM2D->boundaryManager.boundaries.size();

	// 初始化Element_2D vector中的GPUindex
	InitializeOldGPUIndex(num_node, num_element, num_edge);

	// 用FVM_2D数据初始化host
	ConvertNode(num_node);
	ConvertElement(num_element);
	ConvertEdge(num_edge);
	ConvertBoundary(num_boundary);


}

void U2NITS::OldDataConverter::InitializeOldGPUIndex(int num_node, int num_element, int num_edge) {
	// 初始化Element_2D vector中的GPUindex

	for (int i = 0; i < num_node; i++) {
		m_pFVM2D->nodes[i].GPUID = i;
	}
	for (int i = 0; i < num_element; i++) {
		m_pFVM2D->elements[i].GPUID = i;
	}
	for (int i = 0; i < num_edge; i++) {
		m_pFVM2D->edges[i].GPUID = i;
	}
}

void U2NITS::OldDataConverter::ConvertNode(int num_node) {
	//const int num_node = (int)m_pFVM2D->nodes.size();
	for (int i = 0; i < num_node; i++) {
		Node_2D node_i = m_pFVM2D->nodes[i];
		m_pGPUSolver2->node_host.ID[i] = node_i.GPUID;
		m_pGPUSolver2->node_host.xy[0][i] = node_i.x;
		m_pGPUSolver2->node_host.xy[1][i] = node_i.y;
	}
}

void U2NITS::OldDataConverter::ConvertElement(int num_element) {
	/*
	假设是三角形
	*/
	//const int num_element = (int)m_pFVM2D->elements.size();
	GPU::EdgeSoA& edge_host = m_pGPUSolver2->edge_host;
	GPU::NodeSoA& node_host = m_pGPUSolver2->node_host;
	GPU::ElementSoA& element_host = m_pGPUSolver2->element_host;
	GPU::ElementFieldSoA& elementField_host = m_pGPUSolver2->elementField_host;

	// 初始化element_host和elementAdjacent
	FVM_2D* pFVM2D = (FVM_2D*)m_pFVM2D;
	const int nodePerElement = 3;
	const int nValue = 4;
	//#pragma omp parallel for

	for (int i = 0; i < num_element; i++) {

		Element_2D& element_i = pFVM2D->elements[i];
		// ID xy volume
		// 注意volume是用的以前的函数计算，如果要重新读取的话，xy volume需要放在后面，因为
		// node尚未全部初始化
		element_host.ID[i] = element_i.GPUID;
		element_host.xy[0][i] = element_i.x;
		element_host.xy[1][i] = element_i.y;
		element_host.volume[i] = element_i.area;// 20240405新增：在组装edge时计算面积，并调整节点顺序为逆时针
		for (int j = 0; j < nodePerElement; j++) {
			element_host.nodes[j][i] = pFVM2D->getNodeByID(element_i.nodes[j])->GPUID;
			element_host.edges[j][i] = element_i.pEdges[j]->GPUID;
		}
		element_host.nodes[3][i] = -1;// 对三角形单元，其第4个顶点和边不存在，ID取-1
		element_host.edges[3][i] = -1;
		for (int j = 0; j < nValue; j++) {
			elementField_host.U[j][i] = element_i.U[j];
			elementField_host.Ux[j][i] = element_i.Ux[j];
			elementField_host.Uy[j][i] = element_i.Uy[j];
			elementField_host.Flux[j][i] = element_i.Flux[j];
		}
		// 邻居，按edges顺序排列，-1表示无邻居
		for (int j = 0; j < 4; j++) {
			element_host.neighbors[j][i] = -1;
		}
		std::vector<Element_2D*> neighbors_element_i = element_i.findNeighbor();
		size_t num_neighbor = neighbors_element_i.size();// 3
		for (int j = 0; j < num_neighbor; j++) {
			if (neighbors_element_i[j] != nullptr) {
				element_host.neighbors[j][i] = neighbors_element_i[j]->GPUID;
			}
		}
	}


}

void U2NITS::OldDataConverter::ConvertEdge(int num_edge) {
	//const int num_edge = (int)m_pFVM2D->edges.size();
	GPU::EdgeSoA& edge_host = m_pGPUSolver2->edge_host;
	GPU::NodeSoA& node_host = m_pGPUSolver2->node_host;
	// 初始化edge_host
	FVM_2D* pFVM2D = (FVM_2D*)m_pFVM2D;
	//#pragma omp parallel for
	for (int i = 0; i < num_edge; i++) {
		Edge_2D& edge_i = pFVM2D->edges[i];
		edge_host.ID[i] = edge_i.GPUID;
		edge_host.nodes[0][i] = pFVM2D->getNodeByID(edge_i.nodes[0])->GPUID;
		edge_host.nodes[1][i] = pFVM2D->getNodeByID(edge_i.nodes[1])->GPUID;
		edge_host.setID[i] = edge_i.setID;
		edge_host.periodicPair[i] = -1;
		edge_host.elementL[i] = edge_i.pElement_L->GPUID;
		if (edge_i.pElement_R == nullptr) {
			edge_host.elementR[i] = -1;
		}
		else {
			edge_host.elementR[i] = edge_i.pElement_R->GPUID;
		}
		edge_host.length[i] = edge_i.length;
		edge_host.distanceOfElements[i] = edge_i.refLength;

		// 计算边界坐标，等于节点坐标平均值
		int tmpNodeID0 = edge_host.nodes[0][i];// node ID
		int tmpNodeID1 = edge_host.nodes[1][i];
		edge_host.xy[0][i] = (node_host.xy[0][tmpNodeID0] + node_host.xy[0][tmpNodeID1]) / 2.0;
		edge_host.xy[1][i] = (node_host.xy[1][tmpNodeID0] + node_host.xy[1][tmpNodeID1]) / 2.0;
		// 计算normal
		//edge_i.getDirectionN(edge_host.normal[0][i], edge_host.normal[1][i]);
		myfloat dx = node_host.xy[0][tmpNodeID1] - node_host.xy[0][tmpNodeID0];
		myfloat dy = node_host.xy[1][tmpNodeID1] - node_host.xy[1][tmpNodeID0];
		edge_host.normal[0][i] = dy / edge_i.length;
		edge_host.normal[1][i] = -dx / edge_i.length;
	}

}

void U2NITS::OldDataConverter::ConvertBoundary(int num_boundary) {
	/*
	初始化周期edge的elementR
	初始化周期edge的edge_pair
	初始化边界ID类型映射表，即m_pGPUSolver2->boundary_host的<boundary set ID, boundary type>的映射关系
	*/
	/**
	* 周期边界定位对应边(edge pair)的方法
	* 法1：在EdgeSoA新增成员int periodicPair，优点是查找快，缺点是占用空间大。
	* 法2：用map或hashTable存储，但是GPU的数据结构和device函数需要自己写
	*
	* 目前GPU采用法1，CPU采用法2(?)
	*/
	GPU::EdgeSoA& edge_host = m_pGPUSolver2->edge_host;
	BoundaryManager_2D& boundaryManager = m_pFVM2D->boundaryManager;

	// 初始化周期edge的elementR和edge_pair
	for (VirtualBoundarySet_2D& boundary : boundaryManager.boundaries) {
		int bType = boundary.type;
		if (_BC_periodic_0 <= bType && bType <= _BC_periodic_9) {
			// 若某边界集set是周期边界，遍历其edge
			for (Edge_2D* pEdge : boundary.pEdges) {
				// 获取edge和edge pair的ID
				int edgeID = pEdge->GPUID;
				int edgeID_pair = boundaryManager.get_pairEdge_periodic(pEdge)->GPUID;
				// 用pair的elementL更新edge的elementR
				int elementL_pair = edge_host.elementL[edgeID_pair];
				edge_host.elementR[edgeID] = elementL_pair;
				// 用map存储pair
				m_pGPUSolver2->edge_periodic_pair.insert(std::pair<int, int>(edgeID, edgeID_pair));
				// 用EdgeSoA成员存储pair
				edge_host.periodicPair[edgeID] = edgeID_pair;
			}
		}
	}

	// 初始化边界ID类型映射表boundary_host
	for (int i = 0; i < num_boundary; i++) {
		m_pGPUSolver2->boundary_host.type[i] = boundaryManager.boundaries[i].type;
	}

	// 初始化boundary_host_new。
	m_pGPUSolver2->boundary_host_new.edgeSets.resize(boundaryManager.boundaries.size());
	for (int i = 0; i < num_boundary; i++) {
		
		const auto& boundary = boundaryManager.boundaries[i];
		auto& edgeSet = m_pGPUSolver2->boundary_host_new.edgeSets[i];
		edgeSet.ID = boundary.ID;
		edgeSet.name = boundary.name;
		edgeSet.type = boundary.type;
		edgeSet.edge_vector.resize(boundary.pEdges.size());
		for (int ie = 0; ie < boundary.pEdges.size(); ie++) {
			edgeSet.edge_vector[ie] = boundary.pEdges[ie]->GPUID;
		}
		
	}

}
