#include "OldDataConverter.h"

void U2NITS::OldDataConverter::Convert_FVM2D_to_HostData() {
	const int num_node = (int)m_pFVM2D->nodes.size();
	const int num_element = (int)m_pFVM2D->elements.size();
	const int num_edge = (int)m_pFVM2D->edges.size();
	const int num_boundary = (int)m_pFVM2D->boundaryManager.boundaries.size();

	// ��ʼ��Element_2D vector�е�GPUindex
	InitializeOldGPUIndex(num_node, num_element, num_edge);

	// ��FVM_2D���ݳ�ʼ��host
	ConvertNode(num_node);
	ConvertElement(num_element);
	ConvertEdge(num_edge);
	ConvertBoundary(num_boundary);


}

void U2NITS::OldDataConverter::InitializeOldGPUIndex(int num_node, int num_element, int num_edge) {
	// ��ʼ��Element_2D vector�е�GPUindex

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
	������������
	*/
	//const int num_element = (int)m_pFVM2D->elements.size();
	GPU::EdgeSoA& edge_host = m_pGPUSolver2->edge_host;
	GPU::NodeSoA& node_host = m_pGPUSolver2->node_host;
	GPU::ElementSoA& element_host = m_pGPUSolver2->element_host;
	GPU::ElementFieldSoA& elementField_host = m_pGPUSolver2->elementField_host;

	// ��ʼ��element_host��elementAdjacent
	FVM_2D* pFVM2D = (FVM_2D*)m_pFVM2D;
	const int nodePerElement = 3;
	const int nValue = 4;
	//#pragma omp parallel for

	for (int i = 0; i < num_element; i++) {

		Element_2D& element_i = pFVM2D->elements[i];
		// ID xy volume
		// ע��volume���õ���ǰ�ĺ������㣬���Ҫ���¶�ȡ�Ļ���xy volume��Ҫ���ں��棬��Ϊ
		// node��δȫ����ʼ��
		element_host.ID[i] = element_i.GPUID;
		element_host.xy[0][i] = element_i.x;
		element_host.xy[1][i] = element_i.y;
		element_host.volume[i] = element_i.area;// 20240405����������װedgeʱ����������������ڵ�˳��Ϊ��ʱ��
		for (int j = 0; j < nodePerElement; j++) {
			element_host.nodes[j][i] = pFVM2D->getNodeByID(element_i.nodes[j])->GPUID;
			element_host.edges[j][i] = element_i.pEdges[j]->GPUID;
		}
		element_host.nodes[3][i] = -1;// �������ε�Ԫ�����4������ͱ߲����ڣ�IDȡ-1
		element_host.edges[3][i] = -1;
		for (int j = 0; j < nValue; j++) {
			elementField_host.U[j][i] = element_i.U[j];
			elementField_host.Ux[j][i] = element_i.Ux[j];
			elementField_host.Uy[j][i] = element_i.Uy[j];
			elementField_host.Flux[j][i] = element_i.Flux[j];
		}
		// �ھӣ���edges˳�����У�-1��ʾ���ھ�
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
	// ��ʼ��edge_host
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

		// ����߽����꣬���ڽڵ�����ƽ��ֵ
		int tmpNodeID0 = edge_host.nodes[0][i];// node ID
		int tmpNodeID1 = edge_host.nodes[1][i];
		edge_host.xy[0][i] = (node_host.xy[0][tmpNodeID0] + node_host.xy[0][tmpNodeID1]) / 2.0;
		edge_host.xy[1][i] = (node_host.xy[1][tmpNodeID0] + node_host.xy[1][tmpNodeID1]) / 2.0;
		// ����normal
		//edge_i.getDirectionN(edge_host.normal[0][i], edge_host.normal[1][i]);
		myfloat dx = node_host.xy[0][tmpNodeID1] - node_host.xy[0][tmpNodeID0];
		myfloat dy = node_host.xy[1][tmpNodeID1] - node_host.xy[1][tmpNodeID0];
		edge_host.normal[0][i] = dy / edge_i.length;
		edge_host.normal[1][i] = -dx / edge_i.length;
	}

}

void U2NITS::OldDataConverter::ConvertBoundary(int num_boundary) {
	/*
	��ʼ������edge��elementR
	��ʼ������edge��edge_pair
	��ʼ���߽�ID����ӳ�����m_pGPUSolver2->boundary_host��<boundary set ID, boundary type>��ӳ���ϵ
	*/
	/**
	* ���ڱ߽綨λ��Ӧ��(edge pair)�ķ���
	* ��1����EdgeSoA������Աint periodicPair���ŵ��ǲ��ҿ죬ȱ����ռ�ÿռ��
	* ��2����map��hashTable�洢������GPU�����ݽṹ��device������Ҫ�Լ�д
	*
	* ĿǰGPU���÷�1��CPU���÷�2(?)
	*/
	GPU::EdgeSoA& edge_host = m_pGPUSolver2->edge_host;
	BoundaryManager_2D& boundaryManager = m_pFVM2D->boundaryManager;

	// ��ʼ������edge��elementR��edge_pair
	for (VirtualBoundarySet_2D& boundary : boundaryManager.boundaries) {
		int bType = boundary.type;
		if (_BC_periodic_0 <= bType && bType <= _BC_periodic_9) {
			// ��ĳ�߽缯set�����ڱ߽磬������edge
			for (Edge_2D* pEdge : boundary.pEdges) {
				// ��ȡedge��edge pair��ID
				int edgeID = pEdge->GPUID;
				int edgeID_pair = boundaryManager.get_pairEdge_periodic(pEdge)->GPUID;
				// ��pair��elementL����edge��elementR
				int elementL_pair = edge_host.elementL[edgeID_pair];
				edge_host.elementR[edgeID] = elementL_pair;
				// ��map�洢pair
				m_pGPUSolver2->edge_periodic_pair.insert(std::pair<int, int>(edgeID, edgeID_pair));
				// ��EdgeSoA��Ա�洢pair
				edge_host.periodicPair[edgeID] = edgeID_pair;
			}
		}
	}

	// ��ʼ���߽�ID����ӳ���boundary_host
	for (int i = 0; i < num_boundary; i++) {
		m_pGPUSolver2->boundary_host.type[i] = boundaryManager.boundaries[i].type;
	}

	// ��ʼ��boundary_host_new��
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
