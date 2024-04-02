#include "Gradient.h"
#include "../../math/Math.h"
#include "../../output/LogWriter.h"
#include "../../global/GlobalPara.h"

void U2NITS::Space::Gradient::Gradient(GPU::ElementSoA& element_host, GPU::FieldSoA& elementField_host, GPU::NodeSoA& node_host, GPU::EdgeSoA& edge_host) {
	switch (GlobalPara::space::flag_gradient) {
	case _GRA_leastSquare:
		//LogWriter::logAndPrint("GradientLeastSquare_2\n");
		GradientLeastSquare_2(element_host, elementField_host, node_host);
		break;
	case _GRA_greenGauss:
		//LogWriter::logAndPrint("GradientGreenGauss\n");
		GradientGreenGauss_2(element_host, elementField_host, edge_host);
		break;

	default:
		LogWriter::logAndPrintError("invalid graident type.\n");
		exit(-1);
		break;
	}
}

void U2NITS::Space::Gradient::GradientLeastSquare_old(GPU::ElementSoA& element_host, GPU::FieldSoA& elementField_host) {
	// ����device������Ϊhost

	for (int i = 0; i < element_host.num_element; i++) {
		//const int bid = blockIdx.x;
		//const int tid = threadIdx.x;
		//const int id = bid * blockDim.x + tid;
		//const int& i = id;
		//if (i >= element_device.num_element || i < 0) return;

		// ��ȡ��Ч�ھ�
		int nValidNeighbor = 0;
		int validNeighbors[3]{ -1,-1,-1 };// ��Ч�ھӵ�ID
		const int neighborArraySize = 3;// �����ε�Ԫ�����3���ھ�
		for (int j = 0; j < neighborArraySize; j++) {
			if (element_host.neighbors[j][i] != -1) {
				validNeighbors[nValidNeighbor] = element_host.neighbors[j][i];
				nValidNeighbor += 1;
			}
		}

		// �����ݶ�dUdX
		const int nVar = 4;// �غ�����������άΪ4
		const int nX = 2;// �����������άΪ2
		//const int nVNmax = 3;// ����ھӸ���
		REAL dUdX[nX][nVar]{};// �ѳ�ʼ��Ϊ0
		if (nValidNeighbor == 3) {
			// ��ʼ��dX��dU
			const int nVN = 3;// ��Ч�ھӸ���

			REAL dX[nVN][nX]{};
			REAL dU[nVN][nVar]{};
			for (int iVN = 0; iVN < nVN; iVN++) {
				int index = validNeighbors[iVN];// �ھӵ�ID
				dX[iVN][0] = element_host.xy[0][index] - element_host.xy[0][i];
				dX[iVN][1] = element_host.xy[1][index] - element_host.xy[1][i];
				for (int iVar = 0; iVar < nVar; iVar++) {
					//dU[iVN][iVar] = elementField_host.U[index][iVar] - elementField_host.U[i][iVar];
					dU[iVN][iVar] = elementField_host.U[iVar][index] - elementField_host.U[iVar][i];
				}
			}

			// ��С���˷� x=(A'A)^{-1}A'b��dUdX = inv(dXtrans * dX) * dXtrans * dU
			// ! ���º���Ϊdevice��Ӧ��GPUʵ��
			REAL dXtrans[nX][nVN]{};
			REAL invdXdX[nX][nX]{};
			REAL invdXdX_dX[nX][nVN]{};
			U2NITS::Math::Matrix::transpose(nVN, nX, (REAL*)dX, (REAL*)dXtrans);
			U2NITS::Math::Matrix::mul_ixj_jxk(nX, nVN, nX, (REAL*)dXtrans, (REAL*)dX, (REAL*)invdXdX);
			U2NITS::Math::Matrix::inv_2x2(invdXdX);
			U2NITS::Math::Matrix::mul_ixj_jxk(nX, nX, nVN, (REAL*)invdXdX, (REAL*)dXtrans, (REAL*)invdXdX_dX);
			U2NITS::Math::Matrix::mul_ixj_jxk(nX, nVN, nVar, (REAL*)invdXdX_dX, (REAL*)dU, (REAL*)dUdX);
		}
		else if (nValidNeighbor == 2) {
			// ���´���������Ĵ�����ͬ(����nVN��ͬ)���޸�ʱ��һͬ�޸�
			// ��ʼ��dX��dU
			const int nVN = 2;// ��Ч�ھӸ���
			REAL dX[nVN][nX]{};
			REAL dU[nVN][nVar]{};
			for (int iVN = 0; iVN < nVN; iVN++) {
				int index = validNeighbors[iVN];// �ھӵ�ID
				dX[iVN][0] = element_host.xy[0][index] - element_host.xy[0][i];
				dX[iVN][1] = element_host.xy[1][index] - element_host.xy[1][i];
				for (int iVar = 0; iVar < nVar; iVar++) {
					dU[iVN][iVar] = elementField_host.U[iVar][index] - elementField_host.U[iVar][i];
				}
			}

			// ��С���˷� x=(A'A)^{-1}A'b��dUdX = inv(dXtrans * dX) * dXtrans * dU
			// ! ���º���Ϊdevice��Ӧ��GPUʵ��
			REAL dXtrans[nX][nVN]{};
			REAL invdXdX[nX][nX]{};
			REAL invdXdX_dX[nX][nVN]{};
			U2NITS::Math::Matrix::transpose(nVN, nX, (REAL*)dX, (REAL*)dXtrans);
			U2NITS::Math::Matrix::mul_ixj_jxk(nX, nVN, nX, (REAL*)dXtrans, (REAL*)dX, (REAL*)invdXdX);
			U2NITS::Math::Matrix::inv_2x2(invdXdX);
			U2NITS::Math::Matrix::mul_ixj_jxk(nX, nX, nVN, (REAL*)invdXdX, (REAL*)dXtrans, (REAL*)invdXdX_dX);
			U2NITS::Math::Matrix::mul_ixj_jxk(nX, nVN, nVar, (REAL*)invdXdX_dX, (REAL*)dU, (REAL*)dUdX);
		}
		else if (nValidNeighbor == 1) {
			// ���´���������Ĵ�����ͬ(����nVN��ͬ)���޸�ʱ��һͬ�޸�
			// ��ʼ��dX��dU
			const int nVN = 1;// ��Ч�ھӸ���
			REAL dX[nVN][nX]{};
			REAL dU[nVN][nVar]{};
			for (int iVN = 0; iVN < nVN; iVN++) {
				int index = validNeighbors[iVN];// �ھӵ�ID
				dX[iVN][0] = element_host.xy[0][index] - element_host.xy[0][i];
				dX[iVN][1] = element_host.xy[1][index] - element_host.xy[1][i];
				for (int iVar = 0; iVar < nVar; iVar++) {
					dU[iVN][iVar] = elementField_host.U[iVar][index] - elementField_host.U[iVar][i];
				}
			}

			// ��С���˷� x=(A'A)^{-1}A'b��dUdX = inv(dXtrans * dX) * dXtrans * dU
			// ! ���º���Ϊdevice��Ӧ��GPUʵ��
			REAL dXtrans[nX][nVN]{};
			REAL invdXdX[nX][nX]{};
			REAL invdXdX_dX[nX][nVN]{};
			U2NITS::Math::Matrix::transpose(nVN, nX, (REAL*)dX, (REAL*)dXtrans);
			U2NITS::Math::Matrix::mul_ixj_jxk(nX, nVN, nX, (REAL*)dXtrans, (REAL*)dX, (REAL*)invdXdX);
			U2NITS::Math::Matrix::inv_2x2(invdXdX);
			U2NITS::Math::Matrix::mul_ixj_jxk(nX, nX, nVN, (REAL*)invdXdX, (REAL*)dXtrans, (REAL*)invdXdX_dX);
			U2NITS::Math::Matrix::mul_ixj_jxk(nX, nVN, nVar, (REAL*)invdXdX_dX, (REAL*)dU, (REAL*)dUdX);
		}

		// ��dUdX�����Ԫ
		for (int j = 0; j < nVar; j++) {
			elementField_host.Ux[j][i] = dUdX[0][j];
			elementField_host.Uy[j][i] = dUdX[1][j];
		}

		// Barth������ һ�׾��ȸ�ʽ����Ҫ��
		//pE->restructor_in_updateSlope_Barth(f);
		//Limiter::modifySlope_Barth(pE);

		// ����쳣ֵ Ӧ��������
		// ֻ�����Ƿ��С�1e10��������
	}
}

void U2NITS::Space::Gradient::GradientLeastSquare_2(GPU::ElementSoA& element_host, GPU::FieldSoA& elementField_host, GPU::NodeSoA& node_host) {

	for (int iElement = 0; iElement < element_host.num_element; iElement++) {


		// ��ȡ��Ч�ھ�
		int numOfValidNeighbor = 0;
		int validNeighbors[3]{ -1,-1,-1 };// ��Ч�ھӵ�ID
		const int neighborArraySize = 3;// �����ε�Ԫ�����3���ھ�
		for (int j = 0; j < neighborArraySize; j++) {
			if (element_host.neighbors[j][iElement] != -1) {
				validNeighbors[numOfValidNeighbor] = element_host.neighbors[j][iElement];
				numOfValidNeighbor += 1;
			}
		}
		if (numOfValidNeighbor <= 0) {
			LogWriter::logAndPrint("Error: invalid numOfValidNeighbor.\n", LogWriter::Error, LogWriter::Error);
			exit(-1);
		}
		// ���������
		const int nVar = 4;// �غ�����������άΪ4
		const int nX = 2;// �����������άΪ2
		const int& nVN = numOfValidNeighbor;// ��Ч�ھӸ�����ȡֵ1-3
		real* dX = new real[nVN * nX]{};
		real* dUdX = new real[nX * nVar]{};
		real* dU = new real[nVN * nVar]{};
		real* dXtrans = new real[nX * nVN]{};
		real* invdXdX = new real[nX * nX]{};
		real* invdXdX_dX = new real[nX * nVN]{};
		// ��ʼ��dX��dU
		for (int iVN = 0; iVN < nVN; iVN++) {
			int index = validNeighbors[iVN];// �ھӵ�ID
			dX[iVN * nX + 0] = element_host.xy[0][index] - element_host.xy[0][iElement];
			dX[iVN * nX + 1] = element_host.xy[1][index] - element_host.xy[1][iElement];
			for (int iVar = 0; iVar < nVar; iVar++) {
				dU[iVN * nVar + iVar] = elementField_host.U[iVar][index] - elementField_host.U[iVar][iElement];
			}
		}
		// ��С���˷�����dUdX��x=(A'A)^{-1}A'b��dUdX = inv(dXtrans * dX) * dXtrans * dU
		U2NITS::Math::Matrix::transpose(nVN, nX, (real*)dX, (real*)dXtrans);
		U2NITS::Math::Matrix::mul_ixj_jxk(nX, nVN, nX, (real*)dXtrans, (real*)dX, (real*)invdXdX);
		U2NITS::Math::Matrix::inv_2x2(invdXdX);
		U2NITS::Math::Matrix::mul_ixj_jxk(nX, nX, nVN, (real*)invdXdX, (real*)dXtrans, (real*)invdXdX_dX);
		U2NITS::Math::Matrix::mul_ixj_jxk(nX, nVN, nVar, (real*)invdXdX_dX, (real*)dU, (real*)dUdX);

		/*
������
һά�£�minmod��vanleer��������Ŀ���Ǳ�֤�����������������ֵ���������ڵ�Ԫ����������ֵ
�μ�����������ѧ�����������μ��ϼ�P106

���Ƶأ����Ա�֤�����ζ��㴦U���������ڵ�ԪU
���ݵ�ǰб��Ux, Uy�����㶥�����ֵdU_node��ǰ���Ѿ�����ھ��������ֵdU
����ratio = max(abs(dU_node/dU))��������1�����������ratio
		*/
		const int nNode = 3;
		real dU_node[nNode][nVar]{};
		real dX_node[nNode][nX]{};
		// ����dX_node
		for (int iNode = 0; iNode < nNode; iNode++) {
			for (int iVar = 0; iVar < nVar; iVar++) {
				int nodeID = element_host.nodes[iNode][iElement];
				dX_node[iNode][0] = node_host.xy[0][nodeID] - element_host.xy[0][iElement];
				dX_node[iNode][1] = node_host.xy[1][nodeID] - element_host.xy[1][iElement];
			}
		}
		// ����dU_node[nNode,nVar] =  dX_node[nNode,nX] * dUdX[nX,nVar]��Ȼ�����ratio
		U2NITS::Math::Matrix::mul_ixj_jxk(nNode, nX, nVar, (real*)dX_node, dUdX, (real*)dU_node);
		real ratio = 1.0;
		for (int iNode = 0; iNode < nNode; iNode++) {
			for (int iVar = 0; iVar < nVar; iVar++) {
				using namespace U2NITS::Math;
				if (dU[iNode * nVar + iVar] == 0)continue;// ����
				ratio = max(ratio, abs(dU_node[iNode][iVar] / dU[iNode * nVar + iVar]));
			}
		}
		// ratioһ�����ڵ���1.0����˿���ֱ�ӳ���ratio
		U2NITS::Math::Matrix::div_matrix_by_scalar(nX, nVar, dUdX, ratio);
		// ��dUdX�����Ԫ
		for (int jVar = 0; jVar < nVar; jVar++) {
			elementField_host.Ux[jVar][iElement] = dUdX[0 * nVar + jVar];
			elementField_host.Uy[jVar][iElement] = dUdX[1 * nVar + jVar];
		}

		delete[] dX;
		delete[] dUdX;
		delete[] dU;
		delete[] dXtrans;
		delete[] invdXdX;
		delete[] invdXdX_dX;
		
	}
}

void U2NITS::Space::Gradient::GradientGreenGauss_1(GPU::ElementSoA& element_host, GPU::FieldSoA& elementField_host, GPU::EdgeSoA& edge_host) {
	/*
	* �ɰ���ָ�˹�ݶȣ��߽�Ԫ������
	* https://zhuanlan.zhihu.com/p/370586072
	* ����ʸ���ӵ���Ԫʱ��Ӧע�����ǳ��ﻹ�ǳ��⡣
	* ע�����ͨ��ʱ�ų������ڱ߽磬��������
	* ���ڱ߽磬��ͨ��Ӧ����phiC(�μ�01-01-CFD����)
	* ���ڷǷ���(�����ε�Ԫ�ĵ�4����)��edgeID=-1, outElementID=-1, 
	*/
	// ���ݶȵ��ӵ�������
	int subIterationGradient = 3;
	int num_edge = edge_host.num_edge;

	for (int iElement = 0; iElement < element_host.num_element; iElement++) {
		// ��ǰ��ԪID����ǰ��Ԫ�����edgeID��edge�����������ⵥԪID
		int currentElementID = element_host.ID[iElement];// ���鲻Ҫֱ����iElement����ȻĿǰiElement==i(element_host.ID[i]=element_i.GPUID=i)����������
		real volumeC = element_host.volume[currentElementID];// ���
		int edgeID[4]{};
		int edgeSign[4]{};// 1��ʾ�߳��⣬-1��ʾ�߳���
		int edgeType[4]{ -1,-1,-1,-1 };
		for (int i = 0; i < 4; i++) {
			
			edgeID[i] = element_host.edges[i][iElement];
			// �Ƿ�ֵ����Ϊ-1����ֹ�������Խ�硣��������-1������
			if (edgeID[i] <= -1 || edgeID[i] >= num_edge) {
				edgeID[i] = -1;// ���������ε�Ԫ�����4���ߵ�IDֵ�ǲ�ȷ���ġ�������Χ��Ϊ-1
			}
			else if (edge_host.setID[edgeID[i]] != -1) {
				edgeType[i] = edge_host.setID[edgeID[i]];// �洢edge���ͣ�-1��ʾ�ڲ���
				edgeID[i] = -1;// ���ڲ��ߣ�Ĭ���ų��������ų��߽�edge
			}
		}
		int outElementID[4]{ -1,-1,-1,-1 };//�ⵥԪID ע�����ڱ߽�ID!=-1
		for (int i = 0; i < 4; i++) {
			int currentEdgeID = edgeID[i];
			if (currentEdgeID == -1) {
				continue;// Ĭ����-1����edgeID��-1����outElementIDһ����-1
			}
			// ȷ���ⵥԪID�ͳ��򡣶����ڲ��ߣ��޷�ȷ����elementL����elementR
			if (currentElementID != edge_host.elementR[currentEdgeID]) {
				outElementID[i] = edge_host.elementR[currentEdgeID];
				edgeSign[i] = 1;// currentElement=elementL������
			}
			else {
				outElementID[i] = edge_host.elementL[currentEdgeID];
				edgeSign[i] = -1;// currentElement=elementR������
			}
		}
		// �����ʸ��(������)
		real edgeNormal[2][4]{};
		real edgeArea[4]{};
		real edgeSVector[2][4]{};// ���ʸ��
		for (int i = 0; i < 4; i++) {
			int currentEdgeID = edgeID[i];
			if (currentEdgeID == -1) {
				continue;
			}
			edgeNormal[0][i] = edge_host.normal[0][currentEdgeID];
			edgeNormal[1][i] = edge_host.normal[1][currentEdgeID];
			edgeArea[i] = edge_host.length[currentEdgeID];
			edgeSVector[0][i] = edgeNormal[0][i] * edgeArea[i];
			edgeSVector[1][i] = edgeNormal[1][i] * edgeArea[i];
		}

		// ���ÿ���ߣ��󼸺�Ȩ������gc��ʸ��(������)
		real gC_all[4]{};
		real rf_rC_all[4][2]{};
		real rf_rF_all[4][2]{};
		for (int i = 0; i < 4; i++) {
			int currentEdgeID = edgeID[i];
			int currentOutElementID = outElementID[i];
			if (currentEdgeID == -1) {
				continue;
			}

			// ��͵�Ԫ����
			real fx = edge_host.xy[0][currentEdgeID];// face x
			real fy = edge_host.xy[1][currentEdgeID];
			real Cx = element_host.xy[0][currentElementID];// cell x
			real Cy = element_host.xy[1][currentElementID];
			real Fx = element_host.xy[0][currentOutElementID];
			real Fy = element_host.xy[1][currentOutElementID];
			using namespace U2NITS;
			real dFf = Math::distance2D(Fx, Fy, fx, fy);
			real dFC = Math::distance2D(Fx, Fy, Cx, Cy);
			if (dFC < Math::EPSILON)dFC = Math::EPSILON;

			real gC = dFf / dFC;// ����Ȩ������gc
			real rf_rC[2]{ fx - Cx,fy - Cy };// rf-rC
			real rf_rF[2]{ fx - Fx,fy - Fy };// rf-rF

			gC_all[i] = gC;
			rf_rC_all[i][0] = rf_rC[0];
			rf_rC_all[i][1] = rf_rC[1];
			rf_rF_all[i][0] = rf_rF[0];
			rf_rF_all[i][1] = rf_rF[1];
		}

		// ����ÿ���ߵ�����ͨ��phif��Sf��Ȼ����ͣ�����������õ���ԪC���ݶ�
		const int nVar = 4;
		for (int jVar = 0; jVar < nVar; jVar++) {
			// ��͡����泯�ڣ�����ͨ��ȡ�෴��
			real sum_phifSf[2]{};
			for (int i = 0; i < 4; i++) {
				int currentEdgeID = edgeID[i];
				int currentOutElementID = outElementID[i];
				if (currentEdgeID == -1) {
					continue;
				}

				// ��ȡ��Ԫֵ�͵�Ԫ�ݶ�ֵ
				real phiC = elementField_host.U[jVar][currentElementID];
				real phiF = elementField_host.U[jVar][currentOutElementID];
				//real nabla_phiC[2]{// ����ʱ�Ż��õ�
				//	elementField_host.Ux[jVar][currentElementID],
				//	elementField_host.Uy[jVar][currentElementID]
				//};// �ݶ�
				//real nabla_phiF[2]{
				//	elementField_host.Ux[jVar][currentOutElementID],
				//	elementField_host.Uy[jVar][currentOutElementID]
				//};
				// �����ֵ(����)
				real gC = gC_all[i];
				real gC1 = 1.0 - gC;
				real phif = gC * phiC + gC1 * phiF;
				real Sf[2]{ edgeSVector[0][i],edgeSVector[1][i] };
				sum_phifSf[0] += phif * Sf[0] * edgeSign[i];// ���Է��ţ���Ӧ���泯�ڵ�����
				sum_phifSf[1] += phif * Sf[1] * edgeSign[i];
			}

			real nabla_phiC_new[2]{// �ݶ�
				1 / volumeC * sum_phifSf[0],
				1 / volumeC * sum_phifSf[1]
			};
			// ���µ�Ԫ�ݶ�ֵ
			elementField_host.Ux[jVar][currentElementID] = nabla_phiC_new[0];
			elementField_host.Uy[jVar][currentElementID] = nabla_phiC_new[1];
			// ���Ҫ��������Ҫ�����е�Ԫ������ϣ���˵���ѭ��Ҫ����iElementѭ����
		}
	}
}

void U2NITS::Space::Gradient::GradientGreenGauss_2(GPU::ElementSoA& element_host, GPU::FieldSoA& elementField_host, GPU::EdgeSoA& edge_host) {
	/*
	2024-03-28

	* ���ָ�˹�ݶ�
	* https://zhuanlan.zhihu.com/p/370586072
	* ����ʸ���ӵ���Ԫʱ��Ӧע�����ǳ��ﻹ�ǳ��⡣
	* ע�����ͨ��ʱ�ų������ڱ߽磬��������
	* ���ڱ߽磬��ͨ��Ӧ����phiC(�μ�01-01-CFD����)
	* ���ڷǷ���(�����ε�Ԫ�ĵ�4����)��edgeID=-1, outElementID=-1,
	*/
// ���ݶȵ��ӵ�������
	int subIterationGradient = 3;
	int num_edge = edge_host.num_edge;

	for (int iElement = 0; iElement < element_host.num_element; iElement++) {
		// ��ǰ��ԪID����ǰ��Ԫ�����edgeID��edge�����������ⵥԪID
		int currentElementID = element_host.ID[iElement];// ���鲻Ҫֱ����iElement����ȻĿǰiElement==i(element_host.ID[i]=element_i.GPUID=i)����������
		real volumeC = element_host.volume[currentElementID];// ���
		// ��ʼ��edgeID �������鷶Χ�ķǷ�ֵ��Ϊ-1
		int edgeID[4]{ -1,-1,-1,-1 };// ��iElement��element�ĵ�i���ߵ�ID
		for (int i = 0; i < 4; i++) {
			edgeID[i] = element_host.edges[i][iElement];
			if (edgeID[i] <= -1 || edgeID[i] >= num_edge) {
				edgeID[i] = -1;
			}
		}
		// ��ʼ��edgeType
		int edgeType[4]{ -1,-1,-1,-1 };
		for (int i = 0; i < 4; i++) {
			int edgeIDi = edgeID[i];
			if (edgeIDi == -1) {
				continue;
			}

			edgeType[i] = edge_host.setID[edgeIDi];
		}
		// ��ʼ��edgeSign��outElementID
		int edgeSign[4]{ 0,0,0,0 };// 1��ʾ�߳��⣬-1��ʾ�߳���
		int outElementID[4]{ -1,-1,-1,-1 };//�ⵥԪID �ڲ��ߺ����ڱ�
		for (int i = 0; i < 4; i++) {
			int edgeIDi = edgeID[i];
			if (edgeIDi == -1) {
				continue;
			}

			if (currentElementID != edge_host.elementR[edgeIDi]) {
				edgeSign[i] = 1;// currentElement=elementL������
				outElementID[i] = edge_host.elementR[edgeIDi];
			}
			else {
				edgeSign[i] = -1;// currentElement=elementR������
				outElementID[i] = edge_host.elementL[edgeIDi];
			}
		}
		// �����ʸ��edgeSVector(������)
		real edgeNormal[2][4]{};
		real edgeArea[4]{};
		real edgeSVector[2][4]{};// ���ʸ��
		for (int i = 0; i < 4; i++) {
			int edgeIDi = edgeID[i];
			if (edgeIDi == -1) {
				continue;
			}

			edgeNormal[0][i] = edge_host.normal[0][edgeIDi];
			edgeNormal[1][i] = edge_host.normal[1][edgeIDi];
			edgeArea[i] = edge_host.length[edgeIDi];
			edgeSVector[0][i] = edgeNormal[0][i] * edgeArea[i];
			edgeSVector[1][i] = edgeNormal[1][i] * edgeArea[i];
		}

		// �󼸺�Ȩ������gc��ʸ��(������)
		real gC_all[4]{};
		real rf_rC_all[4][2]{};
		real rf_rF_all[4][2]{};
		for (int i = 0; i < 4; i++) {
			int edgeIDi = edgeID[i];
			if (edgeIDi == -1) {
				continue;
			}

			int currentOutElementID = outElementID[i];
			if (currentOutElementID != -1){
				// �ڲ��ߺ����ڱ�
				if (edgeType[i] == -1) {
					// �ڲ���

					// ��͵�Ԫ����
					real fx = edge_host.xy[0][edgeIDi];// face x
					real fy = edge_host.xy[1][edgeIDi];
					real Cx = element_host.xy[0][currentElementID];// cell x
					real Cy = element_host.xy[1][currentElementID];
					real Fx = element_host.xy[0][currentOutElementID];
					real Fy = element_host.xy[1][currentOutElementID];
					// ����Ȩ������gC = |rF-rf|/|rF-rC|���洢��gC_all������
					using namespace U2NITS;
					real dFf = Math::distance2D(Fx, Fy, fx, fy);
					real dFC = Math::distance2D(Fx, Fy, Cx, Cy);
					if (dFC < Math::EPSILON)dFC = Math::EPSILON;
					gC_all[i] = dFf / dFC;
					// ʸ�� rf-rC rf-rF���ڵ���ʱ���õ�
					rf_rC_all[i][0] = fx - Cx;// rf-rC
					rf_rC_all[i][1] = fy - Cy;
					rf_rF_all[i][0] = fx - Fx;// rf-rF
					rf_rF_all[i][1] = fy - Fy;
				}
				else {
					// ���ڱ�
					// ��͵�Ԫ����
					real fx = edge_host.xy[0][edgeIDi];
					real fy = edge_host.xy[1][edgeIDi];
					real Cx = element_host.xy[0][currentElementID];// cell x
					real Cy = element_host.xy[1][currentElementID];
					// ���ڱߵ�outElement������Ҫƽ��
					int pairIDi = edge_host.periodicPair[edgeIDi];
					real fx_pair = edge_host.xy[0][pairIDi];
					real fy_pair = edge_host.xy[1][pairIDi];
					real shiftx = fx - fx_pair;// ��pairָ��ǰλ�õ�ʸ��
					real shifty = fy - fy_pair;
					real Fx = element_host.xy[0][currentOutElementID];
					real Fy = element_host.xy[1][currentOutElementID];
					Fx += shiftx;// �ø�ʸ��ƽ�Ƶ�Ԫ
					Fy += shifty;

					// ������ڲ�����һ����
					// ����Ȩ������gC = |rF-rf|/|rF-rC|���洢��gC_all������
					using namespace U2NITS;
					real dFf = Math::distance2D(Fx, Fy, fx, fy);
					real dFC = Math::distance2D(Fx, Fy, Cx, Cy);
					if (dFC < Math::EPSILON)dFC = Math::EPSILON;
					gC_all[i] = dFf / dFC;
					// ʸ�� rf-rC rf-rF���ڵ���ʱ���õ�
					rf_rC_all[i][0] = fx - Cx;// rf-rC
					rf_rC_all[i][1] = fy - Cy;
					rf_rF_all[i][0] = fx - Fx;// rf-rF
					rf_rF_all[i][1] = fy - Fy;
				}

			}
			else {
				// �߽��(���ڱ߽����)
				// F��Ԫʵ�ʲ����ڣ�����F��Ԫ��C��Ԫ����f�����ĶԳ�
				real fx = edge_host.xy[0][edgeIDi];// face x
				real fy = edge_host.xy[1][edgeIDi];
				real Cx = element_host.xy[0][currentElementID];// cell x
				real Cy = element_host.xy[1][currentElementID];

				using namespace U2NITS;
				real gC = 0.5;// ����Ȩ������gc
				real rf_rC[2]{ fx - Cx,fy - Cy };// rf-rC
				real rf_rF[2]{ -(fx - Cx),-(fy - Cy)};// rf-rF

				gC_all[i] = gC;
				rf_rC_all[i][0] = rf_rC[0];
				rf_rC_all[i][1] = rf_rC[1];
				rf_rF_all[i][0] = rf_rF[0];
				rf_rF_all[i][1] = rf_rF[1];
			}

		}



		// ����ÿ���ߵ�����ͨ��phif��Sf��Ȼ����ͣ�����������õ���ԪC���ݶ�
		const int nVar = 4;
		for (int jVar = 0; jVar < nVar; jVar++) {
			// ��͡����泯�ڣ�����ͨ��ȡ�෴��
			real sum_phifSf[2]{};
			for (int i = 0; i < 4; i++) {
				int currentEdgeID = edgeID[i];
				if (currentEdgeID == -1) {
					continue;
				}
				int currentOutElementID = outElementID[i];

				if (currentOutElementID != -1) {
					// �ڲ��ߺ����ڱ�

					// ��ȡ��Ԫֵ�͵�Ԫ�ݶ�ֵ
					real phiC = elementField_host.U[jVar][currentElementID];
					real phiF = elementField_host.U[jVar][currentOutElementID];
					// �����ֵ(����)
					real gC = gC_all[i];
					real gC1 = 1.0 - gC;
					real phif = gC * phiC + gC1 * phiF;
					real Sf[2]{ edgeSVector[0][i],edgeSVector[1][i] };
					sum_phifSf[0] += phif * Sf[0] * edgeSign[i];// ���Է��ţ���Ӧ���泯�ڵ�����
					sum_phifSf[1] += phif * Sf[1] * edgeSign[i];
				}
				else {
					// �߽��(���ڱ߳���)
					real phiC = elementField_host.U[jVar][currentElementID];
					real phif = phiC;
					real Sf[2]{ edgeSVector[0][i],edgeSVector[1][i] };
					sum_phifSf[0] += phif * Sf[0] * edgeSign[i];// ���Է��ţ���Ӧ���泯�ڵ�����
					sum_phifSf[1] += phif * Sf[1] * edgeSign[i];
				}

			}

			real nabla_phiC_new[2]{// �ݶ�
				1 / volumeC * sum_phifSf[0],
				1 / volumeC * sum_phifSf[1]
			};
			// ���µ�Ԫ�ݶ�ֵ
			elementField_host.Ux[jVar][currentElementID] = nabla_phiC_new[0];
			elementField_host.Uy[jVar][currentElementID] = nabla_phiC_new[1];
			// ���Ҫ��������Ҫ�����е�Ԫ������ϣ���˵���ѭ��Ҫ����iElementѭ����
		}
	}
}


