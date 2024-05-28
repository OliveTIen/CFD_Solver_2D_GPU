//#include <stdio.h>
#include "GradientGPU.h"
#include "../../math/MathGPU.h"
#include "../../output/LogWriter.h"
#include "../../solvers/GPUSolver2.h"
#include "../../global/GlobalPara.h"
#include "../../gpu/GPUGlobalFunction.h"

__global__ void calculateGradientKernel_old(GPU::ElementSoA element_device, GPU::ElementFieldSoA elementField_device) {
	// --- ���㵥Ԫ�ݶȺ˺��� ����� --- 
	// ���룺��Ԫ���ꡢ��ԪU����Ԫ�ھ����ꡢ��Ԫ�ھ�U
	// �������Ԫ�ݶ�
	const int bid = blockIdx.x;
	const int tid = threadIdx.x;
	const int id = bid * blockDim.x + tid;
	const int& i = id;
	if (i >= element_device.num_element || i < 0) return;

	// ��ȡ��Ч�ھ�
	int nValidNeighbor = 0;
	int validNeighbors[3]{ -1,-1,-1 };// ��Ч�ھӵ�ID
	const int neighborArraySize = 3;// �����ε�Ԫ�����3���ھ�
	for (int j = 0; j < neighborArraySize; j++) {
		if (element_device.neighbors[j][i] != -1) {
			nValidNeighbor += 1;
			validNeighbors[nValidNeighbor] = element_device.neighbors[j][i];
		}
	}

	// �����ݶ�dUdX
	const int nVar = 4;// �غ�����������άΪ4
	const int nX = 2;// �����������άΪ2
	//const int nVNmax = 3;// ����ھӸ���
	myfloat dUdX[nX][nVar]{};// �ѳ�ʼ��Ϊ0
	if (nValidNeighbor == 3) {
		// ��ʼ��dX��dU
		const int nVN = 3;// ��Ч�ھӸ���
		myfloat dX[nVN][nX]{};
		myfloat dU[nVN][nVar]{};
		for (int iVN = 0; iVN < nVN; iVN++) {
			int index = validNeighbors[iVN];// �ھӵ�ID
			dX[iVN][0] = element_device.xy[0][index] - element_device.xy[0][i];
			dX[iVN][1] = element_device.xy[1][index] - element_device.xy[1][i];
			for (int iVar = 0; iVar < nVar; iVar++) {
				dU[iVN][iVar] = elementField_device.U[index][iVar] - elementField_device.U[i][iVar];
			}
		}

		// ��С���˷� x=(A'A)^{-1}A'b��dUdX = inv(dXtrans * dX) * dXtrans * dU
		// ! ���º���Ϊdevice��Ӧ��GPUʵ��
		myfloat dXtrans[nX][nVN]{};
		myfloat invdXdX[nX][nX]{};
		myfloat invdXdX_dX[nX][nVN]{};
		GPU::Math::Matrix::transpose(nVN, nX, (myfloat*)dX, (myfloat*)dXtrans);
		GPU::Math::Matrix::mul_ixj_jxk(nX, nVN, nX, (myfloat*)dXtrans, (myfloat*)dX, (myfloat*)invdXdX);
		GPU::Math::Matrix::inv_2x2(invdXdX);
		GPU::Math::Matrix::mul_ixj_jxk(nX, nX, nVN, (myfloat*)invdXdX, (myfloat*)dXtrans, (myfloat*)invdXdX_dX);
		GPU::Math::Matrix::mul_ixj_jxk(nX, nVN, nVar, (myfloat*)invdXdX_dX, (myfloat*)dU, (myfloat*)dUdX);
	}
	else if (nValidNeighbor == 2) {
		// ���´���������Ĵ�����ͬ(����nVN��ͬ)���޸�ʱ��һͬ�޸�
		// ��ʼ��dX��dU
		const int nVN = 2;// ��Ч�ھӸ���
		myfloat dX[nVN][nX]{};
		myfloat dU[nVN][nVar]{};
		for (int iVN = 0; iVN < nVN; iVN++) {
			int index = validNeighbors[iVN];// �ھӵ�ID
			dX[iVN][0] = element_device.xy[0][index] - element_device.xy[0][i];
			dX[iVN][1] = element_device.xy[1][index] - element_device.xy[1][i];
			for (int iVar = 0; iVar < nVar; iVar++) {
				dU[iVN][iVar] = elementField_device.U[index][iVar] - elementField_device.U[i][iVar];
			}
		}

		// ��С���˷� x=(A'A)^{-1}A'b��dUdX = inv(dXtrans * dX) * dXtrans * dU
		// ! ���º���Ϊdevice��Ӧ��GPUʵ��
		myfloat dXtrans[nX][nVN]{};
		myfloat invdXdX[nX][nX]{};
		myfloat invdXdX_dX[nX][nVN]{};
		GPU::Math::Matrix::transpose(nVN, nX, (myfloat*)dX, (myfloat*)dXtrans);
		GPU::Math::Matrix::mul_ixj_jxk(nX, nVN, nX, (myfloat*)dXtrans, (myfloat*)dX, (myfloat*)invdXdX);
		GPU::Math::Matrix::inv_2x2(invdXdX);
		GPU::Math::Matrix::mul_ixj_jxk(nX, nX, nVN, (myfloat*)invdXdX, (myfloat*)dXtrans, (myfloat*)invdXdX_dX);
		GPU::Math::Matrix::mul_ixj_jxk(nX, nVN, nVar, (myfloat*)invdXdX_dX, (myfloat*)dU, (myfloat*)dUdX);
	}
	else if (nValidNeighbor == 1) {
		// ���´���������Ĵ�����ͬ(����nVN��ͬ)���޸�ʱ��һͬ�޸�
		// ��ʼ��dX��dU
		const int nVN = 1;// ��Ч�ھӸ���
		myfloat dX[nVN][nX]{};
		myfloat dU[nVN][nVar]{};
		for (int iVN = 0; iVN < nVN; iVN++) {
			int index = validNeighbors[iVN];// �ھӵ�ID
			dX[iVN][0] = element_device.xy[0][index] - element_device.xy[0][i];
			dX[iVN][1] = element_device.xy[1][index] - element_device.xy[1][i];
			for (int iVar = 0; iVar < nVar; iVar++) {
				dU[iVN][iVar] = elementField_device.U[index][iVar] - elementField_device.U[i][iVar];
			}
		}

		// ��С���˷� x=(A'A)^{-1}A'b��dUdX = inv(dXtrans * dX) * dXtrans * dU
		// ! ���º���Ϊdevice��Ӧ��GPUʵ��
		myfloat dXtrans[nX][nVN]{};
		myfloat invdXdX[nX][nX]{};
		myfloat invdXdX_dX[nX][nVN]{};
		GPU::Math::Matrix::transpose(nVN, nX, (myfloat*)dX, (myfloat*)dXtrans);
		GPU::Math::Matrix::mul_ixj_jxk(nX, nVN, nX, (myfloat*)dXtrans, (myfloat*)dX, (myfloat*)invdXdX);
		GPU::Math::Matrix::inv_2x2(invdXdX);
		GPU::Math::Matrix::mul_ixj_jxk(nX, nX, nVN, (myfloat*)invdXdX, (myfloat*)dXtrans, (myfloat*)invdXdX_dX);
		GPU::Math::Matrix::mul_ixj_jxk(nX, nVN, nVar, (myfloat*)invdXdX_dX, (myfloat*)dU, (myfloat*)dUdX);
	}

	// ��dUdX�����Ԫ
	for (int j = 0; j < nVar; j++) {
		elementField_device.Ux[j][i] = dUdX[0][j];
		elementField_device.Uy[j][i] = dUdX[1][j];
	}

	// Barth������ һ�׾��ȸ�ʽ����Ҫ��
	//pE->restructor_in_updateSlope_Barth(f);
	//Limiter::modifySlope_Barth(pE);

	// ����쳣ֵ Ӧ��������
	// ֻ�����Ƿ��С�1e10��������
}

__global__ void GradientLeastSquareKernel(GPU::ElementSoA element_device, GPU::ElementFieldSoA elementField_device, GPU::NodeSoA node_device) {
	using namespace GPU;
	const int bid = blockIdx.x;
	const int tid = threadIdx.x;
	const int iElement = bid * blockDim.x + tid;

	if (iElement >= element_device.num_element || iElement < 0) return;
	
	// ��ȡ��Ч�ھ�
	int numOfValidNeighbor = 0;
	int validNeighbors[3]{ -1,-1,-1 };// ��Ч�ھӵ�ID
	const int neighborArraySize = 3;// �����ε�Ԫ�����3���ھ�
	for (int j = 0; j < neighborArraySize; j++) {
		if (element_device.neighbors[j][iElement] != -1) {
			validNeighbors[numOfValidNeighbor] = element_device.neighbors[j][iElement];
			numOfValidNeighbor += 1;
		}
	}
	if (numOfValidNeighbor <= 0) {
		//LogWriter::logAndPrint("Error: invalid numOfValidNeighbor.\n", LogWriter::Error, LogWriter::Error);
		printf("Error: invalid numOfValidNeighbor.\n");
		//exit(-1);
	}
	// ���������
	const int nVar = 4;// �غ�����������άΪ4
	const int nX = 2;// �����������άΪ2
	const int& nVN = numOfValidNeighbor;// ��Ч�ھӸ�����ȡֵ1-3
	myfloat* dX = new myfloat[nVN * nX]{};
	myfloat* dUdX = new myfloat[nX * nVar]{};
	myfloat* dU = new myfloat[nVN * nVar]{};
	myfloat* dXtrans = new myfloat[nX * nVN]{};
	myfloat* invdXdX = new myfloat[nX * nX]{};
	myfloat* invdXdX_dX = new myfloat[nX * nVN]{};
	// ��ʼ��dX��dU
	for (int iVN = 0; iVN < nVN; iVN++) {
		int index = validNeighbors[iVN];// �ھӵ�ID
		dX[iVN * nX + 0] = element_device.xy[0][index] - element_device.xy[0][iElement];
		dX[iVN * nX + 1] = element_device.xy[1][index] - element_device.xy[1][iElement];
		for (int iVar = 0; iVar < nVar; iVar++) {
			dU[iVN * nVar + iVar] = elementField_device.U[iVar][index] - elementField_device.U[iVar][iElement];
		}
	}
	// ��С���˷�����dUdX��x=(A'A)^{-1}A'b��dUdX = inv(dXtrans * dX) * dXtrans * dU
	GPU::Math::Matrix::transpose(nVN, nX, (myfloat*)dX, (myfloat*)dXtrans);
	GPU::Math::Matrix::mul_ixj_jxk(nX, nVN, nX, (myfloat*)dXtrans, (myfloat*)dX, (myfloat*)invdXdX);
	GPU::Math::Matrix::inv_2x2(invdXdX);
	GPU::Math::Matrix::mul_ixj_jxk(nX, nX, nVN, (myfloat*)invdXdX, (myfloat*)dXtrans, (myfloat*)invdXdX_dX);
	GPU::Math::Matrix::mul_ixj_jxk(nX, nVN, nVar, (myfloat*)invdXdX_dX, (myfloat*)dU, (myfloat*)dUdX);

	/*
������
һά�£�minmod��vanleer��������Ŀ���Ǳ�֤�����������������ֵ���������ڵ�Ԫ����������ֵ
�μ�����������ѧ�����������μ��ϼ�P106

���Ƶأ����Ա�֤�����ζ��㴦U���������ڵ�ԪU
���ݵ�ǰб��Ux, Uy�����㶥�����ֵdU_node��ǰ���Ѿ�����ھ��������ֵdU
����ratio = max(abs(dU_node/dU))��������1�����������ratio
		*/
	const int nNode = 3;
	myfloat dU_node[nNode][nVar]{};
	myfloat dX_node[nNode][nX]{};
	// ����dX_node
	for (int iNode = 0; iNode < nNode; iNode++) {
		for (int iVar = 0; iVar < nVar; iVar++) {
			int nodeID = element_device.nodes[iNode][iElement];
			dX_node[iNode][0] = node_device.xy[0][nodeID] - element_device.xy[0][iElement];
			dX_node[iNode][1] = node_device.xy[1][nodeID] - element_device.xy[1][iElement];
		}
	}
	// ����dU_node[nNode,nVar] =  dX_node[nNode,nX] * dUdX[nX,nVar]��Ȼ�����ratio
	GPU::Math::Matrix::mul_ixj_jxk(nNode, nX, nVar, (myfloat*)dX_node, dUdX, (myfloat*)dU_node);
	myfloat ratio = 1.0;
	for (int iNode = 0; iNode < nNode; iNode++) {
		for (int iVar = 0; iVar < nVar; iVar++) {
			using namespace GPU::Math;
			if (dU[iNode * nVar + iVar] == 0)continue;// ����
			ratio = GPU::Math::max(ratio, Math::abs(dU_node[iNode][iVar] / dU[iNode * nVar + iVar]));
		}
	}
	// ratioһ�����ڵ���1.0����˿���ֱ�ӳ���ratio
	GPU::Math::Matrix::div_matrix_by_scalar(nX, nVar, dUdX, ratio);
	// ��dUdX�����Ԫ
	for (int jVar = 0; jVar < nVar; jVar++) {
		elementField_device.Ux[jVar][iElement] = dUdX[0 * nVar + jVar];
		elementField_device.Uy[jVar][iElement] = dUdX[1 * nVar + jVar];
	}

	delete[] dX;
	delete[] dUdX;
	delete[] dU;
	delete[] dXtrans;
	delete[] invdXdX;
	delete[] invdXdX_dX;

}

__global__ void GradientGreenGaussKernel(GPU::ElementSoA element, GPU::ElementFieldSoA elementField, GPU::EdgeSoA edge) {
	/*
	���ܴ�element���ã���Ϊ���ñ�����ָ�룬��ָ��洢����host��ַ����device���Ҳ�����Ҳ���ܴ�ֵ����Ϊ�����classĬ�Ϲ��캯�������캯����host�ģ�device�޷�����
	��Ϊstruct������structû��Ĭ�Ϲ�����������������Դ�ֵ
	*/
	using namespace GPU;

	const int iElement = blockIdx.x * blockDim.x + threadIdx.x;
	if (iElement >= element.num_element || iElement < 0) return;

	int num_edge = edge.num_edge;
	// ��ǰ��ԪID����ǰ��Ԫ�����edgeID��edge�����������ⵥԪID
	int currentElementID = element.ID[iElement];
	myfloat volumeC = element.volume[currentElementID];// ���
	// ��ʼ��edgeID �������鷶Χ�ķǷ�ֵ��Ϊ-1
	int edgeID[4]{ -1,-1,-1,-1 };// ��iElement��element�ĵ�i���ߵ�ID
	for (int i = 0; i < 4; i++) {
		edgeID[i] = element.edges[i][iElement];
		if (edgeID[i] <= -1 || edgeID[i] >= num_edge) {
			edgeID[i] = -1;
		}
	}
	// ��ʼ��edgeType��edgeType��setID�������߽��ID��-1��ʾ�ڲ���
	int edgeType[4]{ -1,-1,-1,-1 };
	for (int i = 0; i < 4; i++) {
		int edgeIDi = edgeID[i];
		if (edgeIDi == -1) {
			continue;
		}

		edgeType[i] = edge.setID[edgeIDi];
	}
	// ��ʼ��edgeSign��outElementID
	int edgeSign[4]{ 0,0,0,0 };// 1��ʾ�߳��⣬-1��ʾ�߳���
	int outElementID[4]{ -1,-1,-1,-1 };//�ⵥԪID �ڲ��ߺ����ڱ�
	for (int i = 0; i < 4; i++) {
		int edgeIDi = edgeID[i];
		if (edgeIDi == -1) {
			continue;
		}

		if (currentElementID != edge.elementR[edgeIDi]) {
			edgeSign[i] = 1;// currentElement=elementL������
			outElementID[i] = edge.elementR[edgeIDi];
		}
		else {
			edgeSign[i] = -1;// currentElement=elementR������
			outElementID[i] = edge.elementL[edgeIDi];
		}
	}
	// �����ʸ��edgeSVector(������)
	myfloat edgeNormal[2][4]{};
	myfloat edgeArea[4]{};
	myfloat edgeSVector[2][4]{};// ���ʸ��
	for (int i = 0; i < 4; i++) {
		int edgeIDi = edgeID[i];
		if (edgeIDi == -1) {
			continue;
		}

		edgeNormal[0][i] = edge.normal[0][edgeIDi];
		edgeNormal[1][i] = edge.normal[1][edgeIDi];
		edgeArea[i] = edge.length[edgeIDi];
		edgeSVector[0][i] = edgeNormal[0][i] * edgeArea[i];
		edgeSVector[1][i] = edgeNormal[1][i] * edgeArea[i];
	}

	// �󼸺�Ȩ������gc��ʸ��(������)
	myfloat gC_all[4]{};
	myfloat rf_rC_all[4][2]{};
	myfloat rf_rF_all[4][2]{};
	for (int i = 0; i < 4; i++) {
		int edgeIDi = edgeID[i];
		if (edgeIDi == -1) {
			continue;
		}

		int currentOutElementID = outElementID[i];
		if (currentOutElementID != -1) {
			// �ڲ��ߺ����ڱ�
			if (edgeType[i] == -1) {
				// �ڲ���

				// ��͵�Ԫ����
				myfloat fx = edge.xy[0][edgeIDi];// face x
				myfloat fy = edge.xy[1][edgeIDi];
				myfloat Cx = element.xy[0][currentElementID];// cell x
				myfloat Cy = element.xy[1][currentElementID];
				myfloat Fx = element.xy[0][currentOutElementID];
				myfloat Fy = element.xy[1][currentOutElementID];
				// ����Ȩ������gC = |rF-rf|/|rF-rC|���洢��gC_all������
				myfloat dFf = Math::distance2D(Fx, Fy, fx, fy);
				myfloat dFC = Math::distance2D(Fx, Fy, Cx, Cy);
				if (dFC < U2NITS::Math::EPSILON)dFC = U2NITS::Math::EPSILON;
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
				myfloat fx = edge.xy[0][edgeIDi];
				myfloat fy = edge.xy[1][edgeIDi];
				myfloat Cx = element.xy[0][currentElementID];// cell x
				myfloat Cy = element.xy[1][currentElementID];
				// ���ڱߵ�outElement������Ҫƽ��
				int pairIDi = edge.periodicPair[edgeIDi];
				myfloat fx_pair = edge.xy[0][pairIDi];
				myfloat fy_pair = edge.xy[1][pairIDi];
				myfloat shiftx = fx - fx_pair;// ��pairָ��ǰλ�õ�ʸ��
				myfloat shifty = fy - fy_pair;
				myfloat Fx = element.xy[0][currentOutElementID];
				myfloat Fy = element.xy[1][currentOutElementID];
				Fx += shiftx;// �ø�ʸ��ƽ�Ƶ�Ԫ
				Fy += shifty;

				// ������ڲ�����һ����
				// ����Ȩ������gC = |rF-rf|/|rF-rC|���洢��gC_all������
				myfloat dFf = Math::distance2D(Fx, Fy, fx, fy);
				myfloat dFC = Math::distance2D(Fx, Fy, Cx, Cy);
				if (dFC < U2NITS::Math::EPSILON)dFC = U2NITS::Math::EPSILON;
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
			myfloat fx = edge.xy[0][edgeIDi];// face x
			myfloat fy = edge.xy[1][edgeIDi];
			myfloat Cx = element.xy[0][currentElementID];// cell x
			myfloat Cy = element.xy[1][currentElementID];

			using namespace U2NITS;
			myfloat gC = 0.5;// ����Ȩ������gc
			myfloat rf_rC[2]{ fx - Cx,fy - Cy };// rf-rC
			myfloat rf_rF[2]{ -(fx - Cx),-(fy - Cy) };// rf-rF

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
		myfloat sum_phifSf[2]{};
		for (int i = 0; i < 4; i++) {
			int currentEdgeID = edgeID[i];
			if (currentEdgeID == -1) {
				continue;
			}
			int currentOutElementID = outElementID[i];

			if (currentOutElementID != -1) {
				// �ڲ��ߺ����ڱ�

				// ��ȡ��Ԫֵ�͵�Ԫ�ݶ�ֵ
				myfloat phiC = elementField.U[jVar][currentElementID];
				myfloat phiF = elementField.U[jVar][currentOutElementID];
				// �����ֵ(����)
				myfloat gC = gC_all[i];
				myfloat gC1 = 1.0 - gC;
				myfloat phif = gC * phiC + gC1 * phiF;
				myfloat Sf[2]{ edgeSVector[0][i],edgeSVector[1][i] };
				sum_phifSf[0] += phif * Sf[0] * edgeSign[i];// ���Է��ţ���Ӧ���泯�ڵ�����
				sum_phifSf[1] += phif * Sf[1] * edgeSign[i];
			}
			else {
				// �߽��(���ڱ߳���)
				myfloat phiC = elementField.U[jVar][currentElementID];
				myfloat phif = phiC;
				myfloat Sf[2]{ edgeSVector[0][i],edgeSVector[1][i] };
				sum_phifSf[0] += phif * Sf[0] * edgeSign[i];// ���Է��ţ���Ӧ���泯�ڵ�����
				sum_phifSf[1] += phif * Sf[1] * edgeSign[i];
			}

		}

		myfloat nabla_phiC_new[2]{// �ݶ�
			1 / volumeC * sum_phifSf[0],
			1 / volumeC * sum_phifSf[1]
		};
		// ���µ�Ԫ�ݶ�ֵ
		elementField.Ux[jVar][currentElementID] = nabla_phiC_new[0];
		elementField.Uy[jVar][currentElementID] = nabla_phiC_new[1];
		// ���Ҫ��������Ҫ�����е�Ԫ������ϣ���˵���ѭ��Ҫ����iElementѭ����
	}
}

void GPU::calculateGradient_old(GPU::ElementSoA& element_device, GPU::ElementFieldSoA& elementField_device) {
	int block_size = 512;// �����128 256 512
	// num element = 10216
	int grid_size = (element_device.num_element + block_size - 1) / block_size;
	dim3 block(block_size, 1, 1);
	dim3 grid(grid_size, 1, 1);
	calculateGradientKernel_old <<<grid, block >>> (element_device, elementField_device);

	catchCudaErrorAndExit();
	// cudaDeviceSynchronize();��������
}

__device__ inline void getElementT3ValidNeighbor(int validNeighbors[3], int& numOfValidNeighbor, myint iElement, GPU::ElementSoA& element_device) {
	numOfValidNeighbor = 0;
	for (int j = 0; j < 3; j++) {// �����ε�Ԫ�����3���ھ�
		if (element_device.neighbors[j][iElement] != -1) {
			validNeighbors[numOfValidNeighbor] = element_device.neighbors[j][iElement];
			numOfValidNeighbor += 1;
		}
	}
}

__device__ inline void get_U_linear(myfloat x, myfloat y, myfloat& U_dist, myint i_e, const myfloat* x_e, const myfloat* y_e, const myfloat* U_e, const myfloat* Ux_e, const myfloat* Uy_e) {
	U_dist = U_e[i_e] + Ux_e[i_e] * (x - x_e[i_e]) + Uy_e[i_e] * (y - y_e[i_e]);
}

__global__ void LimiterBarth_device_kernel(GPU::ElementFieldSoA elementField_device, GPU::ElementSoA element_device, GPU::NodeSoA node_device) {
	const int iElement = blockIdx.x * blockDim.x + threadIdx.x;
	if (iElement >= element_device.num_element || iElement < 0) return;

	int validNeighbors[3]{ -1,-1,-1 };
	int numOfValidNeighbor = 0;
	getElementT3ValidNeighbor(validNeighbors, numOfValidNeighbor, iElement, element_device);
	const auto& U_element = elementField_device.U;
	auto& Ux_element = elementField_device.Ux;
	auto& Uy_element = elementField_device.Uy;

	for (int iVar = 0; iVar < 4; iVar++) {
		// ��ȡ�����ھӵ�Ԫ����ֵ���½�U_lower���Ͻ�U_upper
		myfloat U_element_center = U_element[iVar][iElement];
		myfloat U_lower = U_element_center;
		myfloat U_upper = U_element_center;
		for (int j = 0; j < numOfValidNeighbor; j++) {
			myint iNeighbor = validNeighbors[j];
			U_lower = GPU::Math::min(U_lower, U_element[iVar][iNeighbor]);
			U_upper = GPU::Math::max(U_upper, U_element[iVar][iNeighbor]);
		}

		// ����������ϵ��phi��ȡ����������ϵ��phi_node����Сֵ
		myfloat phi = 1.0;
		for (int j = 0; j < 3; j++) {
			myint iNode = element_device.nodes[j][iElement];
			myfloat x_node = node_device.xy[0][iNode];
			myfloat y_node = node_device.xy[1][iNode];
			myfloat U_node = 0.0f;
			get_U_linear(x_node, y_node, U_node, iElement, element_device.xy[0], element_device.xy[1], U_element[iVar], Ux_element[iVar], Uy_element[iVar]);
			myfloat U_node_relative = U_node - U_element_center;
			myfloat U_upper_relative = U_upper - U_element_center;
			myfloat U_lower_relative = U_lower - U_element_center;
			myfloat phi_node = 1.0;
			if (U_node_relative > 0) {
				U_node_relative += U2NITS::Math::EPSILON;
				phi_node = GPU::Math::min(1, U_upper_relative / U_node_relative);
			}
			else if (U_node_relative < 0) {
				U_node_relative -= U2NITS::Math::EPSILON;
				phi_node = GPU::Math::min(1, U_lower_relative / U_node_relative);
			}
			phi = GPU::Math::min(phi, phi_node);
		}

		// �����ݶ�
		Ux_element[iVar][iElement] *= phi;
		Uy_element[iVar][iElement] *= phi;
	}
}

void LimiterBarth_device(GPU::ElementFieldSoA& elementField_device, GPU::ElementSoA& element_device, GPU::NodeSoA& node_device) {
	int block_size = GPU::MY_BLOCK_SIZE;
	int grid_size = (element_device.num_element + block_size - 1) / block_size;
	dim3 block(block_size, 1, 1);
	dim3 grid(grid_size, 1, 1);
	LimiterBarth_device_kernel <<<grid, block >>> (elementField_device, element_device, node_device);
	getLastCudaError("LimiterBarth_device failed.");
}

void GPU::Space::Gradient::Gradient_2(GPU::ElementSoA& element_device, GPU::NodeSoA& node_device, GPU::EdgeSoA& edge_device, GPU::ElementFieldSoA& elementField_device) {
	int block_size = GPU::MY_BLOCK_SIZE;
	int grid_size = (element_device.num_element + block_size - 1) / block_size;
	static bool is_called_first_time = true;
	// ���ݶ�
	switch (GlobalPara::inviscid_flux_method::flag_gradient) {
	case _GRA_leastSquare:
		LogWriter::logAndPrintError("Unimplemented gradient type.\n");
		exit(-1);
		break;
	case _GRA_greenGauss:
		GradientGreenGauss(block_size, grid_size, element_device, elementField_device, edge_device);// �Ѽ��
		break;

	default:
		LogWriter::logAndPrintError("invalid graident type.\n");
		exit(-1);
		break;
	}
	// ������
	switch (GlobalPara::inviscid_flux_method::flux_limiter) {
	case _LIM_none:
		if (is_called_first_time) {
			LogWriter::logAndPrintWarning("No limiter used.\n");
		}
		break;

	case _LIM_barth:
		LimiterBarth_device(elementField_device, element_device, node_device);// �Ѽ��
		break;

	default:
		LogWriter::logAndPrintError("invalid limiter type.\n");
		exit(-1);
		break;
	}

	getLastCudaError("GPU::Space::Gradient::Gradient_2 failed.");
}

void GPU::Space::Gradient::GradientLeastSquare(int block_size, int grid_size, GPU::ElementSoA& element_device, GPU::ElementFieldSoA& elementField_device, GPU::NodeSoA& node_device) {
	dim3 block(block_size, 1, 1);
	dim3 grid(grid_size, 1, 1);
	GradientLeastSquareKernel <<<grid, block >>> (element_device, elementField_device, node_device);
	getLastCudaError("GradientLeastSquare failed.");
}

void GPU::Space::Gradient::GradientGreenGauss(int block_size, int grid_size, GPU::ElementSoA& element_device, GPU::ElementFieldSoA& elementField_device, GPU::EdgeSoA& edge_device) {
	dim3 block(block_size, 1, 1);
	dim3 grid(grid_size, 1, 1);
	GradientGreenGaussKernel <<<grid, block>>> (element_device, elementField_device, edge_device);
	getLastCudaError("GradientGreenGauss failed.");
}
