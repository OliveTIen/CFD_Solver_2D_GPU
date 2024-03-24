//#include <stdio.h>
#include "GradientGPU.h"
#include "../../math/MathGPU.h"
#include "../../output/LogWriter.h"
#include "../../solvers/GPUSolver2.h"

void GPU::calculateGradient_old(GPU::ElementSoA& element_device, GPU::FieldSoA& elementField_device) {
	int block_size = 512;// �����128 256 512
	// num element = 10216
	int grid_size = (element_device.num_element + block_size - 1) / block_size;
	dim3 block(block_size, 1, 1);
	dim3 grid(grid_size, 1, 1);
	calculateGradientKernel_old <<<grid, block >>> (element_device, elementField_device);
	cudaError_t cuda_error = cudaGetLastError();
	if (cuda_error != 0) {
		std::string e = "cudaError=" + std::to_string(cuda_error) + ", " + cudaGetErrorString(cuda_error);
		LogWriter::writeLogAndCout(e, LogWriter::Error, LogWriter::Error);
	}
	// cudaDeviceSynchronize();��������
}

__global__ void GPU::calculateGradientKernel_old(GPU::ElementSoA& element_device, GPU::FieldSoA& elementField_device) {
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
	REAL dUdX[nX][nVar]{};// �ѳ�ʼ��Ϊ0
	if (nValidNeighbor == 3) {
		// ��ʼ��dX��dU
		const int nVN = 3;// ��Ч�ھӸ���
		REAL dX[nVN][nX]{};
		REAL dU[nVN][nVar]{};
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
		REAL dXtrans[nX][nVN]{};
		REAL invdXdX[nX][nX]{};
		REAL invdXdX_dX[nX][nVN]{};
		GPU::Math::Matrix::transpose(nVN, nX, (REAL*)dX, (REAL*)dXtrans);
		GPU::Math::Matrix::mul_ixj_jxk(nX, nVN, nX, (REAL*)dXtrans, (REAL*)dX, (REAL*)invdXdX);
		GPU::Math::Matrix::inv_2x2(invdXdX);
		GPU::Math::Matrix::mul_ixj_jxk(nX, nX, nVN, (REAL*)invdXdX, (REAL*)dXtrans, (REAL*)invdXdX_dX);
		GPU::Math::Matrix::mul_ixj_jxk(nX, nVN, nVar, (REAL*)invdXdX_dX, (REAL*)dU, (REAL*)dUdX);
	}
	else if (nValidNeighbor == 2) {
		// ���´���������Ĵ�����ͬ(����nVN��ͬ)���޸�ʱ��һͬ�޸�
		// ��ʼ��dX��dU
		const int nVN = 2;// ��Ч�ھӸ���
		REAL dX[nVN][nX]{};
		REAL dU[nVN][nVar]{};
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
		REAL dXtrans[nX][nVN]{};
		REAL invdXdX[nX][nX]{};
		REAL invdXdX_dX[nX][nVN]{};
		GPU::Math::Matrix::transpose(nVN, nX, (REAL*)dX, (REAL*)dXtrans);
		GPU::Math::Matrix::mul_ixj_jxk(nX, nVN, nX, (REAL*)dXtrans, (REAL*)dX, (REAL*)invdXdX);
		GPU::Math::Matrix::inv_2x2(invdXdX);
		GPU::Math::Matrix::mul_ixj_jxk(nX, nX, nVN, (REAL*)invdXdX, (REAL*)dXtrans, (REAL*)invdXdX_dX);
		GPU::Math::Matrix::mul_ixj_jxk(nX, nVN, nVar, (REAL*)invdXdX_dX, (REAL*)dU, (REAL*)dUdX);
	}
	else if (nValidNeighbor == 1) {
		// ���´���������Ĵ�����ͬ(����nVN��ͬ)���޸�ʱ��һͬ�޸�
		// ��ʼ��dX��dU
		const int nVN = 1;// ��Ч�ھӸ���
		REAL dX[nVN][nX]{};
		REAL dU[nVN][nVar]{};
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
		REAL dXtrans[nX][nVN]{};
		REAL invdXdX[nX][nX]{};
		REAL invdXdX_dX[nX][nVN]{};
		GPU::Math::Matrix::transpose(nVN, nX, (REAL*)dX, (REAL*)dXtrans);
		GPU::Math::Matrix::mul_ixj_jxk(nX, nVN, nX, (REAL*)dXtrans, (REAL*)dX, (REAL*)invdXdX);
		GPU::Math::Matrix::inv_2x2(invdXdX);
		GPU::Math::Matrix::mul_ixj_jxk(nX, nX, nVN, (REAL*)invdXdX, (REAL*)dXtrans, (REAL*)invdXdX_dX);
		GPU::Math::Matrix::mul_ixj_jxk(nX, nVN, nVar, (REAL*)invdXdX_dX, (REAL*)dU, (REAL*)dUdX);
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

void GPU::Space::Gradient::GradientLeastSquare(GPU::ElementSoA& element_device, GPU::FieldSoA& elementField_device, GPU::NodeSoA& node_device) {
	int block_size = 512;// �����128 256 512
	int grid_size = (element_device.num_element + block_size - 1) / block_size;
	dim3 block(block_size, 1, 1);
	dim3 grid(grid_size, 1, 1);
	GradientLeastSquareKernel <<<grid, block >>> (element_device, elementField_device, node_device);
	cudaError_t cuda_error = cudaGetLastError();
	if (cuda_error != 0) {
		std::string e = "cudaError=" + std::to_string(cuda_error) + ", " + cudaGetErrorString(cuda_error);
		LogWriter::writeLogAndCout(e, LogWriter::Error, LogWriter::Error);
		exit(cuda_error);
	}

}

__global__  void GPU::Space::Gradient::GradientLeastSquareKernel(GPU::ElementSoA& element_host, GPU::FieldSoA& elementField_host, GPU::NodeSoA& node_host) {
	const int bid = blockIdx.x;
	const int tid = threadIdx.x;
	const int iElement = bid * blockDim.x + tid;

	if (iElement >= element_host.num_element || iElement < 0) return;
	
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
		//LogWriter::writeLogAndCout("Error: invalid numOfValidNeighbor.\n", LogWriter::Error, LogWriter::Error);
		printf("Error: invalid numOfValidNeighbor.\n");
		//exit(-1);
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
	GPU::Math::Matrix::transpose(nVN, nX, (real*)dX, (real*)dXtrans);
	GPU::Math::Matrix::mul_ixj_jxk(nX, nVN, nX, (real*)dXtrans, (real*)dX, (real*)invdXdX);
	GPU::Math::Matrix::inv_2x2(invdXdX);
	GPU::Math::Matrix::mul_ixj_jxk(nX, nX, nVN, (real*)invdXdX, (real*)dXtrans, (real*)invdXdX_dX);
	GPU::Math::Matrix::mul_ixj_jxk(nX, nVN, nVar, (real*)invdXdX_dX, (real*)dU, (real*)dUdX);

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
	GPU::Math::Matrix::mul_ixj_jxk(nNode, nX, nVar, (real*)dX_node, dUdX, (real*)dU_node);
	real ratio = 1.0;
	for (int iNode = 0; iNode < nNode; iNode++) {
		for (int iVar = 0; iVar < nVar; iVar++) {
			using namespace GPU::Math;
			if (dU[iNode * nVar + iVar] == 0)continue;// ����
			ratio = max(ratio, abs(dU_node[iNode][iVar] / dU[iNode * nVar + iVar]));
		}
	}
	// ratioһ�����ڵ���1.0����˿���ֱ�ӳ���ratio
	GPU::Math::Matrix::div_matrix_by_scalar(nX, nVar, dUdX, ratio);
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
