#include "../GPUSolver.h"
#include "CalculateGradient.h"
#include "../math/GPUMathKernel.h"

void GPU::calculateGradient(GPU::ElementDataPack& element_device) {
	// --- ���㵥Ԫ�ݶ� --- 
	// ���룺��Ԫ���ꡢ��ԪU����Ԫ�ھ����ꡢ��Ԫ�ھ�U
	// �������Ԫ�ݶ�
	// ���� Reconstructor.cpp
	int block_size = 512;// �����128 256 512
	int grid_size = (element_device.num_element + block_size - 1) / block_size;


	dim3 block(block_size, 1, 1);
	dim3 grid(grid_size, 1, 1);
	calculateGradientKernel <<<grid, block >>> (element_device);
	cudaDeviceSynchronize();


}
void GPU::updateNeighborGradient(GPU::ElementDataPack& element_host, GPU::ElementAdjacent& adjacent) {
	// TODO: �����ھ��ݶ� 1.�����elementAdjacent���飬�洢�ھ�index��Ϣ
	// --- �����ھ��ݶ� ---
	// ���룺��Ԫ�ݶȡ���Ԫ�ھӹ�ϵ����
	// �������Ԫ�ھ��ݶ�
	const int nE = element_host.num_element;
	const int nN = 3;
#pragma omp parallel for
	for (int iE = 0; iE < nE; iE++) {
		for (int iN = 0; iN < nN; iN++) {
			// ��iE����Ԫ�ĵ�iN���ھӵı��
			int index_neighbor = adjacent.neighbors[iN].index[iE];
			// �������Ч�ھ�
			if (index_neighbor >= 0 && index_neighbor < nE) {
				// �ھӾ��ǵ�index_neighbor����Ԫ��ȡ��ֵ
				// ������iE����Ԫ�ĵ�iN���ھӵ�����
				element_host.neighbors[iN].Ux1[iE] = element_host.self.Ux1[index_neighbor];
				element_host.neighbors[iN].Ux2[iE] = element_host.self.Ux2[index_neighbor];
				element_host.neighbors[iN].Ux3[iE] = element_host.self.Ux3[index_neighbor];
				element_host.neighbors[iN].Ux4[iE] = element_host.self.Ux4[index_neighbor];
				element_host.neighbors[iN].Uy1[iE] = element_host.self.Uy1[index_neighbor];
				element_host.neighbors[iN].Uy2[iE] = element_host.self.Uy2[index_neighbor];
				element_host.neighbors[iN].Uy3[iE] = element_host.self.Uy3[index_neighbor];
				element_host.neighbors[iN].Uy4[iE] = element_host.self.Uy4[index_neighbor];

			}
		}
	}
}

__global__ void GPU::calculateGradientKernel(GPU::ElementDataPack& element_device) {
	// --- ���㵥Ԫ�ݶȺ˺��� ����� --- 
	// ���룺��Ԫ���ꡢ��ԪU����Ԫ�ھ����ꡢ��Ԫ�ھ�U
	// �������Ԫ�ݶ�
	const int bid = blockIdx.x;
	const int tid = threadIdx.x;
	const int id = bid * blockDim.x + tid;
	const int& i = id;
	if (i >= element_device.num_element || i < 0) return;

	// �޳���Ч�ھ�
	int nValidNeighbor = 0;
	const int neighborArraySize = 3;
	GPU::Element* validNeighbors[neighborArraySize]{};
	for (int j = 0; j < neighborArraySize; j++) {
		if (!element_device.neighbors[j].isNull[i]) {
			validNeighbors[nValidNeighbor] = &(element_device.neighbors[j]);
			nValidNeighbor += 1;
		}
	}

	// �����ݶ�dUdX
	GPU::Element* self = &(element_device.self);
	const int nVar = 4;// �غ�����������άΪ4
	const int nX = 2;// �����������άΪ2
	const int nVNmax = 3;// ����ھӸ���
	REAL dUdX[nX][nVar]{};// �ѳ�ʼ��Ϊ0
	if (nValidNeighbor == 3) {
		// ��ʼ��dX��dU
		const int nVN = 3;// ��Ч�ھӸ���
		REAL dX[nVN][nX]{};
		REAL dU[nVN][nVar]{};
		for (int iVN = 0; iVN < nVN; iVN++) {
			dX[iVN][0] = validNeighbors[iVN]->x[i] - self->x[i];
			dX[iVN][1] = validNeighbors[iVN]->y[i] - self->y[i];
			dU[iVN][0] = validNeighbors[iVN]->U1 - self->U1;
			dU[iVN][1] = validNeighbors[iVN]->U2 - self->U2;
			dU[iVN][2] = validNeighbors[iVN]->U3 - self->U3;
			dU[iVN][3] = validNeighbors[iVN]->U4 - self->U4;
		}

		// ��С���˷� x=(A'A)^{-1}A'b��dUdX = inv(dXtrans * dX) * dXtrans * dU
		// ! ���º���Ϊdevice��Ӧ��GPUʵ��
		REAL dXtrans[nX][nVN]{};
		REAL invdXdX[nX][nX]{};
		REAL invdXdX_dX[nX][nVN]{};
		GPU::Matrix::transpose(nVN, nX, (REAL*)dX, (REAL*)dXtrans);
		GPU::Matrix::mul_ixj_jxk(nX, nVN, nX, (REAL*)dXtrans, (REAL*)dX, (REAL*)invdXdX);
		GPU::Matrix::inv_2x2(invdXdX);
		GPU::Matrix::mul_ixj_jxk(nX, nX, nVN, (REAL*)invdXdX, (REAL*)dXtrans, (REAL*)invdXdX_dX);
		GPU::Matrix::mul_ixj_jxk(nX, nVN, nVar, (REAL*)invdXdX_dX, (REAL*)dU, (REAL*)dUdX);
	}
	else if (nValidNeighbor == 2) {
		// ���´���������Ĵ�����ͬ(����nVN��ͬ)���޸�ʱ��һͬ�޸�
		// ��ʼ��dX��dU
		const int nVN = 2;
		REAL dX[nVN][nX]{};
		REAL dU[nVN][nVar]{};
		for (int iVN = 0; iVN < nVN; iVN++) {
			dX[iVN][0] = validNeighbors[iVN]->x[i] - self->x[i];
			dX[iVN][1] = validNeighbors[iVN]->y[i] - self->y[i];
			dU[iVN][0] = validNeighbors[iVN]->U1 - self->U1;
			dU[iVN][1] = validNeighbors[iVN]->U2 - self->U2;
			dU[iVN][2] = validNeighbors[iVN]->U3 - self->U3;
			dU[iVN][3] = validNeighbors[iVN]->U4 - self->U4;
		}

		// ��С���˷� x=(A'A)^{-1}A'b��dUdX = inv(dXtrans * dX) * dXtrans * dU
		// ! ���º���Ϊdevice��Ӧ��GPUʵ��
		REAL dXtrans[nX][nVN]{};
		REAL invdXdX[nX][nX]{};
		REAL invdXdX_dX[nX][nVN]{};
		GPU::Matrix::transpose(nVN, nX, (REAL*)dX, (REAL*)dXtrans);
		GPU::Matrix::mul_ixj_jxk(nX, nVN, nX, (REAL*)dXtrans, (REAL*)dX, (REAL*)invdXdX);
		GPU::Matrix::inv_2x2(invdXdX);
		GPU::Matrix::mul_ixj_jxk(nX, nX, nVN, (REAL*)invdXdX, (REAL*)dXtrans, (REAL*)invdXdX_dX);
		GPU::Matrix::mul_ixj_jxk(nX, nVN, nVar, (REAL*)invdXdX_dX, (REAL*)dU, (REAL*)dUdX);
	}
	else if (nValidNeighbor == 1) {
		// ���´���������Ĵ�����ͬ(����nVN��ͬ)���޸�ʱ��һͬ�޸�
		// ��ʼ��dX��dU
		const int nVN = 1;
		REAL dX[nVN][nX]{};
		REAL dU[nVN][nVar]{};
		for (int iVN = 0; iVN < nVN; iVN++) {
			dX[iVN][0] = validNeighbors[iVN]->x[i] - self->x[i];
			dX[iVN][1] = validNeighbors[iVN]->y[i] - self->y[i];
			dU[iVN][0] = validNeighbors[iVN]->U1 - self->U1;
			dU[iVN][1] = validNeighbors[iVN]->U2 - self->U2;
			dU[iVN][2] = validNeighbors[iVN]->U3 - self->U3;
			dU[iVN][3] = validNeighbors[iVN]->U4 - self->U4;
		}

		// ��С���˷� x=(A'A)^{-1}A'b��dUdX = inv(dXtrans * dX) * dXtrans * dU
		// ! ���º���Ϊdevice��Ӧ��GPUʵ��
		REAL dXtrans[nX][nVN]{};
		REAL invdXdX[nX][nX]{};
		REAL invdXdX_dX[nX][nVN]{};
		GPU::Matrix::transpose(nVN, nX, (REAL*)dX, (REAL*)dXtrans);
		GPU::Matrix::mul_ixj_jxk(nX, nVN, nX, (REAL*)dXtrans, (REAL*)dX, (REAL*)invdXdX);
		GPU::Matrix::inv_2x2(invdXdX);
		GPU::Matrix::mul_ixj_jxk(nX, nX, nVN, (REAL*)invdXdX, (REAL*)dXtrans, (REAL*)invdXdX_dX);
		GPU::Matrix::mul_ixj_jxk(nX, nVN, nVar, (REAL*)invdXdX_dX, (REAL*)dU, (REAL*)dUdX);

	}

	// ��dUdX�����Ԫ
	self->Ux1[i] = dUdX[0][0];
	self->Ux2[i] = dUdX[0][1];
	self->Ux3[i] = dUdX[0][2];
	self->Ux4[i] = dUdX[0][3];
	self->Uy1[i] = dUdX[1][0];
	self->Uy2[i] = dUdX[1][1];
	self->Uy3[i] = dUdX[1][2];
	self->Uy4[i] = dUdX[1][3];


	// Barth������ һ�׾��ȸ�ʽ����Ҫ��
	//pE->restructor_in_updateSlope_Barth(f);
	//Limiter::modifySlope_Barth(pE);

	// ����쳣ֵ Ӧ��������
	// ֻ�����Ƿ��С�1e10��������

}
;

__device__ void GPU::Limiter::Barth() {
	// TODO: Barth limiter
	// ���ܣ�������Ԫ�ݶ�
	// ���룺��ԪU���ھ�U���ڵ�U����Ԫ�ݶ�
	// �������Ԫ�ݶ�
}
