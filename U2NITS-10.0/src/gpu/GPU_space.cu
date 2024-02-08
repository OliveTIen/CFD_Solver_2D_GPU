#include "GPUSolver.h"
#include "GPU_space.h"
#include "GPUMathKernel.h"

void GPU::calculateGradient(int block_size, int grid_size, GPU::ElementDataPack& deviceDataPack) {
	// --- ���㵥Ԫ�ݶ� --- 
	// ���룺��Ԫ���ꡢ��ԪU����Ԫ�ھ����ꡢ��Ԫ�ھ�U
	// �������Ԫ�ݶȡ��ھ��ݶ�
	// �μ� Reconstructor.cpp


	for (int i = 0; i < deviceDataPack.num_element; i++) {
		// �޳���Ч�ھ�
		int nValidNeighbor = 0;
		const int neighborArraySize = 3;
		GPU::Element* validNeighbors[neighborArraySize]{};
		for (int j = 0; j < neighborArraySize; j++) {
			if (!deviceDataPack.neighbors[j].isNull[i]) {
				validNeighbors[nValidNeighbor] = &(deviceDataPack.neighbors[j]);
				nValidNeighbor += 1;
			}
		}

		// �����ݶ�dUdX
		GPU::Element* self = &(deviceDataPack.self);
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
	}
}
__global__ void GPU::calculateGradientKernel(GPU::ElementDataPack& deviceDataPack) {
	// --- ���㵥Ԫ�ݶȺ˺��� --- 
	// ���룺��Ԫ���ꡢ��ԪU����Ԫ�ھ����ꡢ��Ԫ�ھ�U
	// �������Ԫ�ݶȡ��ھ��ݶ�
	const unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;


}
;