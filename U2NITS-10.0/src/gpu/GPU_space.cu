#include "GPUSolver.h"
#include "GPU_space.h"
#include "GPUMathKernel.h"

void GPU::calculateGradient(int block_size, int grid_size, GPU::ElementDataPack& deviceDataPack) {
	// --- 计算单元梯度 --- 
	// 输入：单元坐标、单元U、单元邻居坐标、单元邻居U
	// 输出：单元梯度、邻居梯度
	// 参加 Reconstructor.cpp


	for (int i = 0; i < deviceDataPack.num_element; i++) {
		// 剔除无效邻居
		int nValidNeighbor = 0;
		const int neighborArraySize = 3;
		GPU::Element* validNeighbors[neighborArraySize]{};
		for (int j = 0; j < neighborArraySize; j++) {
			if (!deviceDataPack.neighbors[j].isNull[i]) {
				validNeighbors[nValidNeighbor] = &(deviceDataPack.neighbors[j]);
				nValidNeighbor += 1;
			}
		}

		// 计算梯度dUdX
		GPU::Element* self = &(deviceDataPack.self);
		const int nVar = 4;// 守恒量个数，二维为4
		const int nX = 2;// 坐标个数，二维为2
		const int nVNmax = 3;// 最大邻居个数
		REAL dUdX[nX][nVar]{};// 已初始化为0
		if (nValidNeighbor == 3) {
			// 初始化dX、dU
			const int nVN = 3;// 有效邻居个数
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

			// 最小二乘法 x=(A'A)^{-1}A'b，dUdX = inv(dXtrans * dX) * dXtrans * dU
			// ! 以下函数为device，应在GPU实现
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
			// 以下代码与上面的代码相同(仅仅nVN不同)，修改时请一同修改
			// 初始化dX、dU
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

			// 最小二乘法 x=(A'A)^{-1}A'b，dUdX = inv(dXtrans * dX) * dXtrans * dU
			// ! 以下函数为device，应在GPU实现
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
			// 以下代码与上面的代码相同(仅仅nVN不同)，修改时请一同修改
			// 初始化dX、dU
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

			// 最小二乘法 x=(A'A)^{-1}A'b，dUdX = inv(dXtrans * dX) * dXtrans * dU
			// ! 以下函数为device，应在GPU实现
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
	// --- 计算单元梯度核函数 --- 
	// 输入：单元坐标、单元U、单元邻居坐标、单元邻居U
	// 输出：单元梯度、邻居梯度
	const unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;


}
;