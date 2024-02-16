#include "../GPUSolver.h"
#include "CalculateGradient.h"
#include "../math/GPUMathKernel.h"

void GPU::calculateGradient(GPU::ElementDataPack& element_device) {
	// --- 计算单元梯度 --- 
	// 输入：单元坐标、单元U、单元邻居坐标、单元邻居U
	// 输出：单元梯度
	// 参照 Reconstructor.cpp
	int block_size = 512;// 最好是128 256 512
	int grid_size = (element_device.num_element + block_size - 1) / block_size;


	dim3 block(block_size, 1, 1);
	dim3 grid(grid_size, 1, 1);
	calculateGradientKernel <<<grid, block >>> (element_device);
	cudaDeviceSynchronize();


}
void GPU::updateNeighborGradient(GPU::ElementDataPack& element_host, GPU::ElementAdjacent& adjacent) {
	// TODO: 更新邻居梯度 1.已添加elementAdjacent数组，存储邻居index信息
	// --- 更新邻居梯度 ---
	// 输入：单元梯度、单元邻居关系矩阵
	// 输出：单元邻居梯度
	const int nE = element_host.num_element;
	const int nN = 3;
#pragma omp parallel for
	for (int iE = 0; iE < nE; iE++) {
		for (int iN = 0; iN < nN; iN++) {
			// 第iE个单元的第iN个邻居的编号
			int index_neighbor = adjacent.neighbors[iN].index[iE];
			// 如果是有效邻居
			if (index_neighbor >= 0 && index_neighbor < nE) {
				// 邻居就是第index_neighbor个单元，取其值
				// 赋给第iE个单元的第iN个邻居的数据
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
	// --- 计算单元梯度核函数 已完成 --- 
	// 输入：单元坐标、单元U、单元邻居坐标、单元邻居U
	// 输出：单元梯度
	const int bid = blockIdx.x;
	const int tid = threadIdx.x;
	const int id = bid * blockDim.x + tid;
	const int& i = id;
	if (i >= element_device.num_element || i < 0) return;

	// 剔除无效邻居
	int nValidNeighbor = 0;
	const int neighborArraySize = 3;
	GPU::Element* validNeighbors[neighborArraySize]{};
	for (int j = 0; j < neighborArraySize; j++) {
		if (!element_device.neighbors[j].isNull[i]) {
			validNeighbors[nValidNeighbor] = &(element_device.neighbors[j]);
			nValidNeighbor += 1;
		}
	}

	// 计算梯度dUdX
	GPU::Element* self = &(element_device.self);
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

	// 将dUdX存进单元
	self->Ux1[i] = dUdX[0][0];
	self->Ux2[i] = dUdX[0][1];
	self->Ux3[i] = dUdX[0][2];
	self->Ux4[i] = dUdX[0][3];
	self->Uy1[i] = dUdX[1][0];
	self->Uy2[i] = dUdX[1][1];
	self->Uy3[i] = dUdX[1][2];
	self->Uy4[i] = dUdX[1][3];


	// Barth限制器 一阶精度格式不需要加
	//pE->restructor_in_updateSlope_Barth(f);
	//Limiter::modifySlope_Barth(pE);

	// 检查异常值 应放在外面
	// 只需检查是否有±1e10量级的数

}
;

__device__ void GPU::Limiter::Barth() {
	// TODO: Barth limiter
	// 功能：修正单元梯度
	// 输入：单元U、邻居U、节点U、单元梯度
	// 输出：单元梯度
}
