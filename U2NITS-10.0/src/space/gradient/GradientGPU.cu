//#include <stdio.h>
#include "GradientGPU.h"
#include "../../math/MathGPU.h"
#include "../../output/LogWriter.h"
#include "../../solvers/GPUSolver2.h"

void GPU::calculateGradient_old(GPU::ElementSoA& element_device, GPU::FieldSoA& elementField_device) {
	int block_size = 512;// 最好是128 256 512
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
	// cudaDeviceSynchronize();放在外面
}

__global__ void GPU::calculateGradientKernel_old(GPU::ElementSoA& element_device, GPU::FieldSoA& elementField_device) {
	// --- 计算单元梯度核函数 已完成 --- 
	// 输入：单元坐标、单元U、单元邻居坐标、单元邻居U
	// 输出：单元梯度
	const int bid = blockIdx.x;
	const int tid = threadIdx.x;
	const int id = bid * blockDim.x + tid;
	const int& i = id;
	if (i >= element_device.num_element || i < 0) return;

	// 提取有效邻居
	int nValidNeighbor = 0;
	int validNeighbors[3]{ -1,-1,-1 };// 有效邻居的ID
	const int neighborArraySize = 3;// 三角形单元，最多3个邻居
	for (int j = 0; j < neighborArraySize; j++) {
		if (element_device.neighbors[j][i] != -1) {
			nValidNeighbor += 1;
			validNeighbors[nValidNeighbor] = element_device.neighbors[j][i];
		}
	}

	// 计算梯度dUdX
	const int nVar = 4;// 守恒量个数，二维为4
	const int nX = 2;// 坐标个数，二维为2
	//const int nVNmax = 3;// 最大邻居个数
	REAL dUdX[nX][nVar]{};// 已初始化为0
	if (nValidNeighbor == 3) {
		// 初始化dX、dU
		const int nVN = 3;// 有效邻居个数
		REAL dX[nVN][nX]{};
		REAL dU[nVN][nVar]{};
		for (int iVN = 0; iVN < nVN; iVN++) {
			int index = validNeighbors[iVN];// 邻居的ID
			dX[iVN][0] = element_device.xy[0][index] - element_device.xy[0][i];
			dX[iVN][1] = element_device.xy[1][index] - element_device.xy[1][i];
			for (int iVar = 0; iVar < nVar; iVar++) {
				dU[iVN][iVar] = elementField_device.U[index][iVar] - elementField_device.U[i][iVar];
			}
		}

		// 最小二乘法 x=(A'A)^{-1}A'b，dUdX = inv(dXtrans * dX) * dXtrans * dU
		// ! 以下函数为device，应在GPU实现
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
		// 以下代码与上面的代码相同(仅仅nVN不同)，修改时请一同修改
		// 初始化dX、dU
		const int nVN = 2;// 有效邻居个数
		REAL dX[nVN][nX]{};
		REAL dU[nVN][nVar]{};
		for (int iVN = 0; iVN < nVN; iVN++) {
			int index = validNeighbors[iVN];// 邻居的ID
			dX[iVN][0] = element_device.xy[0][index] - element_device.xy[0][i];
			dX[iVN][1] = element_device.xy[1][index] - element_device.xy[1][i];
			for (int iVar = 0; iVar < nVar; iVar++) {
				dU[iVN][iVar] = elementField_device.U[index][iVar] - elementField_device.U[i][iVar];
			}
		}

		// 最小二乘法 x=(A'A)^{-1}A'b，dUdX = inv(dXtrans * dX) * dXtrans * dU
		// ! 以下函数为device，应在GPU实现
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
		// 以下代码与上面的代码相同(仅仅nVN不同)，修改时请一同修改
		// 初始化dX、dU
		const int nVN = 1;// 有效邻居个数
		REAL dX[nVN][nX]{};
		REAL dU[nVN][nVar]{};
		for (int iVN = 0; iVN < nVN; iVN++) {
			int index = validNeighbors[iVN];// 邻居的ID
			dX[iVN][0] = element_device.xy[0][index] - element_device.xy[0][i];
			dX[iVN][1] = element_device.xy[1][index] - element_device.xy[1][i];
			for (int iVar = 0; iVar < nVar; iVar++) {
				dU[iVN][iVar] = elementField_device.U[index][iVar] - elementField_device.U[i][iVar];
			}
		}

		// 最小二乘法 x=(A'A)^{-1}A'b，dUdX = inv(dXtrans * dX) * dXtrans * dU
		// ! 以下函数为device，应在GPU实现
		REAL dXtrans[nX][nVN]{};
		REAL invdXdX[nX][nX]{};
		REAL invdXdX_dX[nX][nVN]{};
		GPU::Math::Matrix::transpose(nVN, nX, (REAL*)dX, (REAL*)dXtrans);
		GPU::Math::Matrix::mul_ixj_jxk(nX, nVN, nX, (REAL*)dXtrans, (REAL*)dX, (REAL*)invdXdX);
		GPU::Math::Matrix::inv_2x2(invdXdX);
		GPU::Math::Matrix::mul_ixj_jxk(nX, nX, nVN, (REAL*)invdXdX, (REAL*)dXtrans, (REAL*)invdXdX_dX);
		GPU::Math::Matrix::mul_ixj_jxk(nX, nVN, nVar, (REAL*)invdXdX_dX, (REAL*)dU, (REAL*)dUdX);
	}

	// 将dUdX存进单元
	for (int j = 0; j < nVar; j++) {
		elementField_device.Ux[j][i] = dUdX[0][j];
		elementField_device.Uy[j][i] = dUdX[1][j];
	}

	// Barth限制器 一阶精度格式不需要加
	//pE->restructor_in_updateSlope_Barth(f);
	//Limiter::modifySlope_Barth(pE);

	// 检查异常值 应放在外面
	// 只需检查是否有±1e10量级的数
}

void GPU::Space::Gradient::GradientLeastSquare(GPU::ElementSoA& element_device, GPU::FieldSoA& elementField_device, GPU::NodeSoA& node_device) {
	int block_size = 512;// 最好是128 256 512
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
	
	// 提取有效邻居
	int numOfValidNeighbor = 0;
	int validNeighbors[3]{ -1,-1,-1 };// 有效邻居的ID
	const int neighborArraySize = 3;// 三角形单元，最多3个邻居
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
	// 申请堆数组
	const int nVar = 4;// 守恒量个数，二维为4
	const int nX = 2;// 坐标个数，二维为2
	const int& nVN = numOfValidNeighbor;// 有效邻居个数，取值1-3
	real* dX = new real[nVN * nX]{};
	real* dUdX = new real[nX * nVar]{};
	real* dU = new real[nVN * nVar]{};
	real* dXtrans = new real[nX * nVN]{};
	real* invdXdX = new real[nX * nX]{};
	real* invdXdX_dX = new real[nX * nVN]{};
	// 初始化dX、dU
	for (int iVN = 0; iVN < nVN; iVN++) {
		int index = validNeighbors[iVN];// 邻居的ID
		dX[iVN * nX + 0] = element_host.xy[0][index] - element_host.xy[0][iElement];
		dX[iVN * nX + 1] = element_host.xy[1][index] - element_host.xy[1][iElement];
		for (int iVar = 0; iVar < nVar; iVar++) {
			dU[iVN * nVar + iVar] = elementField_host.U[iVar][index] - elementField_host.U[iVar][iElement];
		}
	}
	// 最小二乘法计算dUdX，x=(A'A)^{-1}A'b，dUdX = inv(dXtrans * dX) * dXtrans * dU
	GPU::Math::Matrix::transpose(nVN, nX, (real*)dX, (real*)dXtrans);
	GPU::Math::Matrix::mul_ixj_jxk(nX, nVN, nX, (real*)dXtrans, (real*)dX, (real*)invdXdX);
	GPU::Math::Matrix::inv_2x2(invdXdX);
	GPU::Math::Matrix::mul_ixj_jxk(nX, nX, nVN, (real*)invdXdX, (real*)dXtrans, (real*)invdXdX_dX);
	GPU::Math::Matrix::mul_ixj_jxk(nX, nVN, nVar, (real*)invdXdX_dX, (real*)dU, (real*)dUdX);

	/*
限制器
一维下，minmod和vanleer限制器的目的是保证两侧界面上物理量的值不超过相邻单元的物理量的值
参见计算流体力学（李启兵）课件合集P106

类似地，可以保证三角形顶点处U不超过相邻单元U
根据当前斜率Ux, Uy，计算顶点相对值dU_node，前面已经算出邻居中心相对值dU
计算ratio = max(abs(dU_node/dU))，若大于1，则整体除以ratio
		*/
	const int nNode = 3;
	real dU_node[nNode][nVar]{};
	real dX_node[nNode][nX]{};
	// 计算dX_node
	for (int iNode = 0; iNode < nNode; iNode++) {
		for (int iVar = 0; iVar < nVar; iVar++) {
			int nodeID = element_host.nodes[iNode][iElement];
			dX_node[iNode][0] = node_host.xy[0][nodeID] - element_host.xy[0][iElement];
			dX_node[iNode][1] = node_host.xy[1][nodeID] - element_host.xy[1][iElement];
		}
	}
	// 计算dU_node[nNode,nVar] =  dX_node[nNode,nX] * dUdX[nX,nVar]，然后计算ratio
	GPU::Math::Matrix::mul_ixj_jxk(nNode, nX, nVar, (real*)dX_node, dUdX, (real*)dU_node);
	real ratio = 1.0;
	for (int iNode = 0; iNode < nNode; iNode++) {
		for (int iVar = 0; iVar < nVar; iVar++) {
			using namespace GPU::Math;
			if (dU[iNode * nVar + iVar] == 0)continue;// 除零
			ratio = max(ratio, abs(dU_node[iNode][iVar] / dU[iNode * nVar + iVar]));
		}
	}
	// ratio一定大于等于1.0，因此可以直接除以ratio
	GPU::Math::Matrix::div_matrix_by_scalar(nX, nVar, dUdX, ratio);
	// 将dUdX存进单元
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
