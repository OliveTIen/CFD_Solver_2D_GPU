#include "Gradient.h"
#include "../../math/MatrixKernel.h"

void U2NITS::Space::Gradient::GradientLeastSquare(GPU::ElementSoA& element_device, GPU::FieldSoA elementField_device) {
	// 最后把device重命名为host

	for (int i = 0; i < element_device.num_element; i++) {
		//const int bid = blockIdx.x;
		//const int tid = threadIdx.x;
		//const int id = bid * blockDim.x + tid;
		//const int& i = id;
		//if (i >= element_device.num_element || i < 0) return;

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
			U2NITS::Math::Matrix::transpose(nVN, nX, (REAL*)dX, (REAL*)dXtrans);
			U2NITS::Math::Matrix::mul_ixj_jxk(nX, nVN, nX, (REAL*)dXtrans, (REAL*)dX, (REAL*)invdXdX);
			U2NITS::Math::Matrix::inv_2x2(invdXdX);
			U2NITS::Math::Matrix::mul_ixj_jxk(nX, nX, nVN, (REAL*)invdXdX, (REAL*)dXtrans, (REAL*)invdXdX_dX);
			U2NITS::Math::Matrix::mul_ixj_jxk(nX, nVN, nVar, (REAL*)invdXdX_dX, (REAL*)dU, (REAL*)dUdX);
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
			U2NITS::Math::Matrix::transpose(nVN, nX, (REAL*)dX, (REAL*)dXtrans);
			U2NITS::Math::Matrix::mul_ixj_jxk(nX, nVN, nX, (REAL*)dXtrans, (REAL*)dX, (REAL*)invdXdX);
			U2NITS::Math::Matrix::inv_2x2(invdXdX);
			U2NITS::Math::Matrix::mul_ixj_jxk(nX, nX, nVN, (REAL*)invdXdX, (REAL*)dXtrans, (REAL*)invdXdX_dX);
			U2NITS::Math::Matrix::mul_ixj_jxk(nX, nVN, nVar, (REAL*)invdXdX_dX, (REAL*)dU, (REAL*)dUdX);
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
			U2NITS::Math::Matrix::transpose(nVN, nX, (REAL*)dX, (REAL*)dXtrans);
			U2NITS::Math::Matrix::mul_ixj_jxk(nX, nVN, nX, (REAL*)dXtrans, (REAL*)dX, (REAL*)invdXdX);
			U2NITS::Math::Matrix::inv_2x2(invdXdX);
			U2NITS::Math::Matrix::mul_ixj_jxk(nX, nX, nVN, (REAL*)invdXdX, (REAL*)dXtrans, (REAL*)invdXdX_dX);
			U2NITS::Math::Matrix::mul_ixj_jxk(nX, nVN, nVar, (REAL*)invdXdX_dX, (REAL*)dU, (REAL*)dUdX);
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
}