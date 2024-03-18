
#include "CalculateGradientGPU.h"
#include "../../solvers/GPUSolver2.h"
#include "../../math/GPUMatrixKernel.h"
#include "../../output/LogWriter.h"

void GPU::calculateGradient2(GPU::ElementSoA& element_device, GPU::FieldSoA elementField_device) {
	int block_size = 512;// �����128 256 512
	int grid_size = (element_device.num_element + block_size - 1) / block_size;
	dim3 block(block_size, 1, 1);
	dim3 grid(grid_size, 1, 1);
	calculateGradientKernel2 <<<grid, block >>> (element_device, elementField_device);
	cudaError_t cuda_error = cudaGetLastError();
	if (cuda_error != 0) {
		std::string e = "cudaError=" + std::to_string(cuda_error) + ", " + cudaGetErrorString(cuda_error);
		LogWriter::writeLogAndCout(e, LogWriter::Error, LogWriter::Error);
	}
	// cudaDeviceSynchronize();��������
}

__global__ void GPU::calculateGradientKernel2(GPU::ElementSoA& element_device, GPU::FieldSoA elementField_device) {
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
		GPU::Matrix::transpose(nVN, nX, (REAL*)dX, (REAL*)dXtrans);
		GPU::Matrix::mul_ixj_jxk(nX, nVN, nX, (REAL*)dXtrans, (REAL*)dX, (REAL*)invdXdX);
		GPU::Matrix::inv_2x2(invdXdX);
		GPU::Matrix::mul_ixj_jxk(nX, nX, nVN, (REAL*)invdXdX, (REAL*)dXtrans, (REAL*)invdXdX_dX);
		GPU::Matrix::mul_ixj_jxk(nX, nVN, nVar, (REAL*)invdXdX_dX, (REAL*)dU, (REAL*)dUdX);
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
		GPU::Matrix::transpose(nVN, nX, (REAL*)dX, (REAL*)dXtrans);
		GPU::Matrix::mul_ixj_jxk(nX, nVN, nX, (REAL*)dXtrans, (REAL*)dX, (REAL*)invdXdX);
		GPU::Matrix::inv_2x2(invdXdX);
		GPU::Matrix::mul_ixj_jxk(nX, nX, nVN, (REAL*)invdXdX, (REAL*)dXtrans, (REAL*)invdXdX_dX);
		GPU::Matrix::mul_ixj_jxk(nX, nVN, nVar, (REAL*)invdXdX_dX, (REAL*)dU, (REAL*)dUdX);
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
		GPU::Matrix::transpose(nVN, nX, (REAL*)dX, (REAL*)dXtrans);
		GPU::Matrix::mul_ixj_jxk(nX, nVN, nX, (REAL*)dXtrans, (REAL*)dX, (REAL*)invdXdX);
		GPU::Matrix::inv_2x2(invdXdX);
		GPU::Matrix::mul_ixj_jxk(nX, nX, nVN, (REAL*)invdXdX, (REAL*)dXtrans, (REAL*)invdXdX_dX);
		GPU::Matrix::mul_ixj_jxk(nX, nVN, nVar, (REAL*)invdXdX_dX, (REAL*)dU, (REAL*)dUdX);
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