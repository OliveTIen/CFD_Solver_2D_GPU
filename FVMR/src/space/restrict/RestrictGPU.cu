#include "RestrictGPU.h"
#include "../../global/GlobalPara.h"
#include "../../math/PhysicalKernelGPU.h"


__global__ void modifyElementFieldU2d_device_kernel(GPU::ElementSoA element, GPU::ElementFieldSoA elementField, const myfloat gamma) {
    myint iElement = blockDim.x * blockIdx.x + threadIdx.x;
    if (iElement >= element.num_element) return;

    // 守恒量转原始变量
    myfloat U[4]{};
    myfloat ruvp[4]{};
    for (int i = 0; i < 4; i++) {
        U[i] = elementField.U[i][iElement];
    }
    GPU::Math::U2ruvp_device(U, ruvp, gamma);
    // 检测是否需要修正(NaN、越界)
    bool need_modification = false;
    for (int i = 0; i < 4; i++) {
        if (isnan(ruvp[i])) {
            need_modification = true;
        }
    }
    if (GPU::Space::Restrict::outOfRange_device(ruvp))need_modification = true;
    // 若需要修正，则取邻居平均值
    if (need_modification) {
		bool volumeWeight = true;// 体积加权
		if (volumeWeight) {

			// 求邻居个数和邻居U之和
			int numOfValidNeighbor = 0;
			myfloat Usum[4]{};
			myfloat volumeSum = 0.0;
			for (int i = 0; i < 4; i++) {
				// -1表示无邻居；周期边界依然是-1
				int neighborID = element.neighbors[i][iElement];
				if (neighborID == -1) {
					continue;
				}
				numOfValidNeighbor++;
				for (int jVar = 0; jVar < 4; jVar++) {
					myfloat volume = element.volume[neighborID];
					Usum[jVar] += elementField.U[jVar][neighborID] * volume;
					volumeSum += volume;
				}
			}
			if (numOfValidNeighbor == 0) {
				printf("[error] no neighbor @ modifyElementFieldU2d_device_kernel.\n");
				//exit(-1);
				return;
			}
			// 用邻居U平均值更新单元值，体积加权平均
			for (int jVar = 0; jVar < 4; jVar++) {
				elementField.U[jVar][iElement] = Usum[jVar] / volumeSum;
			}
		}
		else {// 算术平均

			// 求邻居个数和邻居U之和
			int numOfValidNeighbor = 0;
			myfloat Usum[4]{};
			for (int i = 0; i < 4; i++) {
				// -1表示无邻居；周期边界依然是-1
				int neighborID = element.neighbors[i][iElement];
				if (neighborID == -1) {
					continue;
				}
				numOfValidNeighbor++;
				for (int jVar = 0; jVar < 4; jVar++) {
					Usum[jVar] += elementField.U[jVar][neighborID];
				}
			}
			if (numOfValidNeighbor == 0) {
				printf("[error] no neighbor @ modifyElementFieldU2d_device_kernel.\n");
				//exit(-1);
				return;
			}
			// 用邻居U平均值更新单元值
			for (int jVar = 0; jVar < 4; jVar++) {
				elementField.U[jVar][iElement] = Usum[jVar] / numOfValidNeighbor;
			}
		}
    }
}

void GPU::Space::Restrict::modifyElementFieldU2d_device(GPU::ElementSoA& element, GPU::ElementFieldSoA& elementField, const myfloat gamma) {
    myint num_element = element.num_element;
    int threadsPerBlock = 256;
    int blocksPerGrid = (num_element + threadsPerBlock - 1) / threadsPerBlock;
    modifyElementFieldU2d_device_kernel <<<blocksPerGrid, threadsPerBlock>>> (element, elementField, gamma);
    getLastCudaError("modifyElementFieldU2d_device failed.");
}
