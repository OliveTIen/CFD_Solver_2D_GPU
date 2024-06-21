#include "RestrictGPU.h"
#include "../../global/GlobalPara.h"
#include "../../math/PhysicalKernelGPU.h"


__global__ void modifyElementFieldU2d_device_kernel(GPU::ElementSoA element, GPU::ElementFieldSoA elementField, const myfloat gamma) {
    myint iElement = blockDim.x * blockIdx.x + threadIdx.x;
    if (iElement >= element.num_element) return;

    // �غ���תԭʼ����
    myfloat U[4]{};
    myfloat ruvp[4]{};
    for (int i = 0; i < 4; i++) {
        U[i] = elementField.U[i][iElement];
    }
    GPU::Math::U2ruvp_device(U, ruvp, gamma);
    // ����Ƿ���Ҫ����(NaN��Խ��)
    bool need_modification = false;
    for (int i = 0; i < 4; i++) {
        if (isnan(ruvp[i])) {
            need_modification = true;
        }
    }
    if (GPU::Space::Restrict::outOfRange_device(ruvp))need_modification = true;
    // ����Ҫ��������ȡ�ھ�ƽ��ֵ
    if (need_modification) {
		bool volumeWeight = true;// �����Ȩ
		if (volumeWeight) {

			// ���ھӸ������ھ�U֮��
			int numOfValidNeighbor = 0;
			myfloat Usum[4]{};
			myfloat volumeSum = 0.0;
			for (int i = 0; i < 4; i++) {
				// -1��ʾ���ھӣ����ڱ߽���Ȼ��-1
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
			// ���ھ�Uƽ��ֵ���µ�Ԫֵ�������Ȩƽ��
			for (int jVar = 0; jVar < 4; jVar++) {
				elementField.U[jVar][iElement] = Usum[jVar] / volumeSum;
			}
		}
		else {// ����ƽ��

			// ���ھӸ������ھ�U֮��
			int numOfValidNeighbor = 0;
			myfloat Usum[4]{};
			for (int i = 0; i < 4; i++) {
				// -1��ʾ���ھӣ����ڱ߽���Ȼ��-1
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
			// ���ھ�Uƽ��ֵ���µ�Ԫֵ
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
