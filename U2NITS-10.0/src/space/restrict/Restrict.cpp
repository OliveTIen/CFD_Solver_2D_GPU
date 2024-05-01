#include "Restrict.h"
#include "../../math/Math.h"
#include "../../output/LogWriter.h"
#include "../../global/GlobalPara.h"

void U2NITS::Space::Restrict::modifyElementFieldU2d(GPU::ElementSoA& element, GPU::ElementFieldSoA& elementField) {
	// ���ڳ�����Χ�����ݣ�ȡ�ھӵ�ƽ��ֵ
	int nElement = element.num_element;
	for (int iElement = 0; iElement < nElement; iElement++) {
		modifyElementFieldUKernel2d(element, elementField, iElement);
	}
}

void U2NITS::Space::Restrict::modifyElementFieldUKernel2d(GPU::ElementSoA& element, GPU::ElementFieldSoA& elementField, int iElement) {
	
	myfloat U[4]{};
	myfloat ruvp[4]{};
	for (int i = 0; i < 4; i++) {
		U[i] = elementField.U[i][iElement];
	}
	Math::U2ruvp_host(U, ruvp, GlobalPara::constant::gamma);

	if (outOfRange(ruvp)) {
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
				LogWriter::logAndPrintError("no neighbor.\n");
				exit(-1);
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
				LogWriter::logAndPrintError("no neighbor.\n");
				exit(-1);
			}
			// ���ھ�Uƽ��ֵ���µ�Ԫֵ
			for (int jVar = 0; jVar < 4; jVar++) {
				elementField.U[jVar][iElement] = Usum[jVar] / numOfValidNeighbor;
			}	
		}


	}
}


