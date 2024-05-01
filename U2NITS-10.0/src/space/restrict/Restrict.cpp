#include "Restrict.h"
#include "../../math/Math.h"
#include "../../output/LogWriter.h"
#include "../../global/GlobalPara.h"

void U2NITS::Space::Restrict::modifyElementFieldU2d(GPU::ElementSoA& element, GPU::ElementFieldSoA& elementField) {
	// 对于超出范围的数据，取邻居的平均值
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
				LogWriter::logAndPrintError("no neighbor.\n");
				exit(-1);
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
				LogWriter::logAndPrintError("no neighbor.\n");
				exit(-1);
			}
			// 用邻居U平均值更新单元值
			for (int jVar = 0; jVar < 4; jVar++) {
				elementField.U[jVar][iElement] = Usum[jVar] / numOfValidNeighbor;
			}	
		}


	}
}


