#include "Gradient.h"
#include "../../math/Math.h"
#include "../../output/LogWriter.h"
#include "../../global/GlobalPara.h"

void U2NITS::Space::Gradient::Gradient(GPU::ElementSoA& element_host, GPU::FieldSoA& elementField_host, GPU::NodeSoA& node_host, GPU::EdgeSoA& edge_host) {
	switch (GlobalPara::space::flag_gradient) {
	case _GRA_leastSquare:
		//LogWriter::logAndPrint("GradientLeastSquare_2\n");
		GradientLeastSquare_2(element_host, elementField_host, node_host);
		break;
	case _GRA_greenGauss:
		//LogWriter::logAndPrint("GradientGreenGauss\n");
		GradientGreenGauss_2(element_host, elementField_host, edge_host);
		break;

	default:
		LogWriter::logAndPrintError("invalid graident type.\n");
		exit(-1);
		break;
	}
}

void U2NITS::Space::Gradient::GradientLeastSquare_old(GPU::ElementSoA& element_host, GPU::FieldSoA& elementField_host) {
	// 最后把device重命名为host

	for (int i = 0; i < element_host.num_element; i++) {
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
			if (element_host.neighbors[j][i] != -1) {
				validNeighbors[nValidNeighbor] = element_host.neighbors[j][i];
				nValidNeighbor += 1;
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
				dX[iVN][0] = element_host.xy[0][index] - element_host.xy[0][i];
				dX[iVN][1] = element_host.xy[1][index] - element_host.xy[1][i];
				for (int iVar = 0; iVar < nVar; iVar++) {
					//dU[iVN][iVar] = elementField_host.U[index][iVar] - elementField_host.U[i][iVar];
					dU[iVN][iVar] = elementField_host.U[iVar][index] - elementField_host.U[iVar][i];
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
				dX[iVN][0] = element_host.xy[0][index] - element_host.xy[0][i];
				dX[iVN][1] = element_host.xy[1][index] - element_host.xy[1][i];
				for (int iVar = 0; iVar < nVar; iVar++) {
					dU[iVN][iVar] = elementField_host.U[iVar][index] - elementField_host.U[iVar][i];
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
				dX[iVN][0] = element_host.xy[0][index] - element_host.xy[0][i];
				dX[iVN][1] = element_host.xy[1][index] - element_host.xy[1][i];
				for (int iVar = 0; iVar < nVar; iVar++) {
					dU[iVN][iVar] = elementField_host.U[iVar][index] - elementField_host.U[iVar][i];
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
			elementField_host.Ux[j][i] = dUdX[0][j];
			elementField_host.Uy[j][i] = dUdX[1][j];
		}

		// Barth限制器 一阶精度格式不需要加
		//pE->restructor_in_updateSlope_Barth(f);
		//Limiter::modifySlope_Barth(pE);

		// 检查异常值 应放在外面
		// 只需检查是否有±1e10量级的数
	}
}

void U2NITS::Space::Gradient::GradientLeastSquare_2(GPU::ElementSoA& element_host, GPU::FieldSoA& elementField_host, GPU::NodeSoA& node_host) {

	for (int iElement = 0; iElement < element_host.num_element; iElement++) {


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
			LogWriter::logAndPrint("Error: invalid numOfValidNeighbor.\n", LogWriter::Error, LogWriter::Error);
			exit(-1);
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
		U2NITS::Math::Matrix::transpose(nVN, nX, (real*)dX, (real*)dXtrans);
		U2NITS::Math::Matrix::mul_ixj_jxk(nX, nVN, nX, (real*)dXtrans, (real*)dX, (real*)invdXdX);
		U2NITS::Math::Matrix::inv_2x2(invdXdX);
		U2NITS::Math::Matrix::mul_ixj_jxk(nX, nX, nVN, (real*)invdXdX, (real*)dXtrans, (real*)invdXdX_dX);
		U2NITS::Math::Matrix::mul_ixj_jxk(nX, nVN, nVar, (real*)invdXdX_dX, (real*)dU, (real*)dUdX);

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
		U2NITS::Math::Matrix::mul_ixj_jxk(nNode, nX, nVar, (real*)dX_node, dUdX, (real*)dU_node);
		real ratio = 1.0;
		for (int iNode = 0; iNode < nNode; iNode++) {
			for (int iVar = 0; iVar < nVar; iVar++) {
				using namespace U2NITS::Math;
				if (dU[iNode * nVar + iVar] == 0)continue;// 除零
				ratio = max(ratio, abs(dU_node[iNode][iVar] / dU[iNode * nVar + iVar]));
			}
		}
		// ratio一定大于等于1.0，因此可以直接除以ratio
		U2NITS::Math::Matrix::div_matrix_by_scalar(nX, nVar, dUdX, ratio);
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
}

void U2NITS::Space::Gradient::GradientGreenGauss_1(GPU::ElementSoA& element_host, GPU::FieldSoA& elementField_host, GPU::EdgeSoA& edge_host) {
	/*
	* 旧版格林高斯梯度，边界元有问题
	* https://zhuanlan.zhihu.com/p/370586072
	* 将面矢量加到单元时，应注意面是朝里还是朝外。
	* 注意计算通量时排除了周期边界，会有隐患
	* 对于边界，其通量应等于phiC(参见01-01-CFD理论)
	* 对于非法边(三角形单元的第4条边)，edgeID=-1, outElementID=-1, 
	*/
	// 求梯度的子迭代步数
	int subIterationGradient = 3;
	int num_edge = edge_host.num_edge;

	for (int iElement = 0; iElement < element_host.num_element; iElement++) {
		// 求当前单元ID、当前单元体积、edgeID、edge朝向正负、外单元ID
		int currentElementID = element_host.ID[iElement];// 建议不要直接用iElement，虽然目前iElement==i(element_host.ID[i]=element_i.GPUID=i)，但有隐患
		real volumeC = element_host.volume[currentElementID];// 体积
		int edgeID[4]{};
		int edgeSign[4]{};// 1表示边朝外，-1表示边朝内
		int edgeType[4]{ -1,-1,-1,-1 };
		for (int i = 0; i < 4; i++) {
			
			edgeID[i] = element_host.edges[i][iElement];
			// 非法值，赋为-1，防止数组访问越界。后面遇到-1会跳过
			if (edgeID[i] <= -1 || edgeID[i] >= num_edge) {
				edgeID[i] = -1;// 对于三角形单元，其第4条边的ID值是不确定的。超出范围则赋为-1
			}
			else if (edge_host.setID[edgeID[i]] != -1) {
				edgeType[i] = edge_host.setID[edgeID[i]];// 存储edge类型，-1表示内部边
				edgeID[i] = -1;// 非内部边，默认排除。用于排除边界edge
			}
		}
		int outElementID[4]{ -1,-1,-1,-1 };//外单元ID 注意周期边界ID!=-1
		for (int i = 0; i < 4; i++) {
			int currentEdgeID = edgeID[i];
			if (currentEdgeID == -1) {
				continue;// 默认是-1。若edgeID是-1，则outElementID一定是-1
			}
			// 确定外单元ID和朝向。对于内部边，无法确定是elementL还是elementR
			if (currentElementID != edge_host.elementR[currentEdgeID]) {
				outElementID[i] = edge_host.elementR[currentEdgeID];
				edgeSign[i] = 1;// currentElement=elementL，朝外
			}
			else {
				outElementID[i] = edge_host.elementL[currentEdgeID];
				edgeSign[i] = -1;// currentElement=elementR，朝内
			}
		}
		// 求面积矢量(几何量)
		real edgeNormal[2][4]{};
		real edgeArea[4]{};
		real edgeSVector[2][4]{};// 面积矢量
		for (int i = 0; i < 4; i++) {
			int currentEdgeID = edgeID[i];
			if (currentEdgeID == -1) {
				continue;
			}
			edgeNormal[0][i] = edge_host.normal[0][currentEdgeID];
			edgeNormal[1][i] = edge_host.normal[1][currentEdgeID];
			edgeArea[i] = edge_host.length[currentEdgeID];
			edgeSVector[0][i] = edgeNormal[0][i] * edgeArea[i];
			edgeSVector[1][i] = edgeNormal[1][i] * edgeArea[i];
		}

		// 针对每条边，求几何权重因子gc、矢径(几何量)
		real gC_all[4]{};
		real rf_rC_all[4][2]{};
		real rf_rF_all[4][2]{};
		for (int i = 0; i < 4; i++) {
			int currentEdgeID = edgeID[i];
			int currentOutElementID = outElementID[i];
			if (currentEdgeID == -1) {
				continue;
			}

			// 面和单元坐标
			real fx = edge_host.xy[0][currentEdgeID];// face x
			real fy = edge_host.xy[1][currentEdgeID];
			real Cx = element_host.xy[0][currentElementID];// cell x
			real Cy = element_host.xy[1][currentElementID];
			real Fx = element_host.xy[0][currentOutElementID];
			real Fy = element_host.xy[1][currentOutElementID];
			using namespace U2NITS;
			real dFf = Math::distance2D(Fx, Fy, fx, fy);
			real dFC = Math::distance2D(Fx, Fy, Cx, Cy);
			if (dFC < Math::EPSILON)dFC = Math::EPSILON;

			real gC = dFf / dFC;// 几何权重因子gc
			real rf_rC[2]{ fx - Cx,fy - Cy };// rf-rC
			real rf_rF[2]{ fx - Fx,fy - Fy };// rf-rF

			gC_all[i] = gC;
			rf_rC_all[i][0] = rf_rC[0];
			rf_rC_all[i][1] = rf_rC[1];
			rf_rF_all[i][0] = rf_rF[0];
			rf_rF_all[i][1] = rf_rF[1];
		}

		// 计算每条边的向量通量phif·Sf，然后求和，除以体积，得到单元C的梯度
		const int nVar = 4;
		for (int jVar = 0; jVar < nVar; jVar++) {
			// 求和。若面朝内，该面通量取相反数
			real sum_phifSf[2]{};
			for (int i = 0; i < 4; i++) {
				int currentEdgeID = edgeID[i];
				int currentOutElementID = outElementID[i];
				if (currentEdgeID == -1) {
					continue;
				}

				// 获取单元值和单元梯度值
				real phiC = elementField_host.U[jVar][currentElementID];
				real phiF = elementField_host.U[jVar][currentOutElementID];
				//real nabla_phiC[2]{// 迭代时才会用到
				//	elementField_host.Ux[jVar][currentElementID],
				//	elementField_host.Uy[jVar][currentElementID]
				//};// 梯度
				//real nabla_phiF[2]{
				//	elementField_host.Ux[jVar][currentOutElementID],
				//	elementField_host.Uy[jVar][currentOutElementID]
				//};
				// 求界面值(近似)
				real gC = gC_all[i];
				real gC1 = 1.0 - gC;
				real phif = gC * phiC + gC1 * phiF;
				real Sf[2]{ edgeSVector[0][i],edgeSVector[1][i] };
				sum_phifSf[0] += phif * Sf[0] * edgeSign[i];// 乘以符号，以应对面朝内的情形
				sum_phifSf[1] += phif * Sf[1] * edgeSign[i];
			}

			real nabla_phiC_new[2]{// 梯度
				1 / volumeC * sum_phifSf[0],
				1 / volumeC * sum_phifSf[1]
			};
			// 更新单元梯度值
			elementField_host.Ux[jVar][currentElementID] = nabla_phiC_new[0];
			elementField_host.Uy[jVar][currentElementID] = nabla_phiC_new[1];
			// 如果要迭代，需要等所有单元更新完毕，因此迭代循环要放在iElement循环外
		}
	}
}

void U2NITS::Space::Gradient::GradientGreenGauss_2(GPU::ElementSoA& element_host, GPU::FieldSoA& elementField_host, GPU::EdgeSoA& edge_host) {
	/*
	2024-03-28

	* 格林高斯梯度
	* https://zhuanlan.zhihu.com/p/370586072
	* 将面矢量加到单元时，应注意面是朝里还是朝外。
	* 注意计算通量时排除了周期边界，会有隐患
	* 对于边界，其通量应等于phiC(参见01-01-CFD理论)
	* 对于非法边(三角形单元的第4条边)，edgeID=-1, outElementID=-1,
	*/
// 求梯度的子迭代步数
	int subIterationGradient = 3;
	int num_edge = edge_host.num_edge;

	for (int iElement = 0; iElement < element_host.num_element; iElement++) {
		// 求当前单元ID、当前单元体积、edgeID、edge朝向正负、外单元ID
		int currentElementID = element_host.ID[iElement];// 建议不要直接用iElement，虽然目前iElement==i(element_host.ID[i]=element_i.GPUID=i)，但有隐患
		real volumeC = element_host.volume[currentElementID];// 体积
		// 初始化edgeID 超出数组范围的非法值赋为-1
		int edgeID[4]{ -1,-1,-1,-1 };// 第iElement个element的第i条边的ID
		for (int i = 0; i < 4; i++) {
			edgeID[i] = element_host.edges[i][iElement];
			if (edgeID[i] <= -1 || edgeID[i] >= num_edge) {
				edgeID[i] = -1;
			}
		}
		// 初始化edgeType
		int edgeType[4]{ -1,-1,-1,-1 };
		for (int i = 0; i < 4; i++) {
			int edgeIDi = edgeID[i];
			if (edgeIDi == -1) {
				continue;
			}

			edgeType[i] = edge_host.setID[edgeIDi];
		}
		// 初始化edgeSign和outElementID
		int edgeSign[4]{ 0,0,0,0 };// 1表示边朝外，-1表示边朝内
		int outElementID[4]{ -1,-1,-1,-1 };//外单元ID 内部边和周期边
		for (int i = 0; i < 4; i++) {
			int edgeIDi = edgeID[i];
			if (edgeIDi == -1) {
				continue;
			}

			if (currentElementID != edge_host.elementR[edgeIDi]) {
				edgeSign[i] = 1;// currentElement=elementL，朝外
				outElementID[i] = edge_host.elementR[edgeIDi];
			}
			else {
				edgeSign[i] = -1;// currentElement=elementR，朝内
				outElementID[i] = edge_host.elementL[edgeIDi];
			}
		}
		// 求面积矢量edgeSVector(几何量)
		real edgeNormal[2][4]{};
		real edgeArea[4]{};
		real edgeSVector[2][4]{};// 面积矢量
		for (int i = 0; i < 4; i++) {
			int edgeIDi = edgeID[i];
			if (edgeIDi == -1) {
				continue;
			}

			edgeNormal[0][i] = edge_host.normal[0][edgeIDi];
			edgeNormal[1][i] = edge_host.normal[1][edgeIDi];
			edgeArea[i] = edge_host.length[edgeIDi];
			edgeSVector[0][i] = edgeNormal[0][i] * edgeArea[i];
			edgeSVector[1][i] = edgeNormal[1][i] * edgeArea[i];
		}

		// 求几何权重因子gc、矢径(几何量)
		real gC_all[4]{};
		real rf_rC_all[4][2]{};
		real rf_rF_all[4][2]{};
		for (int i = 0; i < 4; i++) {
			int edgeIDi = edgeID[i];
			if (edgeIDi == -1) {
				continue;
			}

			int currentOutElementID = outElementID[i];
			if (currentOutElementID != -1){
				// 内部边和周期边
				if (edgeType[i] == -1) {
					// 内部边

					// 面和单元坐标
					real fx = edge_host.xy[0][edgeIDi];// face x
					real fy = edge_host.xy[1][edgeIDi];
					real Cx = element_host.xy[0][currentElementID];// cell x
					real Cy = element_host.xy[1][currentElementID];
					real Fx = element_host.xy[0][currentOutElementID];
					real Fy = element_host.xy[1][currentOutElementID];
					// 几何权重因子gC = |rF-rf|/|rF-rC|，存储于gC_all数组中
					using namespace U2NITS;
					real dFf = Math::distance2D(Fx, Fy, fx, fy);
					real dFC = Math::distance2D(Fx, Fy, Cx, Cy);
					if (dFC < Math::EPSILON)dFC = Math::EPSILON;
					gC_all[i] = dFf / dFC;
					// 矢径 rf-rC rf-rF，在迭代时会用到
					rf_rC_all[i][0] = fx - Cx;// rf-rC
					rf_rC_all[i][1] = fy - Cy;
					rf_rF_all[i][0] = fx - Fx;// rf-rF
					rf_rF_all[i][1] = fy - Fy;
				}
				else {
					// 周期边
					// 面和单元坐标
					real fx = edge_host.xy[0][edgeIDi];
					real fy = edge_host.xy[1][edgeIDi];
					real Cx = element_host.xy[0][currentElementID];// cell x
					real Cy = element_host.xy[1][currentElementID];
					// 周期边的outElement坐标需要平移
					int pairIDi = edge_host.periodicPair[edgeIDi];
					real fx_pair = edge_host.xy[0][pairIDi];
					real fy_pair = edge_host.xy[1][pairIDi];
					real shiftx = fx - fx_pair;// 从pair指向当前位置的矢径
					real shifty = fy - fy_pair;
					real Fx = element_host.xy[0][currentOutElementID];
					real Fy = element_host.xy[1][currentOutElementID];
					Fx += shiftx;// 用该矢径平移单元
					Fy += shifty;

					// 后面和内部边是一样的
					// 几何权重因子gC = |rF-rf|/|rF-rC|，存储于gC_all数组中
					using namespace U2NITS;
					real dFf = Math::distance2D(Fx, Fy, fx, fy);
					real dFC = Math::distance2D(Fx, Fy, Cx, Cy);
					if (dFC < Math::EPSILON)dFC = Math::EPSILON;
					gC_all[i] = dFf / dFC;
					// 矢径 rf-rC rf-rF，在迭代时会用到
					rf_rC_all[i][0] = fx - Cx;// rf-rC
					rf_rC_all[i][1] = fy - Cy;
					rf_rF_all[i][0] = fx - Fx;// rf-rF
					rf_rF_all[i][1] = fy - Fy;
				}

			}
			else {
				// 边界边(周期边界除外)
				// F单元实际不存在，假设F单元和C单元关于f面中心对称
				real fx = edge_host.xy[0][edgeIDi];// face x
				real fy = edge_host.xy[1][edgeIDi];
				real Cx = element_host.xy[0][currentElementID];// cell x
				real Cy = element_host.xy[1][currentElementID];

				using namespace U2NITS;
				real gC = 0.5;// 几何权重因子gc
				real rf_rC[2]{ fx - Cx,fy - Cy };// rf-rC
				real rf_rF[2]{ -(fx - Cx),-(fy - Cy)};// rf-rF

				gC_all[i] = gC;
				rf_rC_all[i][0] = rf_rC[0];
				rf_rC_all[i][1] = rf_rC[1];
				rf_rF_all[i][0] = rf_rF[0];
				rf_rF_all[i][1] = rf_rF[1];
			}

		}



		// 计算每条边的向量通量phif·Sf，然后求和，除以体积，得到单元C的梯度
		const int nVar = 4;
		for (int jVar = 0; jVar < nVar; jVar++) {
			// 求和。若面朝内，该面通量取相反数
			real sum_phifSf[2]{};
			for (int i = 0; i < 4; i++) {
				int currentEdgeID = edgeID[i];
				if (currentEdgeID == -1) {
					continue;
				}
				int currentOutElementID = outElementID[i];

				if (currentOutElementID != -1) {
					// 内部边和周期边

					// 获取单元值和单元梯度值
					real phiC = elementField_host.U[jVar][currentElementID];
					real phiF = elementField_host.U[jVar][currentOutElementID];
					// 求界面值(近似)
					real gC = gC_all[i];
					real gC1 = 1.0 - gC;
					real phif = gC * phiC + gC1 * phiF;
					real Sf[2]{ edgeSVector[0][i],edgeSVector[1][i] };
					sum_phifSf[0] += phif * Sf[0] * edgeSign[i];// 乘以符号，以应对面朝内的情形
					sum_phifSf[1] += phif * Sf[1] * edgeSign[i];
				}
				else {
					// 边界边(周期边除外)
					real phiC = elementField_host.U[jVar][currentElementID];
					real phif = phiC;
					real Sf[2]{ edgeSVector[0][i],edgeSVector[1][i] };
					sum_phifSf[0] += phif * Sf[0] * edgeSign[i];// 乘以符号，以应对面朝内的情形
					sum_phifSf[1] += phif * Sf[1] * edgeSign[i];
				}

			}

			real nabla_phiC_new[2]{// 梯度
				1 / volumeC * sum_phifSf[0],
				1 / volumeC * sum_phifSf[1]
			};
			// 更新单元梯度值
			elementField_host.Ux[jVar][currentElementID] = nabla_phiC_new[0];
			elementField_host.Uy[jVar][currentElementID] = nabla_phiC_new[1];
			// 如果要迭代，需要等所有单元更新完毕，因此迭代循环要放在iElement循环外
		}
	}
}


