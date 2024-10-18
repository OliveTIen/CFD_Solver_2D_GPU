#include "Gradient.h"
#include "../../math/Math.h"
#include "../../output/LogWriter.h"
#include "../../global/GlobalPara.h"

// 函数声明
namespace U2NITS {
	namespace Space {
		namespace Gradient {
			void GradientLeastSquare_2(GPU::ElementSoA& element_host, GPU::ElementFieldSoA& elementField_host, GPU::NodeSoA& node_host);
			void GradientGreenGauss_2(GPU::ElementSoA& element_host, GPU::ElementFieldSoA& elementField_host, GPU::EdgeSoA& edge_host);
		}

	}
}

void gradient_to_Matrix_2x4(myint iElement, myfloat* Ux[4], myfloat* Uy[4], myfloat* dUdX) {
	/*
	将梯度数据转化为2x4矩阵形式
	Ux: 4xn, n = num_element
	Uy: 4xn
	dUdX: 2x4
	*/
	dUdX[0] = Ux[0][iElement];
	dUdX[1] = Ux[1][iElement];
	dUdX[2] = Ux[2][iElement];
	dUdX[3] = Ux[3][iElement];
	dUdX[4] = Uy[0][iElement];
	dUdX[5] = Uy[1][iElement];
	dUdX[6] = Uy[2][iElement];
	dUdX[7] = Uy[3][iElement];

}

inline void getElementT3ValidNeighbor(int validNeighbors[3], int& numOfValidNeighbor, myint iElement, GPU::ElementSoA& element_host) {
	/*
	获取有效邻居
	输出型参数：validNeighbors有效邻居，numOfValidNeighbor有效邻居个数
	*/
	numOfValidNeighbor = 0;
	for (int j = 0; j < 3; j++) {// 三角形单元，最多3个邻居
		if (element_host.neighbors[j][iElement] != -1) {
			validNeighbors[numOfValidNeighbor] = element_host.neighbors[j][iElement];
			numOfValidNeighbor += 1;
		}
	}
}


inline void get_U_linear(myfloat x, myfloat y, myfloat& U_dist, myint i_e, const myfloat* x_e, const myfloat* y_e, const myfloat* U_e, const myfloat* Ux_e, const myfloat* Uy_e) {
	/*
	i_e: 单元index
	输出型参数：U_dist
	*/
	U_dist = U_e[i_e] + Ux_e[i_e] * (x - x_e[i_e]) + Uy_e[i_e] * (y - y_e[i_e]);
}

void LimiterBarth(GPU::ElementFieldSoA& elementField_host, GPU::ElementSoA& element_host, GPU::NodeSoA& node_host) {
	/*
	Barth限制器，用于修正梯度
	对于二维非结构网格，要保证三角形顶点处U不超过相邻单元U
	公式参见论文 The design and application of upwind schemes on unstructured meshes

	创建8个OpenMP线程后，确实CPU占用达到56%，然而计算速度反而下降了，从138.036变为124.603step/s
	*/
//#pragma omp parallel for num_threads(4)
	for (myint iElement = 0; iElement < element_host.num_element; iElement++) {
		int validNeighbors[3]{ -1,-1,-1 };
		int numOfValidNeighbor = 0;
		getElementT3ValidNeighbor(validNeighbors, numOfValidNeighbor, iElement, element_host);
		const auto& U_element = elementField_host.U;
		auto& Ux_element = elementField_host.Ux;
		auto& Uy_element = elementField_host.Uy;

		for (int iVar = 0; iVar < 4; iVar++) {
			// 获取自身及邻居单元中心值的下界U_lower和上界U_upper
			myfloat U_element_center = U_element[iVar][iElement];
			myfloat U_lower = U_element_center;
			myfloat U_upper = U_element_center;
			for (int j = 0; j < numOfValidNeighbor; j++) {
				myint iNeighbor = validNeighbors[j];
				U_lower = U2NITS::Math::min(U_lower, U_element[iVar][iNeighbor]);
				U_upper = U2NITS::Math::max(U_upper, U_element[iVar][iNeighbor]);
			}

			// 计算限制器系数phi，取顶点限制器系数phi_node的最小值
			myfloat phi = 1.0;
			for (int j = 0; j < 3; j++) {
				myint iNode = element_host.nodes[j][iElement];
				myfloat x_node = node_host.xy[0][iNode];
				myfloat y_node = node_host.xy[1][iNode];
				myfloat U_node = 0.0f;
				get_U_linear(x_node, y_node, U_node, iElement, element_host.xy[0], element_host.xy[1], U_element[iVar], Ux_element[iVar], Uy_element[iVar]);
				myfloat U_node_relative = U_node - U_element_center;
				myfloat U_upper_relative = U_upper - U_element_center;
				myfloat U_lower_relative = U_lower - U_element_center;
				myfloat phi_node = 1.0;
				if (U_node_relative > 0) {
					U_node_relative += U2NITS::Math::EPSILON;
					phi_node = U2NITS::Math::min(1, U_upper_relative / U_node_relative);
				}
				else if (U_node_relative < 0) {
					U_node_relative -= U2NITS::Math::EPSILON;
					phi_node = U2NITS::Math::min(1, U_lower_relative / U_node_relative);
				}
				phi = U2NITS::Math::min(phi, phi_node);
			}

			// 修正梯度
			Ux_element[iVar][iElement] *= phi;
			Uy_element[iVar][iElement] *= phi;
		}
	}
}

void U2NITS::Space::Gradient::Gradient(GPU::ElementSoA& element_host, GPU::NodeSoA& node_host, GPU::EdgeSoA& edge_host, GPU::ElementFieldSoA& elementField_host) {
	/*
	虽然对于最小二乘，节点相对值U_node_relative已经计算，限制器中无需再计算
	但对于业务代码，最好把求梯度和限制器的功能拆开，分成两个循环。
	*/

	static bool is_called_first_time = true;
	// 求梯度
	switch (GlobalPara::inviscid_flux_method::flag_gradient) {
	case _GRA_leastSquare:
		GradientLeastSquare_2(element_host, elementField_host, node_host);
		break;
	case _GRA_greenGauss:
		GradientGreenGauss_2(element_host, elementField_host, edge_host);
		break;

	default:
		LogWriter::logAndPrintError("invalid graident type.\n");
		exit(-1);
		break;
	}
	// 限制器
	switch (GlobalPara::inviscid_flux_method::flux_limiter) {
	case _LIM_none:
		if (is_called_first_time) {
			LogWriter::logAndPrintWarning("No limiter used.\n");
		}
		break;

	case _LIM_barth:
		LimiterBarth(elementField_host, element_host, node_host);
		break;

	default:
		LogWriter::logAndPrintError("invalid limiter type.\n");
		exit(-1);
		break;
	}

	if (is_called_first_time) {
		is_called_first_time = false;
	}
}

void U2NITS::Space::Gradient::GradientLeastSquare_2(GPU::ElementSoA& element_host, GPU::ElementFieldSoA& elementField_host, GPU::NodeSoA& node_host) {

	for (myint iElement = 0; iElement < element_host.num_element; iElement++) {


		// 提取有效邻居
		int validNeighbors[3]{ -1,-1,-1 };// 有效邻居的ID
		int numOfValidNeighbor = 0;
		getElementT3ValidNeighbor(validNeighbors, numOfValidNeighbor, iElement, element_host);
		if (numOfValidNeighbor <= 0) {
			LogWriter::logAndPrintError("invalid numOfValidNeighbor.\n");
			exit(-1);
		}
		// 申请堆数组
		const int nVar = 4;// 守恒量个数，二维为4
		const int nX = 2;// 坐标个数，二维为2
		const int& nVN = numOfValidNeighbor;// 有效邻居个数，取值1-3
		myfloat* dX = new myfloat[nVN * nX]{};
		myfloat* dUdX = new myfloat[nX * nVar]{};
		myfloat* dU = new myfloat[nVN * nVar]{};
		myfloat* dXtrans = new myfloat[nX * nVN]{};
		myfloat* invdXdX = new myfloat[nX * nX]{};
		myfloat* invdXdX_dX = new myfloat[nX * nVN]{};
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
		U2NITS::Math::Matrix::transpose(nVN, nX, (myfloat*)dX, (myfloat*)dXtrans);
		U2NITS::Math::Matrix::mul_ixj_jxk(nX, nVN, nX, (myfloat*)dXtrans, (myfloat*)dX, (myfloat*)invdXdX);
		U2NITS::Math::Matrix::inv_2x2(invdXdX);
		U2NITS::Math::Matrix::mul_ixj_jxk(nX, nX, nVN, (myfloat*)invdXdX, (myfloat*)dXtrans, (myfloat*)invdXdX_dX);
		U2NITS::Math::Matrix::mul_ixj_jxk(nX, nVN, nVar, (myfloat*)invdXdX_dX, (myfloat*)dU, (myfloat*)dUdX);

		/*
		限制器
		一维下，minmod和vanleer限制器的目的是保证两侧界面上物理量的值不超过相邻单元的物理量的值
		参见计算流体力学（李启兵）课件合集P106

		类似地，可以保证三角形顶点处U不超过相邻单元U
		根据当前斜率Ux, Uy，计算顶点相对值dU_node，前面已经算出邻居中心相对值dU
		计算ratio = max(abs(dU_node/dU))，若大于1，则整体除以ratio
		*/
		//const int nNode = 3;
		//myfloat dU_node[nNode][nVar]{};
		//myfloat dX_node[nNode][nX]{};
		//// 计算dX_node
		//for (int iNode = 0; iNode < nNode; iNode++) {
		//	for (int iVar = 0; iVar < nVar; iVar++) {
		//		int nodeID = element_host.nodes[iNode][iElement];
		//		dX_node[iNode][0] = node_host.xy[0][nodeID] - element_host.xy[0][iElement];
		//		dX_node[iNode][1] = node_host.xy[1][nodeID] - element_host.xy[1][iElement];
		//	}
		//}
		//// 计算dU_node[nNode,nVar] =  dX_node[nNode,nX] * dUdX[nX,nVar]，然后计算ratio
		//U2NITS::Math::Matrix::mul_ixj_jxk(nNode, nX, nVar, (myfloat*)dX_node, dUdX, (myfloat*)dU_node);
		//myfloat ratio = 1.0;
		//for (int iNode = 0; iNode < nNode; iNode++) {
		//	for (int iVar = 0; iVar < nVar; iVar++) {
		//		using namespace U2NITS::Math;
		//		if (dU[iNode * nVar + iVar] == 0)continue;// 除零
		//		ratio = max(ratio, abs(dU_node[iNode][iVar] / dU[iNode * nVar + iVar]));
		//	}
		//}
		//// ratio一定大于等于1.0，因此可以直接除以ratio
		//U2NITS::Math::Matrix::div_matrix_by_scalar(nX, nVar, dUdX, ratio);

		// 将dUdX存进单元Ux Uy
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

void U2NITS::Space::Gradient::GradientGreenGauss_2(GPU::ElementSoA& element_host, GPU::ElementFieldSoA& elementField_host, GPU::EdgeSoA& edge_host) {
	/*
	2024-03-28

	* 格林高斯梯度
	* https://zhuanlan.zhihu.com/p/370586072
	* 将面矢量加到单元时，应注意面是朝里还是朝外。
	*
	* 对于边界，其通量应等于phiC(参见01-01-CFD理论)
	* 对于非法边(三角形单元的第4条边)，edgeID=-1, outElementID=-1,
	*/
	// 求梯度的子迭代步数
	int subIterationGradient = 3;
	int num_edge = edge_host.num_edge;

	for (int iElement = 0; iElement < element_host.num_element; iElement++) {
		// 求当前单元ID、当前单元体积、edgeID、edge朝向正负、外单元ID
		int currentElementID = element_host.ID[iElement];// 建议不要直接用iElement，虽然目前iElement==i(element_host.ID[i]=element_i.GPUID=i)，但有隐患
		myfloat volumeC = element_host.volume[currentElementID];// 体积
		// 初始化edgeID 超出数组范围的非法值赋为-1
		int edgeID[4]{ -1,-1,-1,-1 };// 第iElement个element的第i条边的ID
		for (int i = 0; i < 4; i++) {
			edgeID[i] = element_host.edges[i][iElement];
			if (!edge_host.has(edgeID[i])) { // !edge_host.has(edgeID[i])  edgeID[i] <= -1 || edgeID[i] >= num_edge
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
		myfloat edgeNormal[2][4]{};
		myfloat edgeArea[4]{};
		myfloat edgeSVector[2][4]{};// 面积矢量
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
		myfloat gC_all[4]{};
		myfloat rf_rC_all[4][2]{};
		myfloat rf_rF_all[4][2]{};
		for (int i = 0; i < 4; i++) {
			int edgeIDi = edgeID[i];
			if (edgeIDi == -1) {
				continue;
			}

			int currentOutElementID = outElementID[i];
			if (currentOutElementID != -1) {
				// 内部边和周期边
				if (edgeType[i] == -1) {
					// 内部边

					// 面和单元坐标
					myfloat fx = edge_host.xy[0][edgeIDi];// face x
					myfloat fy = edge_host.xy[1][edgeIDi];
					myfloat Cx = element_host.xy[0][currentElementID];// cell x
					myfloat Cy = element_host.xy[1][currentElementID];
					myfloat Fx = element_host.xy[0][currentOutElementID];
					myfloat Fy = element_host.xy[1][currentOutElementID];
					// 几何权重因子gC = |rF-rf|/|rF-rC|，存储于gC_all数组中
					using namespace U2NITS;
					myfloat dFf = Math::distance2D(Fx, Fy, fx, fy);
					myfloat dFC = Math::distance2D(Fx, Fy, Cx, Cy);
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
					myfloat fx = edge_host.xy[0][edgeIDi];
					myfloat fy = edge_host.xy[1][edgeIDi];
					myfloat Cx = element_host.xy[0][currentElementID];// cell x
					myfloat Cy = element_host.xy[1][currentElementID];
					// 周期边的outElement坐标需要平移
					int pairIDi = edge_host.periodicPair[edgeIDi];
					myfloat fx_pair = edge_host.xy[0][pairIDi];
					myfloat fy_pair = edge_host.xy[1][pairIDi];
					myfloat shiftx = fx - fx_pair;// 从pair指向当前位置的矢径
					myfloat shifty = fy - fy_pair;
					myfloat Fx = element_host.xy[0][currentOutElementID];
					myfloat Fy = element_host.xy[1][currentOutElementID];
					Fx += shiftx;// 用该矢径平移单元
					Fy += shifty;

					// 后面和内部边是一样的
					// 几何权重因子gC = |rF-rf|/|rF-rC|，存储于gC_all数组中
					using namespace U2NITS;
					myfloat dFf = Math::distance2D(Fx, Fy, fx, fy);
					myfloat dFC = Math::distance2D(Fx, Fy, Cx, Cy);
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
				myfloat fx = edge_host.xy[0][edgeIDi];// face x
				myfloat fy = edge_host.xy[1][edgeIDi];
				myfloat Cx = element_host.xy[0][currentElementID];// cell x
				myfloat Cy = element_host.xy[1][currentElementID];

				using namespace U2NITS;
				myfloat gC = 0.5;// 几何权重因子gc
				myfloat rf_rC[2]{ fx - Cx,fy - Cy };// rf-rC
				myfloat rf_rF[2]{ -(fx - Cx),-(fy - Cy) };// rf-rF

				gC_all[i] = gC;
				rf_rC_all[i][0] = rf_rC[0];
				rf_rC_all[i][1] = rf_rC[1];
				rf_rF_all[i][0] = rf_rF[0];
				rf_rF_all[i][1] = rf_rF[1];
			}

		}


		// 计算每条边的向量通量phif・Sf，然后求和，除以体积，得到单元C的梯度
		const int nVar = 4;
		for (int jVar = 0; jVar < nVar; jVar++) {
			// 求和。若面朝内，该面通量取相反数
			myfloat sum_phifSf[2]{};
			for (int i = 0; i < 4; i++) {
				int currentEdgeID = edgeID[i];
				if (currentEdgeID == -1) {
					continue;
				}
				int currentOutElementID = outElementID[i];

				if (currentOutElementID != -1) {
					// 内部边和周期边

					// 获取单元值和单元梯度值
					myfloat phiC = elementField_host.U[jVar][currentElementID];
					myfloat phiF = elementField_host.U[jVar][currentOutElementID];
					// 求界面值(近似)
					myfloat gC = gC_all[i];
					myfloat gC1 = 1.0 - gC;
					myfloat phif = gC * phiC + gC1 * phiF;
					myfloat Sf[2]{ edgeSVector[0][i],edgeSVector[1][i] };
					sum_phifSf[0] += phif * Sf[0] * edgeSign[i];// 乘以符号，以应对面朝内的情形
					sum_phifSf[1] += phif * Sf[1] * edgeSign[i];
				}
				else {
					// 边界边(周期边除外)
					myfloat phiC = elementField_host.U[jVar][currentElementID];
					myfloat phif = phiC;
					myfloat Sf[2]{ edgeSVector[0][i],edgeSVector[1][i] };
					sum_phifSf[0] += phif * Sf[0] * edgeSign[i];// 乘以符号，以应对面朝内的情形
					sum_phifSf[1] += phif * Sf[1] * edgeSign[i];
				}

			}

			myfloat nabla_phiC_new[2]{// 梯度
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


