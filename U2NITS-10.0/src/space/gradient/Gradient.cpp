#include "Gradient.h"
#include "../../math/Math.h"
#include "../../output/LogWriter.h"
#include "../../global/GlobalPara.h"

// ��������
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
	���ݶ�����ת��Ϊ2x4������ʽ
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
	��ȡ��Ч�ھ�
	����Ͳ�����validNeighbors��Ч�ھӣ�numOfValidNeighbor��Ч�ھӸ���
	*/
	numOfValidNeighbor = 0;
	for (int j = 0; j < 3; j++) {// �����ε�Ԫ�����3���ھ�
		if (element_host.neighbors[j][iElement] != -1) {
			validNeighbors[numOfValidNeighbor] = element_host.neighbors[j][iElement];
			numOfValidNeighbor += 1;
		}
	}
}


inline void get_U_linear(myfloat x, myfloat y, myfloat& U_dist, myint i_e, const myfloat* x_e, const myfloat* y_e, const myfloat* U_e, const myfloat* Ux_e, const myfloat* Uy_e) {
	/*
	i_e: ��Ԫindex
	����Ͳ�����U_dist
	*/
	U_dist = U_e[i_e] + Ux_e[i_e] * (x - x_e[i_e]) + Uy_e[i_e] * (y - y_e[i_e]);
}

void LimiterBarth(GPU::ElementFieldSoA& elementField_host, GPU::ElementSoA& element_host, GPU::NodeSoA& node_host) {
	/*
	����Ŀ�ģ������ݶ�
	���ڶ�ά�ǽṹ����Ҫ��֤�����ζ��㴦U���������ڵ�ԪU
	��ʽ�μ����� The design and application of upwind schemes on unstructured meshes

	OpenMP
	����8���̺߳�ȷʵCPUռ�ôﵽ56%��Ȼ�������ٶȷ����½��ˣ���138.036��Ϊ124.603step/s
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
			// ��ȡ�����ھӵ�Ԫ����ֵ���½�U_lower���Ͻ�U_upper
			myfloat U_element_center = U_element[iVar][iElement];
			myfloat U_lower = U_element_center;
			myfloat U_upper = U_element_center;
			for (int j = 0; j < numOfValidNeighbor; j++) {
				myint iNeighbor = validNeighbors[j];
				U_lower = U2NITS::Math::min(U_lower, U_element[iVar][iNeighbor]);
				U_upper = U2NITS::Math::max(U_upper, U_element[iVar][iNeighbor]);
			}

			// ����������ϵ��phi��ȡ����������ϵ��phi_node����Сֵ
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

			// �����ݶ�
			Ux_element[iVar][iElement] *= phi;
			Uy_element[iVar][iElement] *= phi;
		}
	}
}


void U2NITS::Space::Gradient::Gradient(GPU::ElementSoA& element_host, GPU::NodeSoA& node_host, GPU::EdgeSoA& edge_host, GPU::ElementFieldSoA& elementField_host) {
	/*
	��Ȼ������С���ˣ��ڵ����ֵU_node_relative�Ѿ����㣬�������������ټ���
	������ҵ����룬��ð����ݶȺ��������Ĺ��ܲ𿪣��ֳ�����ѭ����
	*/

	/*
	�������ڵľ�̬�������������������������������Ǻ����塣���ڵ�һ�ε��ú���ʱ����ʼ����������ú���ʱ��ʹ��֮ǰ��ֵ��
	�����ڱ��溯�����״̬���μ� https://www.geeksforgeeks.org/static-keyword-cpp/
	*/
	static bool is_called_first_time = true;
	
	// ���ݶ�
	switch (GlobalPara::inviscid_flux_method::flag_gradient) {
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
	// ������
	switch (GlobalPara::inviscid_flux_method::flux_limiter) {
	case _LIM_none:
		if (is_called_first_time) {
			LogWriter::logAndPrintWarning("No limiter used.\n");
		}
		break;

	case _LIM_barth:
		LimiterBarth(elementField_host, element_host, node_host);
		break;

	//case _LIM_minmod:

	//	break;

	//case _LIM_vanleer:

	//	break;

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


		// ��ȡ��Ч�ھ�
		int validNeighbors[3]{ -1,-1,-1 };// ��Ч�ھӵ�ID
		int numOfValidNeighbor = 0;
		getElementT3ValidNeighbor(validNeighbors, numOfValidNeighbor, iElement, element_host);
		if (numOfValidNeighbor <= 0) {
			LogWriter::logAndPrintError("invalid numOfValidNeighbor.\n");
			exit(-1);
		}
		// ���������
		const int nVar = 4;// �غ�����������άΪ4
		const int nX = 2;// �����������άΪ2
		const int& nVN = numOfValidNeighbor;// ��Ч�ھӸ�����ȡֵ1-3
		myfloat* dX = new myfloat[nVN * nX]{};
		myfloat* dUdX = new myfloat[nX * nVar]{};
		myfloat* dU = new myfloat[nVN * nVar]{};
		myfloat* dXtrans = new myfloat[nX * nVN]{};
		myfloat* invdXdX = new myfloat[nX * nX]{};
		myfloat* invdXdX_dX = new myfloat[nX * nVN]{};
		// ��ʼ��dX��dU
		for (int iVN = 0; iVN < nVN; iVN++) {
			int index = validNeighbors[iVN];// �ھӵ�ID
			dX[iVN * nX + 0] = element_host.xy[0][index] - element_host.xy[0][iElement];
			dX[iVN * nX + 1] = element_host.xy[1][index] - element_host.xy[1][iElement];
			for (int iVar = 0; iVar < nVar; iVar++) {
				dU[iVN * nVar + iVar] = elementField_host.U[iVar][index] - elementField_host.U[iVar][iElement];
			}
		}
		// ��С���˷�����dUdX��x=(A'A)^{-1}A'b��dUdX = inv(dXtrans * dX) * dXtrans * dU
		U2NITS::Math::Matrix::transpose(nVN, nX, (myfloat*)dX, (myfloat*)dXtrans);
		U2NITS::Math::Matrix::mul_ixj_jxk(nX, nVN, nX, (myfloat*)dXtrans, (myfloat*)dX, (myfloat*)invdXdX);
		U2NITS::Math::Matrix::inv_2x2(invdXdX);
		U2NITS::Math::Matrix::mul_ixj_jxk(nX, nX, nVN, (myfloat*)invdXdX, (myfloat*)dXtrans, (myfloat*)invdXdX_dX);
		U2NITS::Math::Matrix::mul_ixj_jxk(nX, nVN, nVar, (myfloat*)invdXdX_dX, (myfloat*)dU, (myfloat*)dUdX);

		/*
		������
		һά�£�minmod��vanleer��������Ŀ���Ǳ�֤�����������������ֵ���������ڵ�Ԫ����������ֵ
		�μ�����������ѧ�����������μ��ϼ�P106

		���Ƶأ����Ա�֤�����ζ��㴦U���������ڵ�ԪU
		���ݵ�ǰб��Ux, Uy�����㶥�����ֵdU_node��ǰ���Ѿ�����ھ��������ֵdU
		����ratio = max(abs(dU_node/dU))��������1�����������ratio
		*/
		//const int nNode = 3;
		//myfloat dU_node[nNode][nVar]{};
		//myfloat dX_node[nNode][nX]{};
		//// ����dX_node
		//for (int iNode = 0; iNode < nNode; iNode++) {
		//	for (int iVar = 0; iVar < nVar; iVar++) {
		//		int nodeID = element_host.nodes[iNode][iElement];
		//		dX_node[iNode][0] = node_host.xy[0][nodeID] - element_host.xy[0][iElement];
		//		dX_node[iNode][1] = node_host.xy[1][nodeID] - element_host.xy[1][iElement];
		//	}
		//}
		//// ����dU_node[nNode,nVar] =  dX_node[nNode,nX] * dUdX[nX,nVar]��Ȼ�����ratio
		//U2NITS::Math::Matrix::mul_ixj_jxk(nNode, nX, nVar, (myfloat*)dX_node, dUdX, (myfloat*)dU_node);
		//myfloat ratio = 1.0;
		//for (int iNode = 0; iNode < nNode; iNode++) {
		//	for (int iVar = 0; iVar < nVar; iVar++) {
		//		using namespace U2NITS::Math;
		//		if (dU[iNode * nVar + iVar] == 0)continue;// ����
		//		ratio = max(ratio, abs(dU_node[iNode][iVar] / dU[iNode * nVar + iVar]));
		//	}
		//}
		//// ratioһ�����ڵ���1.0����˿���ֱ�ӳ���ratio
		//U2NITS::Math::Matrix::div_matrix_by_scalar(nX, nVar, dUdX, ratio);

		// ��dUdX�����ԪUx Uy
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

	* ���ָ�˹�ݶ�
	* https://zhuanlan.zhihu.com/p/370586072
	* ����ʸ���ӵ���Ԫʱ��Ӧע�����ǳ��ﻹ�ǳ��⡣
	*
	* ���ڱ߽磬��ͨ��Ӧ����phiC(�μ�01-01-CFD����)
	* ���ڷǷ���(�����ε�Ԫ�ĵ�4����)��edgeID=-1, outElementID=-1,
	*/
	// ���ݶȵ��ӵ�������
	int subIterationGradient = 3;
	int num_edge = edge_host.num_edge;

	for (int iElement = 0; iElement < element_host.num_element; iElement++) {
		// ��ǰ��ԪID����ǰ��Ԫ�����edgeID��edge�����������ⵥԪID
		int currentElementID = element_host.ID[iElement];// ���鲻Ҫֱ����iElement����ȻĿǰiElement==i(element_host.ID[i]=element_i.GPUID=i)����������
		myfloat volumeC = element_host.volume[currentElementID];// ���
		// ��ʼ��edgeID �������鷶Χ�ķǷ�ֵ��Ϊ-1
		int edgeID[4]{ -1,-1,-1,-1 };// ��iElement��element�ĵ�i���ߵ�ID
		for (int i = 0; i < 4; i++) {
			edgeID[i] = element_host.edges[i][iElement];
			if (!edge_host.has(edgeID[i])) { // !edge_host.has(edgeID[i])  edgeID[i] <= -1 || edgeID[i] >= num_edge
				edgeID[i] = -1;
			}
		}
		// ��ʼ��edgeType
		int edgeType[4]{ -1,-1,-1,-1 };
		for (int i = 0; i < 4; i++) {
			int edgeIDi = edgeID[i];
			if (edgeIDi == -1) {
				continue;
			}

			edgeType[i] = edge_host.setID[edgeIDi];
		}
		// ��ʼ��edgeSign��outElementID
		int edgeSign[4]{ 0,0,0,0 };// 1��ʾ�߳��⣬-1��ʾ�߳���
		int outElementID[4]{ -1,-1,-1,-1 };//�ⵥԪID �ڲ��ߺ����ڱ�
		for (int i = 0; i < 4; i++) {
			int edgeIDi = edgeID[i];
			if (edgeIDi == -1) {
				continue;
			}

			if (currentElementID != edge_host.elementR[edgeIDi]) {
				edgeSign[i] = 1;// currentElement=elementL������
				outElementID[i] = edge_host.elementR[edgeIDi];
			}
			else {
				edgeSign[i] = -1;// currentElement=elementR������
				outElementID[i] = edge_host.elementL[edgeIDi];
			}
		}
		// �����ʸ��edgeSVector(������)
		myfloat edgeNormal[2][4]{};
		myfloat edgeArea[4]{};
		myfloat edgeSVector[2][4]{};// ���ʸ��
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

		// �󼸺�Ȩ������gc��ʸ��(������)
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
				// �ڲ��ߺ����ڱ�
				if (edgeType[i] == -1) {
					// �ڲ���

					// ��͵�Ԫ����
					myfloat fx = edge_host.xy[0][edgeIDi];// face x
					myfloat fy = edge_host.xy[1][edgeIDi];
					myfloat Cx = element_host.xy[0][currentElementID];// cell x
					myfloat Cy = element_host.xy[1][currentElementID];
					myfloat Fx = element_host.xy[0][currentOutElementID];
					myfloat Fy = element_host.xy[1][currentOutElementID];
					// ����Ȩ������gC = |rF-rf|/|rF-rC|���洢��gC_all������
					using namespace U2NITS;
					myfloat dFf = Math::distance2D(Fx, Fy, fx, fy);
					myfloat dFC = Math::distance2D(Fx, Fy, Cx, Cy);
					if (dFC < Math::EPSILON)dFC = Math::EPSILON;
					gC_all[i] = dFf / dFC;
					// ʸ�� rf-rC rf-rF���ڵ���ʱ���õ�
					rf_rC_all[i][0] = fx - Cx;// rf-rC
					rf_rC_all[i][1] = fy - Cy;
					rf_rF_all[i][0] = fx - Fx;// rf-rF
					rf_rF_all[i][1] = fy - Fy;
				}
				else {
					// ���ڱ�
					// ��͵�Ԫ����
					myfloat fx = edge_host.xy[0][edgeIDi];
					myfloat fy = edge_host.xy[1][edgeIDi];
					myfloat Cx = element_host.xy[0][currentElementID];// cell x
					myfloat Cy = element_host.xy[1][currentElementID];
					// ���ڱߵ�outElement������Ҫƽ��
					int pairIDi = edge_host.periodicPair[edgeIDi];
					myfloat fx_pair = edge_host.xy[0][pairIDi];
					myfloat fy_pair = edge_host.xy[1][pairIDi];
					myfloat shiftx = fx - fx_pair;// ��pairָ��ǰλ�õ�ʸ��
					myfloat shifty = fy - fy_pair;
					myfloat Fx = element_host.xy[0][currentOutElementID];
					myfloat Fy = element_host.xy[1][currentOutElementID];
					Fx += shiftx;// �ø�ʸ��ƽ�Ƶ�Ԫ
					Fy += shifty;

					// ������ڲ�����һ����
					// ����Ȩ������gC = |rF-rf|/|rF-rC|���洢��gC_all������
					using namespace U2NITS;
					myfloat dFf = Math::distance2D(Fx, Fy, fx, fy);
					myfloat dFC = Math::distance2D(Fx, Fy, Cx, Cy);
					if (dFC < Math::EPSILON)dFC = Math::EPSILON;
					gC_all[i] = dFf / dFC;
					// ʸ�� rf-rC rf-rF���ڵ���ʱ���õ�
					rf_rC_all[i][0] = fx - Cx;// rf-rC
					rf_rC_all[i][1] = fy - Cy;
					rf_rF_all[i][0] = fx - Fx;// rf-rF
					rf_rF_all[i][1] = fy - Fy;
				}

			}
			else {
				// �߽��(���ڱ߽����)
				// F��Ԫʵ�ʲ����ڣ�����F��Ԫ��C��Ԫ����f�����ĶԳ�
				myfloat fx = edge_host.xy[0][edgeIDi];// face x
				myfloat fy = edge_host.xy[1][edgeIDi];
				myfloat Cx = element_host.xy[0][currentElementID];// cell x
				myfloat Cy = element_host.xy[1][currentElementID];

				using namespace U2NITS;
				myfloat gC = 0.5;// ����Ȩ������gc
				myfloat rf_rC[2]{ fx - Cx,fy - Cy };// rf-rC
				myfloat rf_rF[2]{ -(fx - Cx),-(fy - Cy) };// rf-rF

				gC_all[i] = gC;
				rf_rC_all[i][0] = rf_rC[0];
				rf_rC_all[i][1] = rf_rC[1];
				rf_rF_all[i][0] = rf_rF[0];
				rf_rF_all[i][1] = rf_rF[1];
			}

		}



		// ����ÿ���ߵ�����ͨ��phif��Sf��Ȼ����ͣ�����������õ���ԪC���ݶ�
		const int nVar = 4;
		for (int jVar = 0; jVar < nVar; jVar++) {
			// ��͡����泯�ڣ�����ͨ��ȡ�෴��
			myfloat sum_phifSf[2]{};
			for (int i = 0; i < 4; i++) {
				int currentEdgeID = edgeID[i];
				if (currentEdgeID == -1) {
					continue;
				}
				int currentOutElementID = outElementID[i];

				if (currentOutElementID != -1) {
					// �ڲ��ߺ����ڱ�

					// ��ȡ��Ԫֵ�͵�Ԫ�ݶ�ֵ
					myfloat phiC = elementField_host.U[jVar][currentElementID];
					myfloat phiF = elementField_host.U[jVar][currentOutElementID];
					// �����ֵ(����)
					myfloat gC = gC_all[i];
					myfloat gC1 = 1.0 - gC;
					myfloat phif = gC * phiC + gC1 * phiF;
					myfloat Sf[2]{ edgeSVector[0][i],edgeSVector[1][i] };
					sum_phifSf[0] += phif * Sf[0] * edgeSign[i];// ���Է��ţ���Ӧ���泯�ڵ�����
					sum_phifSf[1] += phif * Sf[1] * edgeSign[i];
				}
				else {
					// �߽��(���ڱ߳���)
					myfloat phiC = elementField_host.U[jVar][currentElementID];
					myfloat phif = phiC;
					myfloat Sf[2]{ edgeSVector[0][i],edgeSVector[1][i] };
					sum_phifSf[0] += phif * Sf[0] * edgeSign[i];// ���Է��ţ���Ӧ���泯�ڵ�����
					sum_phifSf[1] += phif * Sf[1] * edgeSign[i];
				}

			}

			myfloat nabla_phiC_new[2]{// �ݶ�
				1 / volumeC * sum_phifSf[0],
				1 / volumeC * sum_phifSf[1]
			};
			// ���µ�Ԫ�ݶ�ֵ
			elementField_host.Ux[jVar][currentElementID] = nabla_phiC_new[0];
			elementField_host.Uy[jVar][currentElementID] = nabla_phiC_new[1];
			// ���Ҫ��������Ҫ�����е�Ԫ������ϣ���˵���ѭ��Ҫ����iElementѭ����
		}
	}
}


