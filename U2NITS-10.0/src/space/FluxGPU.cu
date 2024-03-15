#include "FluxGPU.h"
#include "../FVM_2D.h"
#include "../math/PhysicalConvertKernel.h"

void GPU::Space::Flux::calculateFluxDevice(ElementSoA& element_device, FieldSoA& elementField_device, EdgeSoA& edge_device, BoundarySetMap& boundary_device, DInfPara& infPara_device) {
	// TODO: Implement this function.
    resetElementFlux(elementField_device);
    cudaDeviceSynchronize();
    calculateFlux(element_device, elementField_device, edge_device, boundary_device, infPara_device);


}

void GPU::Space::Flux::resetElementFlux(FieldSoA& elementField_device) {
    // ��ʼ����Ԫ��ֵͨ��
    int block_size = 512;// �����128 256 512
    int grid_size = (elementField_device.num + block_size - 1) / block_size;
    dim3 block(block_size, 1, 1);
    dim3 grid(grid_size, 1, 1);
    resetElementFluxKernel <<<grid, block>>> (elementField_device);
}

__global__ void GPU::Space::Flux::resetElementFluxKernel(FieldSoA& elementField_device) {
    // --- �˺��� --- 

    // ��ȡid���ж��Ƿ���Ч
    const int bid = blockIdx.x;
    const int tid = threadIdx.x;
    const int id = bid * blockDim.x + tid;
    const int& i = id;
    if (i >= elementField_device.num || i < 0) return;

    // ����Ԫ��ֵͨ����ֵ
    for (int jValue = 0; jValue < 4; jValue++) {
        elementField_device.Flux[jValue][i] = 0;//��Ԫ��ֵͨ��Flux���㣬Ϊ����Ӽ���׼��
        // f->elements[ie].deltaeig = 0;//ÿһ��deltaeig���� ���Roe��ʽ Ŀǰ����Ҫ
        // �ݶ��Ѿ����㣬��˲���Ҫ�ٴμ���
    }
}

void GPU::Space::Flux::calculateFlux(ElementSoA& element_device, FieldSoA& elementField_device, EdgeSoA& edge_device, BoundarySetMap& boundary_device, DInfPara& infPara_device) {

    int block_size = 512;
    int grid_size = (edge_device.num_edge + block_size - 1) / block_size;
    dim3 block(block_size, 1, 1);
    dim3 grid(grid_size, 1, 1);
    calculateFluxKernel <<<grid, block >>> (element_device, elementField_device, edge_device, boundary_device, infPara_device);
}

__global__ void GPU::Space::Flux::calculateFluxKernel(ElementSoA& element_device, FieldSoA& elementField_device, EdgeSoA& edge_device, BoundarySetMap& boundary_device, DInfPara& infPara_device) {
    // TODO: Implement this kernel.
    // δ���

    // ���ڱ߽磬EdgeSoA��periodicPair��Ա��������ɣ��Ѿ�������edge_device�У������ٴ���
    // ��ȡ�߽�����ʱ������ʹ��boundary_host

    // ��ȡid���ж��Ƿ���Ч
    const int bid = blockIdx.x;
    const int tid = threadIdx.x;
    const int id = bid * blockDim.x + tid;
    const int& iedge = id;
    if (iedge >= edge_device.num_edge || iedge < 0) return;

    // ÿ���߼�����ճ��ֵͨ��

    //FVM_2D* f = FVM_2D::getInstance();
    int elementL = edge_device.elementL[iedge];
    int elementR = edge_device.elementR[iedge];
    int setID = edge_device.setID[iedge];// �ڲ�edge�������κ�set�����setIDΪ-1 setID�ĳ�ʼ����readContinueFile
    //int boundaryNum = f->boundaryManager.boundaries.size();
    int boundaryNum = boundary_device.size;
    int bType = -1;
    if (setID != -1) {
        bType = boundary_device.type[setID - 1];
    }
    double flux[4]{};

    switch (bType) {
        //�Գơ���ŷ�����̣��൱����ճ�̱�
    case _BC_symmetry:
        GPU::Space::Flux::getEdgeFlux_wallNonViscous(element_device, elementField_device, edge_device, iedge, flux);
        break;
        //��ճ�̱�
    case _BC_wall_nonViscous:
        GPU::Space::Flux::getEdgeFlux_wallNonViscous(element_device, elementField_device, edge_device, iedge, flux);
        break;
        //�޻��ƾ���
    case _BC_wall_adiabat:
        // TODO: ʵ���޻��ƾ��ȱ߽������ļ���
        GPU::Space::Flux::getEdgeFlux_wall_adiabat(element_device, elementField_device, edge_device, iedge, flux);
        break;
        //���
    case _BC_inlet:
        GPU::Space::Flux::getEdgeFlux_farField(element_device, elementField_device, edge_device, iedge,
            infPara_device.ruvp_inlet->ptr, flux);
        break;
        //����
    case _BC_outlet:
        GPU::Space::Flux::getEdgeFlux_farField(element_device, elementField_device, edge_device, iedge,
            infPara_device.ruvp_outlet->ptr, flux);
        break;
        //Զ��
    case _BC_inf:
        GPU::Space::Flux::getEdgeFlux_farField(element_device, elementField_device, edge_device, iedge,
            infPara_device.ruvp_inf->ptr, flux);
        break;
    default:// �ڲ���bType=-1
        if (elementR != -1) {
            // ���ں��ڲ� ͳһ����
            GPU::Space::Flux::getEdgeFlux_inner_and_periodic(element_device, elementField_device, edge_device, iedge, flux);
        }
    }

    // ���µ�Ԫ��ֵͨ��
    for (int j = 0; j < 4; j++) {
        elementField_device.Flux[j][elementL] += flux[j];
        // �ڲ��߽�
        if (bType == -1) {
            // �˴�����elementR==-1����ֹ���ڱ߽���2��
            elementField_device.Flux[j][elementR] -= flux[j];
        }
    }
}

__device__ void GPU::Space::Flux::getEdgeFlux_wallNonViscous(ElementSoA& element_host, FieldSoA& elementField_host, EdgeSoA& edge_host, long idx, double* flux) {
    δ�����;
}
__device__ void GPU::Space::Flux::getEdgeFlux_wall_adiabat(ElementSoA& element_host, FieldSoA& elementField_host, EdgeSoA& edge_host, long idx, double* flux){

}
__device__ void GPU::Space::Flux::getEdgeFlux_farField(ElementSoA& element_host, FieldSoA& elementField_host, EdgeSoA& edge_host, long idx, REAL* ruvp_inf, REAL* flux){

}
__device__ void GPU::Space::Flux::modify_ruvpL_farField(const REAL nx, const REAL ny, REAL* ruvp, const REAL* ruvp_inf){}
__device__ void GPU::Space::Flux::getEdgeFlux_inner_and_periodic(ElementSoA& element_host, FieldSoA& elementField_host, EdgeSoA& edge_host, long idx, REAL* flux){}
__device__ void GPU::Space::Flux::getUByXYandElementID(ElementSoA& element_host, FieldSoA& elementField_host, REAL x, REAL y, int elementID, REAL* U) {}

__device__ void GPU::Space::Flux::RiemannSolve(const REAL* UL, const REAL* UR, const REAL nx, const REAL ny, const REAL length, REAL* flux, const int scheme) {

}