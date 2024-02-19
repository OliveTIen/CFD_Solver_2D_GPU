#include "GPUSolver.h"
#include "../FVM_2D.h"
#include "space/CalculateGradient.h"
#include "space/CalculateFlux.h"
#include "datatype/EdgeSoA.h"

void GPUSolver::initialze() {
	/*
	��FVM_2D����ת��ΪElement�ṹ����ʼ��
	��ʼ���������У�
	    ��Ԫ���ꡢU����Ԫ�ھ����ꡢ��Ԫ�ھ�U
	��ʱδ��ʼ����
		��Ԫ�ڵ����ꡢ��

	��ȡ�ļ��󣬽�FVM_2D�����ݽṹת��ΪElement�ṹ��
	1. ��ʼ��selfxy, node1, node2, node3 �������������⣬self����CPU�ϼ���ģ��Ժ��ת�Ƶ�GPU�ϣ�
	2. ��ʼ��selfValue.U1~U4, neighbor1.U1~U4, neighbor2.U1~U4, neighbor3.U1~U4
	3. �ݶȳ�ʼ��Ϊ0


	cudaMallocӦ����ѭ���⣬�Լ��ٿ��� �μ���P364
	���Ǽ��㵥Ԫ�ͼ���������ײ�ͬ��С�����飬�����Ҫ��������
	�������Ƶ��豸

	��Ҫ��ʼ����
		��Ԫ����
		��Ԫ�ڵ�����

	*/
	// ����GPU
	int iDeviceCount = 0;
	cudaError_t error = cudaGetDeviceCount(&iDeviceCount);
	if (error != cudaSuccess || iDeviceCount <= 0) {
		printf("No GPU found\n");
		throw error;
	}
	printf("Num of GPU: %d\n", iDeviceCount);

	int iDev = 0;
	error = cudaSetDevice(iDev);
	if (error != cudaSuccess) {
		printf("Fail to set GPU %d\n", iDev);
		throw error;
	}
	printf("Activate GPU: %d\n", iDev);

	// �����ڴ�
	FVM_2D* pFVM2D = FVM_2D::pFVM2D;
	const int num_element = (int) pFVM2D->elements.size();
	const int num_neighbor = 3;
	const int num_edge = pFVM2D->edges.size();
	this->element_host.alloc(num_element);
	this->elementAdjacent.alloc(num_neighbor, num_element);
	this->edge_host.alloc(num_edge);

	// ����device�ڴ�
	try{ 
		element_device.cuda_alloc(num_element); 
		edge_device.cuda_alloc(num_edge);
	}
	catch (const char* e) {
		// ! �쳣����δ���
		std::cout << e << std::endl;
		cudaError_t error = cudaError_t::cudaErrorDeviceUninitialized;
		throw error;
	}

	// ��ʼ��Element_2D vector�е�GPUindex
	#pragma omp parallel for
	for (int i = 0; i < num_element; i++) {
		pFVM2D->elements[i].GPUindex = i;
	}
	#pragma omp barrier
	#pragma omp parallel for
	for (int i = 0; i < num_edge; i++) {
		pFVM2D->edges[i].GPUID = i;
	}
	#pragma omp barrier
	
	// ��ʼ��element_host��elementAdjacent
	initialize_elementHost_elementAdjacent(pFVM2D, num_element, num_neighbor);
	// ��ʼ��edge_host
	initialize_edgeHost(pFVM2D, num_edge);

	// ����host���ݵ�device
	element_host.cuda_copy_to_device(&element_device);
	GPU::EdgeSoA::cuda_memcpy(&edge_device, &edge_host, cudaMemcpyHostToDevice);
}

void GPUSolver::iteration() {
	/*
    
	1. ��ʼ��
		��Ԫ��ֵͨ����Ϊ0
		���ݵ�ԪU��ʼ����Ԫ�ھ�U

    2. ���㵥Ԫ�ݶ� - GPU
	    2.1
        �������ꡢU�������ݶ�
        �����������ݶ�
        ����쳣ֵ
        ���룺��Ԫ���ꡢ��ԪU����Ԫ�ھ����ꡢ��Ԫ�ھ�U
        �������Ԫ�ݶ�

	    2.2. ��Ԫ�ݶȸ��Ƶ�CPU��Ȼ������ھ��ݶ�
		
    4. ����߽���ֵͨ��
        ����߷�����nx ny��������edgex, edgey���߳�edgeLength
		    ��Ҫ�������ڵ�����
        ��������ҵ�Ԫ�ڱ߽��Uֵ(���Ҽ���UL��UR)
        ���ݱ����ͣ�����UL��UR
        ������������
            ��Ҫ�Ĳ�����UL UR nx ny edgeLength
            ������߽����ֵͨ��
    5. ���㵥Ԫ��ֵͨ��
        ���߽���ֵͨ���ӵ���Ԫ��ֵͨ���� - ��Լ������OpenMP
		(ע��ÿ��ѭ����Ҫ����Ԫ��ֵͨ����ʼ��Ϊ0)
    6. ��ʽʱ���ƽ�

 
	*/


	// --- ���㵥Ԫ�ݶ� --- 
	// ���룺��Ԫ���ꡢ��ԪU����Ԫ�ھ����ꡢ��Ԫ�ھ�U
	// �������Ԫ�ݶ�
	GPU::calculateGradient(this->element_device);

	// �����ھ��ݶ� �˴�ֱ�Ӷ�device���� ������P312 ���Գ���������ӳ��
	GPU::updateNeighborGradient(this->element_device, this->elementAdjacent);

	/*
	// ����device���ݵ�host
	element_device.cuda_copy_to_host(&element_host);
	// ������ Ŀǰ����Ҫ��
	
	// �����ھ��ݶ� 1.�����elementAdjacent���飬�洢�ھ�index��Ϣ
	GPU::updateNeighborGradient(this->element_host, this->elementAdjacent);

	// �����ݶ����ݵ�device
	for (int i = 0; i < 3; i++) {
		cudaMemcpy(element_device.neighbors[i].Ux1, element_host.neighbors[i].Ux1, sizeof(REAL) * element_device.num_element, cudaMemcpyHostToDevice);
		cudaMemcpy(element_device.neighbors[i].Ux2, element_host.neighbors[i].Ux2, sizeof(REAL) * element_device.num_element, cudaMemcpyHostToDevice);
		cudaMemcpy(element_device.neighbors[i].Ux3, element_host.neighbors[i].Ux3, sizeof(REAL) * element_device.num_element, cudaMemcpyHostToDevice);
		cudaMemcpy(element_device.neighbors[i].Ux4, element_host.neighbors[i].Ux4, sizeof(REAL) * element_device.num_element, cudaMemcpyHostToDevice);
		cudaMemcpy(element_device.neighbors[i].Uy1, element_host.neighbors[i].Uy1, sizeof(REAL) * element_device.num_element, cudaMemcpyHostToDevice);
		cudaMemcpy(element_device.neighbors[i].Uy2, element_host.neighbors[i].Uy2, sizeof(REAL) * element_device.num_element, cudaMemcpyHostToDevice);
		cudaMemcpy(element_device.neighbors[i].Uy3, element_host.neighbors[i].Uy3, sizeof(REAL) * element_device.num_element, cudaMemcpyHostToDevice);
		cudaMemcpy(element_device.neighbors[i].Uy4, element_host.neighbors[i].Uy4, sizeof(REAL) * element_device.num_element, cudaMemcpyHostToDevice);
	}
	*/


	// --- ����߽���ֵͨ�� --- 
	// ���룺��Ԫ
	// ������߽���ֵͨ��
	GPU::calculateFlux(this->element_device);
}

// inline����ֻ����ͷ�ļ��ж��壬����ᱨ��LNK2001
void GPUSolver::finalize() {
	this->element_host.free(); 
	this->element_device.cuda_free();
	this->elementAdjacent.free();
}

void GPUSolver::initialize_elementHost_elementAdjacent(void* _pFVM2D_, int num_element, int num_neighbor) {
	// ��ʼ��element_host��elementAdjacent
	FVM_2D* pFVM2D = (FVM_2D*)_pFVM2D_;
#pragma omp parallel for
	for (int i = 0; i < num_element; i++) {

		Element_2D& element_i = pFVM2D->elements[i];

		// ��ʼ����Ԫ���ꡢ��ԪU����Ԫ�ݶ�
		element_host.self.isNull[i] = false;
		element_host.self.x[i] = element_i.x;
		element_host.self.y[i] = element_i.y;
		element_host.self.U1[i] = element_i.U[0];
		element_host.self.U2[i] = element_i.U[1];
		element_host.self.U3[i] = element_i.U[2];
		element_host.self.U4[i] = element_i.U[3];
		element_host.self.Ux1[i] = 0;
		element_host.self.Uy1[i] = 0;
		element_host.self.Ux2[i] = 0;
		element_host.self.Uy2[i] = 0;
		element_host.self.Ux3[i] = 0;
		element_host.self.Uy3[i] = 0;
		element_host.self.Ux4[i] = 0;
		element_host.self.Uy4[i] = 0;

		// ��ʼ���ھ����ꡢ�ھ�U���ھ��ݶ�
		std::vector<Element_2D*> neighbors_element_i = element_i.findNeighbor();// �ھ�ָ��
		// �ж��ھ��Ƿ�Ϊnullptr������ֵ���˴�Ĭ��Ϊ3���ھ� ! ����������
		for (int j = 0; j < num_neighbor; j++) {
			if (neighbors_element_i[j] == nullptr) {
				element_host.neighbors[j].isNull[i] = true;
				elementAdjacent.neighbors[j].index[i] = -1;
			}
			else {
				element_host.neighbors[j].isNull[i] = false;
				element_host.neighbors[j].x[i] = neighbors_element_i[j]->x;
				element_host.neighbors[j].y[i] = neighbors_element_i[j]->y;
				element_host.neighbors[j].U1[i] = neighbors_element_i[j]->U[0];
				element_host.neighbors[j].U2[i] = neighbors_element_i[j]->U[1];
				element_host.neighbors[j].U3[i] = neighbors_element_i[j]->U[2];
				element_host.neighbors[j].U4[i] = neighbors_element_i[j]->U[3];

				element_host.neighbors[j].Ux1[i] = 0;
				element_host.neighbors[j].Uy1[i] = 0;
				element_host.neighbors[j].Ux2[i] = 0;
				element_host.neighbors[j].Uy2[i] = 0;
				element_host.neighbors[j].Ux3[i] = 0;
				element_host.neighbors[j].Uy3[i] = 0;
				element_host.neighbors[j].Ux4[i] = 0;
				element_host.neighbors[j].Uy4[i] = 0;

				// ��ʼ��elementAdjacent
				elementAdjacent.neighbors[j].index[i] = neighbors_element_i[j]->GPUindex;
			}


		}
	}

}

void GPUSolver::initialize_edgeHost(void* _pFVM2D_, int num_edge) {
	FVM_2D* pFVM2D = (FVM_2D*)_pFVM2D_;
#pragma omp parallel for
	for (int i = 0; i < num_edge; i++) {
		Edge_2D& edge_i = pFVM2D->edges[i];
		edge_host.ID[i]=edge_i.GPUID;
		edge_host.nodes[0][i] = edge_i.nodes[0];
		edge_host.nodes[1][i] = edge_i.nodes[1];
		edge_host.setID[i] = edge_i.setID;
		edge_host.elementL[i] = edge_i.pElement_L->GPUindex;
		edge_host.elementR[i] = edge_i.pElement_R->GPUindex;
		edge_host.length[i] = edge_i.length;
		edge_host.distanceOfElements[i] = edge_i.refLength;
	}

	/*
		int* ID;
		int* nodes[2];
		int* setID;
		int* elementL;
		int* elementR;
		REAL* length;
		REAL* distanceOfElements;// ���൥Ԫ���ľ���
	
	*/
}
