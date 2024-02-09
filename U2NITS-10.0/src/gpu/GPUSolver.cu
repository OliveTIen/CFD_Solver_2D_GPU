#include "GPUSolver.h"
#include "../FVM_2D.h"
#include "GPU_space.h"


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
	*/
	/*
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
	int num_element = pFVM2D->elements.size();
	this->element_host.alloc(num_element);
	
	try{ 
		element_device.cuda_alloc(num_element); 
	}
	catch (const char* e) {
		// ! �쳣����δ���
		std::cout << e << std::endl;
		cudaError_t error = cudaError_t::cudaErrorDeviceUninitialized;
		throw error;
	}
	
	// ��ʼ��host����
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
		const int num_neighbor = 3;
		std::vector<Element_2D*> neighbors_element_i = element_i.findNeighbor();// �ж��ھ��Ƿ�Ϊnullptr������ֵ���˴�Ĭ��Ϊ3���ھ� ! ����������
		for (int j = 0; j < num_neighbor; j++) {
			if (neighbors_element_i[j] == nullptr) {
				element_host.neighbors[j].isNull[i] = true;
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
			}
		}
	}

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
	// ����host���ݵ�device
	element_host.cuda_copy_to_device(&element_device);

	// --- ���㵥Ԫ�ݶ� --- 
	// ���룺��Ԫ���ꡢ��ԪU����Ԫ�ھ����ꡢ��Ԫ�ھ�U
	// �������Ԫ�ݶ�
	GPU::calculateGradient(this->element_device);
	// ����device���ݵ�host
	element_device.cuda_copy_to_host(&element_host);
	// ������ Ŀǰ����Ҫ��
	
	// TODO: �����ھ��ݶ� 1.��ElementDataPack�����neighborIndex����

	// --- ����߽���ֵͨ�� --- 

}

inline void GPUSolver::finalize() {
	this->element_host.free(); 
	this->element_device.cuda_free();
}
