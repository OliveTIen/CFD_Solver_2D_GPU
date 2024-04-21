#ifndef GPU_SOLVER_2_H
#define GPU_SOLVER_2_H

#include "../gpu/datatype/Datatype.h"
#include <map>

namespace GPU {
	class GPUSolver2 {
	private:
		bool hostMemoryAllocated = false;
		bool hasSetGPUDevice = false;// �豸�������Ѿ���
		bool deviceMemoryAllocated = false;
		bool iterationStarted = false;
		bool hostDataInitialized = false;
		bool deviceDataInitialized = false;

	public:
		// host/device˫������
		GPU::NodeSoA node_host;
		GPU::NodeSoA node_device;
		GPU::ElementSoA element_host;
		GPU::ElementSoA element_device;
		GPU::FieldSoA elementField_host;
		GPU::FieldSoA elementField_device;
		GPU::EdgeSoA edge_host;
		GPU::EdgeSoA edge_device;
		GPU::BoundarySetMap boundary_host;// ���Ա���boundary��setIDӳ�䵽type
		GPU::BoundarySetMap boundary_device;

		// host����
		REAL* element_vruvp[4]{};// �洢Uת��Ϊ��ruvp�����ڼ���ڵ�ruvp�Լ�����ʱ�䲽��Dt
		REAL* element_U_old[4]{};// �洢��һ��U���ڼ���в�, 4xn����
		const int residualVectorSize = 4;
		real residualVector[4]{ 1,1,1,1 };// �вresidualVectorSize x 1����
		std::map<int, int> edge_periodic_pair;// �洢���ڱ߽�ı߽��
		GPU::OutputNodeFieldSoA outputNodeField;

		// device����
		GPU::DGlobalPara* infPara_device;// Զ��������GPU�ĸ��������ڴ���[RAII����]
		SDevicePara sDevicePara;

	public:
		// �����ڴ档Ӧ�ڶ�ȡfield file��ʹ�ã���ΪҪȷ���ж��ٵ�Ԫ����Ҫ����ڴ�
		void allocateMemory();
		// ��FVM_2D��ʼ������
		void initializeData_byOldData();
		// ����
		void iteration(real& t, real T);
		// ���½ڵ㳡
		void updateOutputNodeField();
		// ���²в�
		void updateResidual();
		// ��ȡ�в�
		void getResidual();
		// �ͷ���Դ
		void freeMemory();

		bool isIterationStarted() { return iterationStarted; }

		// device���ݸ��µ�host
		void device_to_host();
		// host���ݸ��µ�device
		void host_to_device();
	private:
		void setGPUDevice();
		// ��FVM_2D�����ݳ�ʼ��hostData
		void initializeHostData_byOldData();
		// ��hostData��ʼ��deviceData
		void initializeDeviceData_byHostData();
		void initialize_nodeHost(void* _pFVM2D_, int num_node);
		// ��ʼ��element��elementField��host����
		void initialize_elementHost(void* _pFVM2D_, int num_element);
		// ��ʼ��edge��host����
		void initialize_edgeHost(void* _pFVM2D_, int num_edge);
		// ��ʼ���߽����� �����ڱ߽��Ϊ�ڲ��߽�
		void initialize_boundary(void* _pFVM2D_);


		void allocateMemory(const int num_node, const int num_element, const int num_edge, const int num_boundary);
		void allocateHostMemory(const int num_node, const int num_element, const int num_edge, const int num_boundary);
		void allocateDeviceMemory(const int num_node, const int num_element, const int num_edge, const int num_boundary);// ��Ҫ��setGPUDevice֮��
		void freeHostMemory();
		void freeDeviceMemory();

	};
}

#endif // GPU_SOLVER_2_H