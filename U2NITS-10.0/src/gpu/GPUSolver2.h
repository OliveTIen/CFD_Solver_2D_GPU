#ifndef GPU_SOLVER_2_H
#define GPU_SOLVER_2_H

//#include "dataType/ElementDataPack.h"
//#include "dataType/ElementAdjacent.h"
#include "datatype/NodeSoA.h"
#include "dataType/ElementSoA.h"
#include "datatype/FieldSoA.h"
#include "dataType/EdgeSoA.h"
#include <map>

namespace GPU {
	class GPUSolver2 {
	public:
		GPU::NodeSoA node_host;
		GPU::NodeSoA node_device;
		GPU::ElementSoA element_host;
		GPU::ElementSoA element_device;
		GPU::FieldSoA elementField_host;
		GPU::FieldSoA elementField_device;
		GPU::EdgeSoA edge_host;
		GPU::EdgeSoA edge_device;

		REAL* element_vruvp[4]{};// �洢Uת��Ϊ��ruvp�����ڼ���ڵ�ruvp�Լ�����ʱ�䲽��Dt
		REAL* element_U_old[4]{};// �洢��һ��U���ڼ���в�
		std::map<int, int> edge_periodic_pair;// �洢���ڱ߽�ı߽��
		GPU::OutputNodeFieldSoA outputNodeField;
		//GPU::PeriodicBoundary periodicBoundary;// ����ʱ��;���Ժ���Ż��� 
		/* ռ�ݴ����ռ䣬��Ȼ���Ҹ��Ӷ�ΪO(1)���������0Ԫ��
		*/

	public:
		// ������Դ����ʼ��GPUID��host���Ƶ�device
		void initialze();
		void iteration(REAL& t, REAL T);
		// ���½ڵ㳡
		void updateOutputNodeField();
		// �ͷ���Դ
		void finalize();
		
	private:
		void initialize_nodeHost(void* _pFVM2D_, int num_node);
		// ��ʼ��element��elementField��host����
		void initialize_elementHost(void* _pFVM2D_, int num_element);
		// ��ʼ��edge��host����
		void initialize_edgeHost(void* _pFVM2D_, int num_edge);
		// ��ʼ���߽����� �����ڱ߽��Ϊ�ڲ��߽�
		void initialize_boundary(void* _pFVM2D_);

		// ��host�˵ĸ��� ��U����ruvp U_old
		void update_ruvp_Uold();
		// ����ʱ�䲽��
		REAL calculateDt(REAL t, REAL T);
		// GPU iteration�� ��������GPU kernel
		void iterationDevice(REAL dt);
		// device���ݸ��µ�host
		void device_to_host();
		// ��host iteration
		void iterationHost(REAL dt);
		// host���ݸ��µ�device
		void host_to_device();
	};
}

#endif // GPU_SOLVER_2_H