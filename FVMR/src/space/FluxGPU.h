#ifndef _CALCULATE_FLUX_2_CUH_
#define _CALCULATE_FLUX_2_CUH_

#include "../gpu/datatype/Datatype.h"
#include <map>
// __host__: �����������������������á���Ĭ��ֵ
// __global__: �豸�������������������ã��ҵ���ʱҪָ��<<<>>>
// __device__: �豸���������豸��������

namespace GPU {
	namespace Space {
		// δ���
		namespace Flux {
			void calculateFluxDevice_2(ElementSoA& element_device, EdgeSoA& edge_device, ElementFieldSoA& elementField_device);
		}
	}
}


#endif