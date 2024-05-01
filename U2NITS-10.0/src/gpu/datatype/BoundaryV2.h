#ifndef GPU_BOUNDARY_V2_H
#define GPU_BOUNDARY_V2_H
#include "DefineType.h"
#include <vector>
#include <string>

namespace GPU {
	/*
	��2��߽�
	20240426
	
	*/
	struct BoundaryV2 {
		// ���ڶ��� �߼�
		struct EdgeSet {
			int type = -1;// ���ͣ���ʼ��ʱ��������ȷ��
			int ID = -1;// ��������ָ��+1��Ĭ��ֵ��-1��Ϊ���ݾɰ�����ã���ò�Ҫʹ�ã���ֱ���������±���ʡ�
			std::string name;// ���ƣ�����far, wall, foil
			std::vector<myint> edge_vector;// edge index�����飬�洢GPUID
		};

		std::vector<EdgeSet> edgeSets;// �߼�������
	};
}

#endif // !GPU_BOUNDARY_V2_H
