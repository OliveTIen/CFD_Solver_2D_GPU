#ifndef GPU_BOUNDARY_V2_H
#define GPU_BOUNDARY_V2_H
#include "DefineType.h"
#include <vector>
#include <string>

namespace GPU {
	/*
	第2版边界
	20240426
	
	*/
	struct BoundaryV2 {
		// 类内定义 边集
		struct EdgeSet {
			int type = -1;// 类型，初始化时根据名称确定
			int ID = -1;// 等于数组指标+1，默认值是-1。为兼容旧版而设置，最好不要使用，请直接用数组下标访问。
			std::string name;// 名称，例如far, wall, foil
			std::vector<myint> edge_vector;// edge index，数组，存储GPUID
		};

		std::vector<EdgeSet> edgeSets;// 边集合数组
	};
}

#endif // !GPU_BOUNDARY_V2_H
