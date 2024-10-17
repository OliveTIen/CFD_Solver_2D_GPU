#pragma once
#include <vector_types.h>
#include <cuda_runtime_api.h>

namespace GPU {
	struct ColorMap {
	public:
		// data: {position, r, g, b}, {position, r, g, b}, ...
		int num_control_point = 0;
		float4* data = nullptr;

		void alloc(int _num_point) {
			num_control_point = _num_point;
			data = new float4[num_control_point];
		}

		void free() {
			delete[] data;
		}

		void cuda_alloc(int _num_point) {
			num_control_point = _num_point;
			cudaMalloc((void**) & data, num_control_point * sizeof(float4));
		}

		void cuda_free() {
			cudaFree(data);
		}

		static void cuda_memcpy(ColorMap* dist, const ColorMap* src, cudaMemcpyKind kind);

		bool has(int point_index) {
			return(point_index >= 0 && point_index < num_control_point);
		}
	};
}