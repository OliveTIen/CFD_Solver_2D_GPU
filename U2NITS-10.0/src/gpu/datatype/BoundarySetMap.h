
#ifndef _BOUNDARY_SET_MAP_H_
#define _BOUNDARY_SET_MAP_H_
#include "Define.h"
namespace GPU {
	class BoundarySetMap {
	public:
		int size = 0;
		int* type{};
	public:
		void alloc(int _size) {
			size = _size;
			type = new int[size];
		}
		void free() {
			delete[] type;
		}
		void cuda_alloc(int _size) {
			size = _size;
			cudaMalloc(&type, size * sizeof(int));
		}

		void cuda_free() {
			cudaFree(type);
		}

		static void cuda_memcpy(BoundarySetMap* dist, const BoundarySetMap* src, cudaMemcpyKind kind);

	};
}

#endif // !_BOUNDARY_SET_MAP_H_
