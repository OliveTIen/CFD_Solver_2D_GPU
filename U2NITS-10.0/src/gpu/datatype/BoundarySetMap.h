
#ifndef _BOUNDARY_SET_MAP_H_
#define _BOUNDARY_SET_MAP_H_
#include "Define.h"
namespace GPU {
	struct BoundarySetMap {
	public:
		int size = 0;
		int* _size_;// GPU π”√

		int* type{};
	public:
		void alloc(int _size) {
			size = _size;
			_size_ = new int;
			*_size_ = size;

			type = new int[size];
		}
		void free() {
			delete _size_;

			delete[] type;
		}
		void cuda_alloc(int _size) {
			size = _size;
			cudaMalloc(&_size_, 1 * sizeof(int));
			cudaMemcpy(_size_, &size, 1 * sizeof(int), ::cudaMemcpyHostToDevice);

			cudaMalloc(&type, size * sizeof(int));
		}

		void cuda_free() {
			cudaFree(_size_);

			cudaFree(type);
		}

		static void cuda_memcpy(BoundarySetMap* dist, const BoundarySetMap* src, cudaMemcpyKind kind);

	};
}

#endif // !_BOUNDARY_SET_MAP_H_
