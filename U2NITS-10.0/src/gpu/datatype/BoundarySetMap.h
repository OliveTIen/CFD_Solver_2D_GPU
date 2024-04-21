
#ifndef _BOUNDARY_SET_MAP_H_
#define _BOUNDARY_SET_MAP_H_
#include "DefineType.h"
#include "../Env.h"
namespace GPU {
	struct BoundarySetMap {
	public:
		integer size = 0;

		integer* type = nullptr;

		// ������������ɡ��洢�����߽�����Щedge�����ڱ����߽�ߣ���������������
		//integer* numOfEdge = nullptr;
		//integer** edge = nullptr;

		// ��Ҫ�����Ҫ��ָ�룬���αȽ��鷳����Ҫ��num of boundary set����Ҫ��ÿ��set��edge������������
		// �Ҳ�����������edge�������ڴ棬��ȵ�����ʱ�ٴ���
		// Ŀǰ����Э��ʽ������BoundaryManager_2D
	public:
		void alloc(integer _size) {
			size = _size;

			type = new integer[size];
			//numOfEdge = new integer[size];
			//for (integer i = 0; i < size; i++) {
			//	numOfEdge[i]=new integer
			//}
		}
		void free() {

			delete[] type;
			//delete[] numOfEdge;
		}
		void cuda_alloc(integer _size) {
			size = _size;

			cudaMalloc(&type, size * sizeof(integer));
			//cudaMalloc(&numOfEdge, size * sizeof(integer));
		}

		void cuda_free() {

			cudaFree(type);
			//cudaFree(numOfEdge);
		}

		static void cuda_memcpy(BoundarySetMap* dist, const BoundarySetMap* src, cudaMemcpyKind kind);

	};
}

#endif // !_BOUNDARY_SET_MAP_H_
