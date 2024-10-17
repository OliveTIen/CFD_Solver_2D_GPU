
#ifndef _BOUNDARY_SET_MAP_H_
#define _BOUNDARY_SET_MAP_H_
#include "DefineType.h"
#include "../Env.h"
namespace GPU {
	/*
	�洢�ӱ߽�ID���߽�����type��ӳ�䡣���洢��������Щedge
	*/
	struct BoundarySetMap {
	public:
		myint size = 0;

		myint* type = nullptr;

		// ������������ɡ��洢�����߽�����Щedge�����ڱ����߽�ߣ���������������
		//myint* numOfEdge = nullptr;
		//myint** edge = nullptr;

		// ��Ҫ�����Ҫ��ָ�룬���αȽ��鷳����Ҫ��num of boundary set����Ҫ��ÿ��set��edge������������
		// �Ҳ�����������edge�������ڴ棬��ȵ�����ʱ�ٴ���
		// Ŀǰ����Э��ʽ������BoundaryManager_2D

		/*
		����
		BoundarySetMap�洢std::vector����


		���⣺�����std::vector������ͨ��ʱ����֪edge����λ�ȡ�����߽����ͣ�
		�𰸣��Ѿ���edge_host�д洢��setID��
			���ò�ѯedge_host��ȡsetID��
			int setID = edge_host.setID[iedge];
			����setID��ȡ�߽����ͣ�
			bType = f->boundaryManager.boundaries[setID - 1].type;
		
		���⣺�����std::vector��GPU��δ��Σ�
		�𰸣�GPU�п�����thrust::vector
			�Ȳ���GPU��ֻ��CPU
		*/
	public:
		void alloc(myint _size) {
			size = _size;

			type = new myint[size];
			//numOfEdge = new myint[size];
			//for (myint i = 0; i < size; i++) {
			//	numOfEdge[i]=new myint
			//}
		}
		void free() {

			delete[] type;
			//delete[] numOfEdge;
		}
		void cuda_alloc(myint _size) {
			size = _size;

			cudaMalloc(&type, size * sizeof(myint));
			//cudaMalloc(&numOfEdge, size * sizeof(myint));
		}

		void cuda_free() {

			cudaFree(type);
			//cudaFree(numOfEdge);
		}

		static void cuda_memcpy(BoundarySetMap* dist, const BoundarySetMap* src, cudaMemcpyKind kind);

	};
}

#endif // !_BOUNDARY_SET_MAP_H_
