#ifndef _OLDDATACONVERTER_H_
#define _OLDDATACONVERTER_H_

#include "GPUSolver2.h"
#include "../FVM_2D.h"

namespace U2NITS {
	/**
	* 旧数据类型到新数据的转换
	* 
	*/
	class OldDataConverter {
	private:
		FVM_2D* m_pFVM2D = nullptr;
		GPU::GPUSolver2* m_pGPUSolver2 = nullptr;

	public:
		OldDataConverter(FVM_2D* pF, GPU::GPUSolver2* pG) {
			m_pFVM2D = pF;
			m_pGPUSolver2 = pG;
		}
		void Convert_FVM2D_to_HostData();

	private:
		void InitializeOldGPUIndex(int num_node, int num_element, int num_edge);
		void ConvertNode(int num_node);
		void ConvertElement(int num_element);
		void ConvertEdge(int num_edge);
		void ConvertBoundary(int num_boundary);
	};
}

#endif // !_OLDDATACONVERTER_H_
