#ifndef _CINPUT_H_
#define _CINPUT_H_
#include "../global/SGlobalPara.h"

namespace U2NITS {
	class CInput {
	public:
		//SGlobalPara::SGlobalPara para;

		void readConfig();
		void readField_1();

		// [todo]ֱ�Ӷ�ȡ��GPU::ElementSoA�ȣ�������FVM_2D
		void readField_2_unused();

	private:
		
	};
}

#endif // !_CINPUT_H_
