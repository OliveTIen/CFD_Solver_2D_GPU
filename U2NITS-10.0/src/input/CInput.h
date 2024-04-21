#ifndef _CINPUT_H_
#define _CINPUT_H_
#include "../global/SGlobalPara.h"

namespace U2NITS {
	class CInput {
	public:
		//SGlobalPara::SGlobalPara para;

		void readConfig();
		void readField_1();

		// [todo]直接读取到GPU::ElementSoA等，不经过FVM_2D
		void readField_2_unused();

	private:
		
	};
}

#endif // !_CINPUT_H_
