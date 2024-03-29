#ifndef _CINPUT_H_
#define _CINPUT_H_

namespace U2NITS {
	class CInput {
	public:

		void readConfig();
		void readField_old();

		// [todo]直接读取到GPU::ElementSoA等，不经过FVM_2D
		void readField_2();

	private:
		
	};
}

#endif // !_CINPUT_H_
