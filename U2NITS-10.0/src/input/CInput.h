#ifndef _CINPUT_H_
#define _CINPUT_H_

namespace U2NITS {
	class CInput {
	public:

		void readConfig();
		void readField_old();

		// [todo]ֱ�Ӷ�ȡ��GPU::ElementSoA�ȣ�������FVM_2D
		void readField_2();

	private:
		
	};
}

#endif // !_CINPUT_H_
