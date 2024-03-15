#ifndef _GPU_CDRIVER_H_
#define _GPU_CDRIVER_H_
// �������FVM_2D���ֹ���
namespace U2NITS {
	class CDriver {
	public:
		static void run_GPU();
		static void solve_GPU();
		static void run_GPU_2();

	private:
		enum MySignal {
			_NoSignal,
			_EscPressed,
			_TimeReached,
			_StableReached
		};
		static void printSignalInfo(MySignal signal);
	};
}

#endif