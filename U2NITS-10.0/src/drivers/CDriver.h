#ifndef _GPU_CDRIVER_H_
#define _GPU_CDRIVER_H_
// 用于替代FVM_2D部分功能
namespace U2NITS {
	class CDriver {
	public:
		static void run_current();
		static void solve_current();
		static void run_2();

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