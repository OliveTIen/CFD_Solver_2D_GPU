#ifndef _GPU_MAIN_H_
#define _GPU_MAIN_H_
// �������FVM_2D���ֹ���
namespace GPU {
	class GPUMain {
	public:
		static void run_GPU();
		static void solve_GPU();

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