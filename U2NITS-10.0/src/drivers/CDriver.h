#ifndef _GPU_CDRIVER_H_
#define _GPU_CDRIVER_H_
// 用于替代FVM_2D部分功能
#include "../gpu/datatype/Define.h"

namespace U2NITS {
	class CDriver {
	public:
		enum PauseSignal {
			_NoSignal,// 无暂停信号
			_EscPressed,
			_TimeReached,
			_StableReached,
			_NanDetected
		};
		class SignalPack {
		public:
			bool b_writeContinue = false;
			bool b_writeTecplot = false;
			bool b_nanDetected = false;
			bool b_writeHist = false;
			bool b_print = false;
			PauseSignal pauseSignal = _NoSignal;// 暂停信号
		};

		static void run_current();
		
	private:
		static SignalPack emitSignalPack(int istep, int maxIteration, real t, real T, real residualRho);
		static void onSignalPack(const SignalPack& sp);
	};
}

#endif