#ifndef _GPU_CDRIVER_H_
#define _GPU_CDRIVER_H_
// 用于替代FVM_2D部分功能
#include "../solvers/GPUSolver2.h"
#include "../input/CInput.h"
#include "../output/COutput.h"
#include "../gpu/datatype/DefineType.h"

namespace U2NITS {
	class CDriver {
	private:
		static CDriver* pCDriver;
		CInput input;
		COutput out;
		GPU::GPUSolver2 solver;
		int m_last_istep = 0;
		double m_last_t = 0.0;
		double m_last_snapshot_time = 0.0;// 上次保存快照的时间，默认为0
		double m_snapshot_cold_time = 1800.0;// 保存快照的冷却时间[秒]，默认为30min=1800s

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
			bool b_writeContinue = false;// 续算
			bool b_writeTecplot = false;// 流场
			bool b_writeRecovery = false;// [已弃用]因为不如保存续算文件，所以没用到
			bool b_nanDetected = false;
			bool b_writeHist = false;// 残差
			bool b_print = false;// 屏幕
			PauseSignal pauseSignal = _NoSignal;// 暂停信号
		};

	public:
		static CDriver* getInstance();
		// 主循环
		void start();
		static void saveAndExit(int _Code);
		
	private:
		CDriver() {};
		void modifySignalPack1_output(SignalPack& s, int istep, int maxIteration, double current_time);
		void modifySignalPack2_pause(SignalPack& s, int istep, int maxIteration, myfloat t, myfloat T, myfloat residualRho);
		// 存储当前istep和t
		void updateOldData_istep_t(int istep, double t) {
			m_last_istep = istep;
			m_last_t = t;
		}
		void m_saveAndExit(int _Code);
	};
}

#endif