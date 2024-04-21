#ifndef _GPU_CDRIVER_H_
#define _GPU_CDRIVER_H_
// �������FVM_2D���ֹ���
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

	public:
		enum PauseSignal {
			_NoSignal,// ����ͣ�ź�
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
			PauseSignal pauseSignal = _NoSignal;// ��ͣ�ź�
		};

	public:
		static CDriver* getInstance();
		void run();
		static void saveAndExit(int _Code);
		
	private:
		CDriver() {};
		SignalPack emitSignalPack(int istep, int maxIteration, real t, real T, real residualRho);
		void onSignalPack(const SignalPack& sp);
		void updateOldData(int istep, double t) {
			m_last_istep = istep;
			m_last_t = t;
		}
		void m_saveAndExit(int _Code);
	};
}

#endif