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
		double m_last_snapshot_time = 0.0;// �ϴα�����յ�ʱ�䣬Ĭ��Ϊ0
		double m_snapshot_cold_time = 1800.0;// ������յ���ȴʱ��[��]��Ĭ��Ϊ30min=1800s

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
			bool b_writeContinue = false;// ����
			bool b_writeTecplot = false;// ����
			bool b_writeRecovery = false;// [������]��Ϊ���籣�������ļ�������û�õ�
			bool b_nanDetected = false;
			bool b_writeHist = false;// �в�
			bool b_print = false;// ��Ļ
			PauseSignal pauseSignal = _NoSignal;// ��ͣ�ź�
		};

	public:
		static CDriver* getInstance();
		// ��ѭ��
		void start();
		static void saveAndExit(int _Code);
		
	private:
		CDriver() {};
		void modifySignalPack1_output(SignalPack& s, int istep, int maxIteration, double current_time);
		void modifySignalPack2_pause(SignalPack& s, int istep, int maxIteration, myfloat t, myfloat T, myfloat residualRho);
		// �洢��ǰistep��t
		void updateOldData_istep_t(int istep, double t) {
			m_last_istep = istep;
			m_last_t = t;
		}
		void m_saveAndExit(int _Code);
	};
}

#endif