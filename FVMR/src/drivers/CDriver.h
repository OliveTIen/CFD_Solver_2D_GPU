#ifndef _GPU_CDRIVER_H_
#define _GPU_CDRIVER_H_
// 用于替代FVM_2D部分功能
#include "../solvers/GPUSolver2.h"
#include "../input/CInput.h"
#include "../output/COutput.h"
#include "../gpu/datatype/DefineType.h"
#include "../gui/WinGui.h"


namespace U2NITS {
	/**
	 * \class CDriver
	 * \brief 
	 * Program routine controller. 
	 * \details @ref CDriver::start controls the whole routine,
	 * and @ref CDriver::gui controls the physics simulation scene.<br>
	 * TODO:<br>
	 * to make it a game engine, it's suggested to attach the solver to a CudaGrid Object
	 */
	class CDriver {
	private:
		static CDriver* pCDriver;
		CInput input;
		COutput out;
		GPU::GPUSolver2 solver;
		WinGui gui;///< controls the physics simulation scene.
		int m_last_istep = 0;
		double m_last_t = 0.0;
		double m_last_snapshot_time = 0.0;
		double m_snapshot_cold_time = 1800.0;///< cold time for snapshot, default 30 minutes

	public:
		enum PauseSignal {
			_NoSignal,
			_EscPressed,
			_TimeReached,
			_StableReached,
			_NanDetected
		};
		/// <summary>
		/// struct storing types of pause signals
		/// </summary>
		class SignalPack {
		public:
			bool b_writeContinue = false;///< 
			bool b_writeTecplot = false;///< 
			bool b_writeRecovery = false;///< [unused]
			bool b_nanDetected = false;
			bool b_writeHist = false;///< hist file such as residual
			bool b_print = false;///< 
			PauseSignal pauseSignal = _NoSignal;///< 
		};

		bool show_gui = true;///< Show graphic interface. Read from input.toml
		bool enable_print = true;///< Enable console output. Read from input.toml
		bool enable_write_file = true;///< Enable field file output (it doesn't affect log file). Read from input.toml

		static CDriver* getInstance();
		void start();///< \brief entry point
		static void saveAndExit(int _Code);
		
	private:
		CDriver() {};
		void modifySignalPack1_output(SignalPack& s, int istep, int maxIteration, double current_time);
		void modifySignalPack2_pause(SignalPack& s, int istep, int maxIteration, myfloat t, myfloat T, myfloat residualRho);
		void updateOldData_istep_t(int istep, double t) {
			m_last_istep = istep;
			m_last_t = t;
		}
		void m_saveAndExit(int _Code);
	};
}

#endif