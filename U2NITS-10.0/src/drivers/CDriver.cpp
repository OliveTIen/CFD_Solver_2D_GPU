#include "CDriver.h"
#include "../FVM_2D.h"
#include "../math/Common.h"
#include "../global/CExit.h"

U2NITS::CDriver* U2NITS::CDriver::pCDriver = nullptr;

U2NITS::CDriver* U2NITS::CDriver::getInstance() {
	if (pCDriver == nullptr) {
		pCDriver = new CDriver;
	}
	return pCDriver;
}



void U2NITS::CDriver::run() {

	out.console.printWelcomeInterface();
	input.readConfig();
	input.readField_1();
	solver.allocateMemory();
	solver.initializeData_byOldData();
	// ���ֿ��Ʋ��������á���Ҫ��readConfig��ʹ�ã�����δ��ʼ��
	const int istep_start = GlobalPara::time::istep_previous;// �Ѽ��㲽
	const int maxIteration = GlobalPara::output::maxIteration;// Ŀ�경
	int istep = istep_start;
	double t = GlobalPara::time::t_previous;// ��ǰ����ʱ��
	const double T = GlobalPara::time::max_physical_time;

	// ����
	out.hist.setFilePath(out.getDir() + GlobalPara::basic::filename + "_hist.dat");
	out.timer.start();// ��ʱ
	out.console.setClearStartPosition();
	while (true) {// ѭ����ֹ������signal
		istep++;
		// �����ļ���
		out.updateFileName(istep);
		solver.iteration(t, T);
		std::cout << ".";// �����.����ʾ���ڼ��㣬û�п���

		SignalPack sp = emitSignalPack(istep, maxIteration, t, T, solver.residualVector[0]);

		// ��ӡ����
		if (sp.b_print) {
			solver.getResidual();
			out.console.clear();
			out.console.drawProgressBar((double)istep, (double)maxIteration, 45);
			out.console.setSolveInfo(istep_start, istep, maxIteration, out.field.getNumTecplotFileWritten(), out.timer.getIime(), t, T, solver.residualVector);
			out.console.print(out.console.getSolveInfo());
		}
		// д�ļ�����
		if (sp.b_writeContinue || sp.b_writeTecplot) {
			solver.updateOutputNodeField();// ���ݵ�ԪU���½ڵ�ruvp
			if (sp.b_nanDetected) {
				out.console.printInfo(out.console.InfoType::type_nan_detected);
				out.continueFilePath = out.continueFilePath_nan;
			}
			if (sp.b_writeContinue) {
				out.field.writeContinueFile_1(
					istep, t, out.continueFilePath,
					solver.node_host, solver.element_host, solver.elementField_host.U
				);
			}
			if (sp.b_writeTecplot) {
				out.field.writeTecplotFile_1(t, out.tecplotFilePath, "title", solver.node_host, solver.element_host, solver.outputNodeField);
			}
		}
		if (sp.b_writeHist) {
			// ��ȡ�в���������residual_vector��
			solver.getResidual();
			out.hist.writeHistFile_2(istep, solver.residualVector, solver.residualVectorSize);
		}
		// ��ֹ����
		if (sp.pauseSignal != _NoSignal) {
			out.log.log(out.console.getSolveInfo());
			out.log.logAndPrintSignalInfo(sp.pauseSignal);// ��ʾ�� ����д�ļ�֮�󣬷�ֹ����
			break;// ���ɷŽ����������У���ΪҪ��ֹѭ��
		}
		updateOldData(istep, t);
	}

	// �ͷ���Դ
	solver.freeMemory();

	// ֹͣ
	CExit::pressAnyKeyToExit();

	/*
���裺���Ż���


������ʼ
	���㣨�����Ͳв CSolver.iterate() �в�������ʱ����
	���������Ϣ�������ļ�������һ�ι��λ�á�������ȣ�COutput.update()����(COutput.updateData COutput.printScreen())
	����ļ����������вCOutput.writeFile() ����COutput.writeField() COutput.writeHist()(����CSolver.updateResidual() )
	�ж��Ƿ�����ѭ������ֹ��this->checkStop() ��Ϣ���У�
	���������Ϣ����ʾ�ʣ�COutput.updateScreen()
�ͷ���Դ CSolver.finalize()

	*/

}

void U2NITS::CDriver::saveAndExit(int _Code) {
	getInstance()->m_saveAndExit(_Code);
}

void U2NITS::CDriver::m_saveAndExit(int _Code) {
	if (!solver.isIterationStarted()) {// ����δ��ʼ����������ļ�
		exit(_Code);
	}

	out.updateFileName(m_last_istep);
	out.field.writeContinueFile_1(
		m_last_istep, m_last_t, out.continueFilePath,
		solver.node_host, solver.element_host, solver.element_U_old
	);
	exit(_Code);
}


U2NITS::CDriver::SignalPack U2NITS::CDriver::emitSignalPack(int istep, int maxIteration, real t, real T, real residualRho) {
	SignalPack s;

	// �����ź�
	const int offset = 0;
	const double big_value = 1e8;

	//PauseSignal pauseSignal = _NoSignal;// ��ͣ�ź� ��ʾ�ʷ���
	if (istep % GlobalPara::output::step_per_print == offset) {// �����������
		s.b_print = true;
	}
	if (istep % GlobalPara::output::step_per_output_field == offset) {// �����������
		s.b_writeTecplot = true;
		s.b_print = true;
	}
	if (istep % GlobalPara::output::step_per_output_hist == offset) {// ��������в�
		s.b_writeHist = true;
	}
	if (residualRho > big_value) {// �в������ֹ
		s.b_writeContinue = true;
		s.b_writeTecplot = true;
		s.b_nanDetected = true;
		s.pauseSignal = _NanDetected;
	}
	else if (_kbhit()) {// ��esc��ֹ
		if (_getch() == 27) {
			s.b_writeContinue = true;
			s.pauseSignal = _EscPressed;
			std::cout << "Stopping...\n";
		}
	}
	else if (t > T || istep >= maxIteration) {// �ﵽ�涨����ʱ��������������ֹ
		s.b_writeContinue = true;
		s.b_writeTecplot = true;
		s.pauseSignal = _TimeReached;
	}
	else if (residualRho <= U2NITS::Math::EPSILON) {// �в��㹻С����ֹ
		s.b_writeContinue = true;
		s.b_writeTecplot = true;
		s.pauseSignal = _StableReached;
	}

	return s;
}

void U2NITS::CDriver::onSignalPack(const SignalPack& sp) {

}

