#include "CDriver.h"
#include "../FVM_2D.h"
#include "../math/Common.h"
#include "../global/CExit.h"
#include "RoutineController.h"
#include "../time/CIteration.h"

U2NITS::CDriver* U2NITS::CDriver::pCDriver = nullptr;

U2NITS::CDriver* U2NITS::CDriver::getInstance() {
	if (pCDriver == nullptr) {
		pCDriver = new CDriver;
	}
	return pCDriver;
}


void U2NITS::CDriver::run_20240517() {
	/*
	GPU����ƿ�����ڿ���ռ��50%ʱ��
	��Ҫÿ����device_to_host��ֻ����Ҫд�ļ�ʱ��device_to_host
	�����Ҫ��device_to_host����solver.updateResidualVector();֮ǰ
	*/

	out.console.printWelcomeInterface();
	input.readConfig();
	input.readField_1();
	solver.allocateMemory();
	solver.initializeData_byOldData();
	// ���ֿ��Ʋ��������á���Ҫ��readConfig��ʹ�ã�����δ��ʼ��
	const int istep_start = GlobalPara::time::istep_previous;// �Ѽ��㲽
	const int maxIteration = GlobalPara::output::maxIteration;// Ŀ�경
	int istep = istep_start;
	myfloat t = GlobalPara::time::t_previous;// ��ǰ����ʱ��
	const myfloat T = GlobalPara::time::max_physical_time;

	// ����
	out.hist.setFilePath(out.getDir() + GlobalPara::basic::filename + "_hist.dat");
	out.set_tecplot_hist_path(out.getDir() + GlobalPara::basic::filename + "_hist.dat");
	out.timer.start();// ��ʱ
	out.console.setClearStartPosition();
	while (true) {// ѭ����ֹ������signal
		istep++;
		// �����ļ���
		out.updateFileName(istep);

		if (GlobalPara::basic::useGPU) {
			CIteration::iteration_device_20240517(t, T, &solver);
		}
		else {
			CIteration::iteration_host(t, T, &solver);
		}

		if (!solver.isIterationStarted()) {
			solver.setIterationStarted(true);
		}

		std::cout << ".";// �����.����ʾ���ڼ��㣬û�п���

		SignalPack sp;
		modifySignalPack1_output(sp, istep, maxIteration);
		modifySignalPack2_pause(sp, istep, maxIteration, t, T, solver.residualVector[0]);

		/*
		�´��� ����Ҫ��ӡ������ļ����򿽱��豸������
		��Ҫ���������ݣ�
		elementField_host ����в�residualVector������ڵ����ֵnodeField_host(����tecplot)�����continuefile/recoveryfile
		*/ 
		if (GlobalPara::basic::useGPU) {
			if (sp.b_print ||
				sp.b_writeContinue || sp.b_writeTecplot || sp.b_writeRecovery || sp.b_writeHist) {
				GPU::ElementFieldSoA::cuda_memcpy(&solver.elementField_host, &solver.elementField_device, cudaMemcpyDeviceToHost);
			}
		}

		// ��ӡ����
		if (sp.b_print) {
			solver.updateResidualVector();// ����elementField_host����residualVector
			// ����CFL�����������dtʱ���õ�
			RoutineController::getInstance()->applyStrategy(solver.residualVector[0], GlobalPara::time::CFL, istep);
			out.console.clear();
			out.console.drawProgressBar((double)istep, (double)maxIteration, 45);
			out.console.setSolveInfo(istep_start, istep, maxIteration, out.getFieldWriter()->getNumTecplotFileWritten(), out.timer.getTime(), t, T, solver.residualVector, GlobalPara::time::CFL);
			out.console.print(out.console.getSolveInfo());
		}
		// д�ļ�����
		if (sp.b_writeContinue || sp.b_writeTecplot || sp.b_writeRecovery || sp.b_writeHist) {
			out.console.logSpeed(istep, istep_start, out.timer.getTime());
			solver.updateOutputNodeField();// ����elementField_host����nodeField_host
			solver.updateResidualVector();
			if (sp.b_nanDetected) {
				out.console.printInfo(out.console.InfoType::type_nan_detected);
				out.continueFilePath = out.continueFilePath_nan;
			}
			if (sp.b_writeContinue) {
				out.getFieldWriter()->writeContinueFile_1(istep, t, out.continueFilePath, solver.node_host, solver.element_host, solver.elementField_host.U);
			}
			if (sp.b_writeTecplot) {
				out.getFieldWriter()->write_tecplot_volume_file(t, out.tecplotFilePath, "title", solver.node_host, solver.element_host, solver.outputNodeField);
				out.getFieldWriter()->write_tecplot_boundary_file(t, out.tecplotBoundaryFilePath, solver.node_host, solver.edge_host, solver.outputNodeField, solver.boundary_host_new);
			}
			if (sp.b_writeRecovery) {
				out.getFieldWriter()->writeContinueFile_1(istep, t, out.recoveryFilePath, solver.node_host, solver.element_host, solver.elementField_host.U);
			}
			if (sp.b_writeHist) {
				out.getFieldWriter()->write_tecplot_hist_file(out.tecplot_hist_path, istep, t, solver.residualVector, solver.edge_host, solver.outputNodeField, solver.boundary_host_new);
			}
		}
		// ��ֹ����
		if (sp.pauseSignal != _NoSignal) {
			out.log.log(out.console.getSolveInfo());
			out.log.logAndPrintSignalInfo(sp.pauseSignal);// ��ʾ�� ����д�ļ�֮�󣬷�ֹ����
			break;// ���ɷŽ����������У���ΪҪ��ֹѭ��
		}
		updateOldData_istep_t(istep, t);
	}

	// �ͷ���Դ
	solver.freeMemory();

	// ֹͣ
	CExit::pressAnyKeyToExit();

}

void U2NITS::CDriver::saveAndExit(int _Code) {
	getInstance()->m_saveAndExit(_Code);
}

void U2NITS::CDriver::m_saveAndExit(int _Code) {
	if (!solver.isIterationStarted()) {// ����δ��ʼ����������ļ�
		exit(_Code);
	}

	out.updateFileName(m_last_istep);
	out.getFieldWriter()->writeContinueFile_1(
		m_last_istep, m_last_t, out.continueFilePath,
		solver.node_host, solver.element_host, solver.elementField_host.Uold
	);
	exit(_Code);
}

void U2NITS::CDriver::modifySignalPack1_output(SignalPack& s, int istep, int maxIteration) {
	/*
	
	*/

	// �����ź�
	const int offset = 0;

	if (istep % GlobalPara::output::step_per_print == offset) {
		// ��Ļ
		s.b_print = true;
		// ��ԭ��
		s.b_writeRecovery = true;
	}
	if (istep % GlobalPara::output::step_per_output_field == offset &&
		istep >= GlobalPara::output::start_output_field
		) {
		// ��Ļ
		s.b_print = true;
		// ����
		s.b_writeTecplot = true;
	}
	if (istep % GlobalPara::output::step_per_output_hist == offset) {
		// �в�
		s.b_writeHist = true;
	}

	// ��һ��ǿ�����
	if (istep == 1) {
		s.b_print = true;
		s.b_writeTecplot = true;
		s.b_writeHist = true;
	}
}

void U2NITS::CDriver::modifySignalPack2_pause(SignalPack& s, int istep, int maxIteration, myfloat t, myfloat T, myfloat residualRho) {

	
	if (isnan(residualRho)) {// ��ɢ����ֹ
		s.b_writeContinue = true;
		s.b_writeTecplot = true;
		s.b_nanDetected = true;
		s.b_writeRecovery = false;// ��ֹԭ������recovery������nan��recovery����
		s.pauseSignal = _NanDetected;
	}
	else if (_kbhit()) {// ��esc��ֹ
		if (_getch() == 27) {
			s.b_writeContinue = true;
			s.pauseSignal = _EscPressed;
			std::cout << "Stopping at step " << istep << "...\n";
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

}

void U2NITS::CDriver::onSignalPack(const SignalPack& sp) {

}

