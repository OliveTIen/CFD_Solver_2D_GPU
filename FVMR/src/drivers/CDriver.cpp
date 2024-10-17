#include "CDriver.h"
#include "../legacy/FVM_2D.h"
#include "../math/CommonValue.h"
#include "../global/CExit.h"
#include "RoutineController.h"
#include "../time/CIteration.h"
#include <chrono>
#include <sstream>

U2NITS::CDriver* U2NITS::CDriver::pCDriver = nullptr;

U2NITS::CDriver* U2NITS::CDriver::getInstance() {
	if (pCDriver == nullptr) {
		pCDriver = new CDriver;
	}
	return pCDriver;
}


void U2NITS::CDriver::start() {
	/*
	The bottle neck of GPU code is memcpy, which takes up about 50% of 
	the total calculation time.

	GPU����ƿ�����ڿ���ռ��50%ʱ��
	��Ҫÿ����device_to_host��ֻ����Ҫд�ļ�ʱ��device_to_host
	�����Ҫ��device_to_host����solver.updateResidualVector();֮ǰ
	*/

	out.console.printWelcomeInterface();
	input.readConfig();
	input.readField_1();
	solver.allocateMemory();
	solver.initializeData();
	// should be put after reading config ���ֿ��Ʋ��������á���Ҫ��readConfig��ʹ�ã�����δ��ʼ��
	const int istep_start = GlobalPara::time::istep_previous;
	const int maxIteration = GlobalPara::output::maxIteration;
	int istep = istep_start;// current step
	myfloat t = GlobalPara::time::t_previous;// current physical time
	const myfloat T = GlobalPara::time::max_physical_time;
	double current_time = 0.0;// program running time (in seconds)��ǰʱ��[��](�Ӽ�ʱ����ʼ����)
	double pure_calculate_speed = 0.0;// time calculation costs (in seconds) ���������㲿������ʱ��
	double mean_speed = 0.0;

	gui.init();
	// iteration
	out.hist.setFilePath(out.getDir() + GlobalPara::basic::filename + "_hist.dat");
	out.set_tecplot_hist_path(out.getDir() + GlobalPara::basic::filename + "_hist.dat");
	out.timer.start();
	current_time = out.timer.getTime();
	out.console.setClearStartPosition();
	while (true) {
		gui.render();
		getLastCudaError("gui.render() failed.");

		istep++;

		out.updateFileName(istep);

		current_time = out.timer.getTime();
		auto iteration_start = std::chrono::system_clock::now();

		if (GlobalPara::basic::useGPU) {
			CIteration::iteration_device_20240517(t, T, &solver);
		}
		else {
			CIteration::iteration_host(t, T, &solver);
		}

		auto iteration_end = std::chrono::system_clock::now();
		auto duration = std::chrono::duration_cast<std::chrono::microseconds>(iteration_end - iteration_start);
		double duration_seconds = double(duration.count()) * std::chrono::microseconds::period::num / std::chrono::microseconds::period::den;
		pure_calculate_speed = 1.0 / duration_seconds;

		if (!solver.isIterationStarted()) {
			solver.setIterationStarted(true);
		}

		if(enable_print) std::cout << ".";// to show that it doesn't stucks �����.����ʾ���ڼ��㣬û�п���

		double new_time = out.timer.getTime();
		//pure_calculate_speed = 1.0 / (new_time - current_time);
		current_time = new_time;
		mean_speed = out.console.getMeanSpeed(istep, istep_start, current_time);

		SignalPack sp;
		modifySignalPack1_output(sp, istep, maxIteration, current_time);
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

		// print ��ӡ����
		if (enable_print && sp.b_print) {
			solver.updateResidualVector();// ����elementField_host����residualVector
			// ����CFL�����������dtʱ���õ�
			bool apply_strategy_success = false;
			RoutineController::getInstance()->applyStrategy(solver.residualVector[0], GlobalPara::time::CFL, istep, apply_strategy_success);
			if (apply_strategy_success) {
				std::stringstream ss;
				ss << "step = " << istep << "\tresidual_rho = " << solver.residualVector[0] << "\tCFL = " << GlobalPara::time::CFL 
					<< "\tmean speed = "<< mean_speed << " step/s"
					<< "\tpure calculate speed = "<< pure_calculate_speed << " step/s" 
					<< "\n";
				LogWriter::log(ss.str());
			}
			out.console.clear();
			out.console.drawProgressBar((double)istep, (double)maxIteration, 45);
			out.console.setSolveInfo(istep, maxIteration, out.getFieldWriter()->getNumTecplotFileWritten(), current_time, t, T, solver.residualVector, GlobalPara::time::CFL, mean_speed, pure_calculate_speed);
			out.console.print(out.console.getSolveInfo());
		}
		// show info ��¼����ٶ�
		if (sp.b_writeContinue || sp.b_writeTecplot) {
			out.log.log("write field file. \n");
			std::stringstream ss;
			ss << "step = " << istep << "\tresidual_rho = " << solver.residualVector[0] << "\tCFL = " << GlobalPara::time::CFL
				<< "\tmean speed = " << mean_speed << " step/s"
				<< "\tpure calculate speed = " << pure_calculate_speed << " step/s"
				<< "\n";
			LogWriter::log(ss.str());

		}


		// write file
		if (enable_write_file && (sp.b_writeContinue || sp.b_writeTecplot || sp.b_writeRecovery || sp.b_writeHist)) {
			
			solver.updateOutputNodeField();// calculate nodeField_host by elementField_host
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
		// jump out of the loop
		if (sp.pauseSignal != _NoSignal) {
			out.log.log(out.console.getSolveInfo());
			out.log.logAndPrintSignalInfo(sp.pauseSignal);// show prompt words after stopping
			break;
		}
		updateOldData_istep_t(istep, t);
	}
	while (!gui.windowShouldClose()) {
		gui.render();
	}
	gui.cleanup();
	solver.freeMemory();
	CExit::pressAnyKeyToExit();

}

void U2NITS::CDriver::saveAndExit(int _Code) {
	getInstance()->m_saveAndExit(_Code);
}

void U2NITS::CDriver::m_saveAndExit(int _Code) {
	if (!solver.isIterationStarted()) {
		exit(_Code);
	}

	out.updateFileName(m_last_istep);
	out.getFieldWriter()->writeContinueFile_1(
		m_last_istep, m_last_t, out.continueFilePath,
		solver.node_host, solver.element_host, solver.elementField_host.Uold
	);
	exit(_Code);
}

void U2NITS::CDriver::modifySignalPack1_output(SignalPack& s, int istep, int maxIteration, double current_time) {

	// force to write file in the first step ��һ��ǿ�����
	if (istep == 1) {
		s.b_print = true;
		s.b_writeTecplot = true;
		s.b_writeHist = true;
	}

	// write file every certain steps ÿ��һ���������
	const int offset = 0;
	if (istep % GlobalPara::output::step_per_print == offset) {
		s.b_print = true;
	}
	if (istep % GlobalPara::output::step_per_output_field == offset &&
		istep >= GlobalPara::output::start_output_field
		) {
		s.b_print = true;
		s.b_writeTecplot = true;
	}
	if (istep % GlobalPara::output::step_per_output_hist == offset) {
		s.b_writeHist = true;
	}

	// write file every certain time ÿ��һ��ʱ�����
	double duration_since_last_snapshot = current_time - m_last_snapshot_time;
	if (duration_since_last_snapshot > m_snapshot_cold_time) {
		s.b_writeContinue = true;
		m_last_snapshot_time = current_time;
	}
}

void U2NITS::CDriver::modifySignalPack2_pause(SignalPack& s, int istep, int maxIteration, myfloat t, myfloat T, myfloat residualRho) {

	
	if (isnan(residualRho)) {
		s.b_writeContinue = true;
		s.b_writeTecplot = true;
		s.b_nanDetected = true;
		s.b_writeRecovery = false;// do NOT write abnormal recovery file ��ֹԭ������recovery������nan��recovery����
		s.pauseSignal = _NanDetected;
	}
	else if ((_kbhit() && _getch() == 27) || gui.windowShouldClose()) {// on "ESC" pressed
		s.b_writeContinue = true;
		s.pauseSignal = _EscPressed;
		std::cout << "Stopping at step " << istep << "...\n";
	}
	else if (t > T || istep >= maxIteration) {
		s.b_writeContinue = true;
		s.b_writeTecplot = true;
		s.pauseSignal = _TimeReached;
	}
	else if (residualRho <= U2NITS::Math::EPSILON) {
		s.b_writeContinue = true;
		s.b_writeTecplot = true;
		s.pauseSignal = _StableReached;
	}

}

