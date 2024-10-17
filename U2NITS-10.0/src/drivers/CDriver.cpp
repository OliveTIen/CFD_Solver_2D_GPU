#include "CDriver.h"
#include "../FVM_2D.h"
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
	GPU代码瓶颈在于拷贝占用50%时间
	不要每步都device_to_host，只有需要写文件时才device_to_host
	因此需要把device_to_host放在solver.updateResidualVector();之前
	*/

	out.console.printWelcomeInterface();
	input.readConfig();
	input.readField_1();
	solver.allocateMemory();
	solver.initializeData_byOldData();
	// 部分控制参数的引用。需要在readConfig后使用，否则未初始化
	const int istep_start = GlobalPara::time::istep_previous;// 起始步
	const int maxIteration = GlobalPara::output::maxIteration;// 目标步
	int istep = istep_start;// 当前步
	myfloat t = GlobalPara::time::t_previous;// 当前物理时间
	const myfloat T = GlobalPara::time::max_physical_time;// 目标物理时间
	double current_time = 0.0;// 当前时间[秒](从计时器开始算起)
	double pure_calculate_speed = 0.0;// 迭代纯计算部分所用时间[s]
	double mean_speed = 0.0;

	// 迭代
	out.hist.setFilePath(out.getDir() + GlobalPara::basic::filename + "_hist.dat");
	out.set_tecplot_hist_path(out.getDir() + GlobalPara::basic::filename + "_hist.dat");
	out.timer.start();// 计时
	current_time = out.timer.getTime();
	out.console.setClearStartPosition();
	while (true) {// 循环终止条件见signal
		istep++;
		// 更新文件名
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

		std::cout << ".";// 输出“.”表示正在计算，没有卡顿

		double new_time = out.timer.getTime();
		//pure_calculate_speed = 1.0 / (new_time - current_time);
		current_time = new_time;
		mean_speed = out.console.getMeanSpeed(istep, istep_start, current_time);

		SignalPack sp;
		modifySignalPack1_output(sp, istep, maxIteration, current_time);
		modifySignalPack2_pause(sp, istep, maxIteration, t, T, solver.residualVector[0]);

		/*
		新代码 若需要打印或输出文件，则拷贝设备到主机
		需要拷贝的数据：
		elementField_host 计算残差residualVector、计算节点输出值nodeField_host(用于tecplot)、输出continuefile/recoveryfile
		*/ 
		if (GlobalPara::basic::useGPU) {
			if (sp.b_print ||
				sp.b_writeContinue || sp.b_writeTecplot || sp.b_writeRecovery || sp.b_writeHist) {
				GPU::ElementFieldSoA::cuda_memcpy(&solver.elementField_host, &solver.elementField_device, cudaMemcpyDeviceToHost);
			}
		}

		// 打印操作
		if (sp.b_print) {
			solver.updateResidualVector();// 根据elementField_host计算residualVector
			// 调整CFL数，后面计算dt时会用到
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
		// 记录求解速度
		if (sp.b_writeContinue || sp.b_writeTecplot) {
			out.log.log("write field file. \n");
			std::stringstream ss;
			ss << "step = " << istep << "\tresidual_rho = " << solver.residualVector[0] << "\tCFL = " << GlobalPara::time::CFL
				<< "\tmean speed = " << mean_speed << " step/s"
				<< "\tpure calculate speed = " << pure_calculate_speed << " step/s"
				<< "\n";
			LogWriter::log(ss.str());

		}


		// 写文件操作
		if (sp.b_writeContinue || sp.b_writeTecplot || sp.b_writeRecovery || sp.b_writeHist) {
			
			solver.updateOutputNodeField();// 根据elementField_host计算nodeField_host
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
		// 终止操作
		if (sp.pauseSignal != _NoSignal) {
			out.log.log(out.console.getSolveInfo());
			out.log.logAndPrintSignalInfo(sp.pauseSignal);// 提示词 放在写文件之后，防止卡顿
			break;// 不可放进单独函数中，因为要终止循环
		}
		updateOldData_istep_t(istep, t);
	}

	// 释放资源
	solver.freeMemory();

	// 停止
	CExit::pressAnyKeyToExit();

}

void U2NITS::CDriver::saveAndExit(int _Code) {
	getInstance()->m_saveAndExit(_Code);
}

void U2NITS::CDriver::m_saveAndExit(int _Code) {
	if (!solver.isIterationStarted()) {// 迭代未开始，不用输出文件
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

	// 第一步强制输出
	if (istep == 1) {
		s.b_print = true;
		s.b_writeTecplot = true;
		s.b_writeHist = true;
	}

	// 每隔一定步数输出
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

	// 每隔一段时间输出
	double duration_since_last_snapshot = current_time - m_last_snapshot_time;
	if (duration_since_last_snapshot > m_snapshot_cold_time) {
		s.b_writeContinue = true;// 每10分钟保存一次快照
		m_last_snapshot_time = current_time;
	}
}

void U2NITS::CDriver::modifySignalPack2_pause(SignalPack& s, int istep, int maxIteration, myfloat t, myfloat T, myfloat residualRho) {

	
	if (isnan(residualRho)) {// 发散，终止
		s.b_writeContinue = true;
		s.b_writeTecplot = true;
		s.b_nanDetected = true;
		s.b_writeRecovery = false;// 防止原正常的recovery被含有nan的recovery覆盖
		s.pauseSignal = _NanDetected;
	}
	else if (_kbhit()) {// 按esc终止
		if (_getch() == 27) {
			s.b_writeContinue = true;
			s.pauseSignal = _EscPressed;
			std::cout << "Stopping at step " << istep << "...\n";
		}
	}
	else if (t > T || istep >= maxIteration) {// 达到规定迭代时间或迭代步数，终止
		s.b_writeContinue = true;
		s.b_writeTecplot = true;
		s.pauseSignal = _TimeReached;
	}
	else if (residualRho <= U2NITS::Math::EPSILON) {// 残差足够小，终止
		s.b_writeContinue = true;
		s.b_writeTecplot = true;
		s.pauseSignal = _StableReached;
	}

}

