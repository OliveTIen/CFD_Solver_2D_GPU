#include "CDriver.h"
#include "../FVM_2D.h"
#include "../math/Common.h"
#include "../global/CExit.h"
#include "RoutineController.h"

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
	// 部分控制参数的引用。需要在readConfig后使用，否则未初始化
	const int istep_start = GlobalPara::time::istep_previous;// 已计算步
	const int maxIteration = GlobalPara::output::maxIteration;// 目标步
	int istep = istep_start;
	double t = GlobalPara::time::t_previous;// 当前物理时间
	const double T = GlobalPara::time::max_physical_time;

	// 迭代
	out.hist.setFilePath(out.getDir() + GlobalPara::basic::filename + "_hist.dat");
	out.set_tecplot_hist_path(out.getDir() + GlobalPara::basic::filename + "_hist.dat");
	out.timer.start();// 计时
	out.console.setClearStartPosition();
	while (true) {// 循环终止条件见signal
		istep++;
		// 更新文件名
		out.updateFileName(istep);

		solver.iteration(t, T);
		std::cout << ".";// 输出“.”表示正在计算，没有卡顿

		SignalPack sp = emitSignalPack(istep, maxIteration, t, T, solver.residualVector[0]);

		// 打印操作
		if (sp.b_print) {
			solver.updateResidualVector();
			// 调整CFL数，后面计算dt时会用到
			RoutineController::getInstance()->applyStrategy(solver.residualVector[0], GlobalPara::time::CFL, istep);
			out.console.clear();
			out.console.drawProgressBar((double)istep, (double)maxIteration, 45);
			out.console.setSolveInfo(istep_start, istep, maxIteration, out.getFieldWriter()->getNumTecplotFileWritten(), out.timer.getIime(), t, T, solver.residualVector, GlobalPara::time::CFL);
			out.console.print(out.console.getSolveInfo());
		}
		// 写文件操作
		if (sp.b_writeContinue || sp.b_writeTecplot || sp.b_writeRecovery || sp.b_writeHist) {
			solver.updateOutputNodeField();
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
		updateOldData(istep, t);
	}

	// 释放资源
	solver.freeMemory();

	// 停止
	CExit::pressAnyKeyToExit();

	/*
步骤：（优化后）


迭代开始
	计算（流场和残差） CSolver.iterate() 残差仅需输出时计算
	更新输出信息（包括文件名、下一次光标位置、输出进度）COutput.update()调用(COutput.updateData COutput.printScreen())
	输出文件（流场、残差）COutput.writeFile() 包括COutput.writeField() COutput.writeHist()(调用CSolver.updateResidual() )
	判断是否跳出循环（终止）this->checkStop() 消息队列？
	输出结束信息（提示词）COutput.updateScreen()
释放资源 CSolver.finalize()

	*/

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
		solver.node_host, solver.element_host, solver.element_U_old
	);
	exit(_Code);
}


U2NITS::CDriver::SignalPack U2NITS::CDriver::emitSignalPack(int istep, int maxIteration, myfloat t, myfloat T, myfloat residualRho) {
	SignalPack s;

	// 发射信号
	const int offset = 0;
	const double big_value = 1e8;

	if (istep % GlobalPara::output::step_per_print == offset) {
		// 屏幕
		s.b_print = true;
		// 还原点
		s.b_writeRecovery = true;
	}
	if (istep % GlobalPara::output::step_per_output_field == offset) {
		// 屏幕
		s.b_print = true;
		// 流场
		s.b_writeTecplot = true;
	}
	if (istep % GlobalPara::output::step_per_output_hist == offset) {
		// 残差
		s.b_writeHist = true;
	}

	if (residualRho > big_value) {// 残差过大，终止
		s.b_writeContinue = true;
		s.b_writeTecplot = true;
		s.b_nanDetected = true;
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

	return s;
}

void U2NITS::CDriver::onSignalPack(const SignalPack& sp) {

}

