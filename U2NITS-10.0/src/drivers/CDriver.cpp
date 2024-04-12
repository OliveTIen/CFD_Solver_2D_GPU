#include "CDriver.h"
#include "../FVM_2D.h"
#include "../solvers/GPUSolver2.h"
#include "../input/CInput.h"
#include "../output/COutput.h"


void U2NITS::CDriver::run_current() {

/*
[开发中]以下为GPU代码

*/
	
	CInput input;
	COutput out;
	GPU::GPUSolver2 solver;

	out.console.printWelcome_2();
	input.readConfig();
	input.readField_old();

	solver.initialize();
	// 部分控制参数的引用。需要在readConfig后使用，否则未初始化
	const int istep_previous = GlobalPara::time::istep_previous;// 已计算步
	const int maxIteration = GlobalPara::output::maxIteration;// 目标步
	double t = GlobalPara::time::t_previous;// 当前物理时间
	const double T = GlobalPara::time::max_physical_time;

	// 迭代
	out.hist.setFilePath(out.getDir() + GlobalPara::basic::filename + "_hist.dat");
	out.timer.start();// 计时
	out.console.setClearStartPosition();
	for (int istep = istep_previous; ; istep++) {// 循环终止条件见signal
		// 更新文件名
		out.updateFileName(istep);
		solver.iteration(t, T);
		std::cout << ".";// 输出“.”表示正在计算，没有卡顿

		SignalPack sp = emitSignalPack(istep, maxIteration, t, T, solver.residualVector[0]);

		// 打印操作
		if (sp.b_print) {
			solver.getResidual();
			out.console.clear();
			out.console.drawProgressBar(int(double(istep) / double(maxIteration) * 100.0));

			double calculateTime = out.timer.getIime();
			double calculateSpeed = ((istep - istep_previous) / calculateTime);
			out.console.assemblySolveInfo(calculateTime, istep, maxIteration, calculateSpeed, out.field.getNumTecplotFileWritten(), t, T, solver.residualVector);
			out.console.print(out.console.getSolveInfo());
		}
		// 写文件操作
		if (sp.b_writeContinue || sp.b_writeTecplot) {
			solver.updateOutputNodeField();// 根据单元U更新节点ruvp
			if (sp.b_nanDetected) {
				out.console.printInfo(out.console.InfoType::type_nan_detected);
				out.continueFilePath = out.continueFilePath_nan;
			}
			if (sp.b_writeContinue) {
				out.field.writeContinueFile_GPU(
					istep, t, out.continueFilePath,
					solver.node_host, solver.element_host, solver.elementField_host
				);
			}
			if (sp.b_writeTecplot) {
				out.field.writeTecplotFile_GPU(t, out.tecplotFilePath, "title", solver.node_host, solver.element_host, solver.outputNodeField);
			}
		}
		if (sp.b_writeHist) {
			// 获取残差，结果保存在residual_vector中
			solver.getResidual();
			out.hist.writeHistFile_2(istep, solver.residualVector, solver.residualVectorSize);
		}
		// 终止操作
		if (sp.pauseSignal != _NoSignal) {
			out.log.log(out.console.getSolveInfo());
			out.log.logAndPrintSignalInfo(sp.pauseSignal);// 提示词 放在写文件之后，防止卡顿
			break;// 不可放进单独函数中，因为要终止循环
		}
	}

	// 释放资源
	solver.freeMemory();

	// 停止
#ifdef _WIN32
	system("pause");
#elif defined __linux__
	getchar();
#endif

	/*
步骤：（优化后）
读取配置（控制参数） CInput.readConfig()
初始化流场：读取续算文件/读取网格文件+初始化 CInput.readField()
输出初始信息（程序信息、边界参数等）COutput.printScreen()
输出残差文件头 COutput.writeHistFileHead() 该判断放到函数中去
计时器开始 CTimer.start()
初始化资源（GPU内存）CSolver.initialize()
迭代开始
	计算（流场和残差） CSolver.iterate() 残差仅需输出时计算
	更新输出信息（包括文件名、下一次光标位置、输出进度）COutput.update()调用(COutput.updateData COutput.printScreen())
	输出文件（流场、残差）COutput.writeFile() 包括COutput.writeField() COutput.writeHist()(调用CSolver.updateResidual() )
	判断是否跳出循环（终止）this->checkStop() 消息队列？
	输出结束信息（提示词）COutput.updateScreen()
释放资源 CSolver.finalize()

为什么要多线程：
磁盘IO操作比较耗时，此时CPU空闲，可以提高效率
一次IO操作耗时相当于2~3次迭代；网格量增加时，IO和计算时间都会增加
掌握多线程知识
用消息传递的方式输出信息

	*/

}

U2NITS::CDriver::SignalPack U2NITS::CDriver::emitSignalPack(int istep, int maxIteration, real t, real T, real residualRho) {
	SignalPack s;

	// 发射信号
	const int offset = 0;
	const double big_value = 1e8;

	//PauseSignal pauseSignal = _NoSignal;// 暂停信号 提示词分类
	if (istep % GlobalPara::output::step_per_print == offset) {// 定期输出进度
		s.b_print = true;
	}
	if (istep % GlobalPara::output::step_per_output_field == offset) {// 定期输出流场
		s.b_writeTecplot = true;
		s.b_print = true;
	}
	if (istep % GlobalPara::output::step_per_output_hist == offset) {// 定期输出残差
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
			std::cout << "Stopping...\n";
		}
	}
	else if (t > T || istep >= maxIteration) {// 达到规定迭代时间或迭代步数，终止
		s.b_writeContinue = true;
		s.b_writeTecplot = true;
		s.pauseSignal = _TimeReached;
	}
	else if (residualRho <= GlobalPara::constant::epsilon) {// 残差足够小，终止
		s.b_writeContinue = true;
		s.b_writeTecplot = true;
		s.pauseSignal = _StableReached;
	}

	return s;
}

void U2NITS::CDriver::onSignalPack(const SignalPack& sp) {

}
