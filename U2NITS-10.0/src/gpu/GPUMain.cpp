#include "GPUMain.h"
#include "../GlobalPara.h"
#include "../FVM_2D.h"
#include "../global/FilePathManager.h"
#include <ctime>
#include "../output/ConsolePrinter.h"
#include "../output/FieldWriter.h"
#include "../output/HistWriter.h"
#include "GPUSolver2.h"
#include "../output/ResidualCalculator.h"
#include "../output/LogWriter.h"
#include "../input/SU2MeshReader.h"
#include "../input/InpMeshReader.h"
#include "../global/StringProcessor.h"

void GPU::GPUMain::run_GPU() {
	FVM_2D* pFVM2D = FVM_2D::pFVM2D;
	//// 初始化
	// 尝试读取续算文件。若读取失败或者_continue==false则从头开始
	bool startFromZero = false;
	if (GlobalPara::basic::_continue) {
		int flag_readContinue = pFVM2D->readContinueFile();
		if (flag_readContinue == -1) {
			std::string error_msg = "Info: Fail to read previous mesh(pause_*.dat). ";
			error_msg += "Will try to start from zero again.\n";
			LogWriter::writeLogAndCout(error_msg, LogWriter::Info);
			startFromZero = true;
			// 防止后面histWriter不写文件头
			GlobalPara::basic::_continue = false;
		}
	}
	else {
		startFromZero = true;
	}
	// 从头开始。读取网格、初始化
	if (startFromZero) {
		const std::string& type = GlobalPara::basic::meshFileType;
		if (type == "inp") {
			std::string dir = FilePathManager::getInstance()->getInputDirectory();
			int flag_readMesh = InpMeshReader::readGmeshFile(dir + GlobalPara::basic::filename + ".inp");
			if (flag_readMesh == -1) {
				std::string error_msg = "Error: Fail to read mesh. Program will exit.\n";
				LogWriter::writeLogAndCout(error_msg);
				return;//退出
			}
			pFVM2D->setInitialCondition();
		}
		else if (type == "su2") {
			std::string dir = FilePathManager::getInstance()->getInputDirectory();
			int flag_readMesh = SU2MeshReader::readFile(dir + GlobalPara::basic::filename + ".su2", true);
			if (flag_readMesh == -1) {
				std::string error_msg = "Error: Fail to read mesh. Program will exit.\n";
				LogWriter::writeLogAndCout(error_msg);
				return;//退出
			}
			pFVM2D->setInitialCondition();
		}
		else {
			std::string error_msg = "Invalid mesh file type: " + type + ". Program will exit.\n";
			LogWriter::writeLogAndCout(error_msg);
			return;
		}

	}

	// 日志记录边界参数
	std::string str;
	str += "BoundaryCondition:\n";
	str += "inlet::ruvp\t" + StringProcessor::doubleArray_2_string(GlobalPara::boundaryCondition::_2D::inlet::ruvp, 4)
		+ "\noutlet::ruvp\t" + StringProcessor::doubleArray_2_string(GlobalPara::boundaryCondition::_2D::outlet::ruvp, 4)
		+ "\ninf::ruvp\t" + StringProcessor::doubleArray_2_string(GlobalPara::boundaryCondition::_2D::inf::ruvp, 4)
		+ "\n";
	LogWriter::writeLog(str);

	//// 求解
	if (GlobalPara::physicsModel::equation == _EQ_euler) {
		solve_GPU();
	}
	else {
		std::cout << "Error: Invalid equation type";
	}

}

void GPU::GPUMain::solve_GPU() {
	/*
	[开发中]以下为GPU代码
	TODO:
	*/
	FVM_2D* pFVM2D = FVM_2D::pFVM2D;
	const int residual_vector_size = 4;
	const int istep_previous = GlobalPara::time::istep_previous;
	const double previous_t = GlobalPara::time::t_previous;
	const double big_value = 1e8;
	const int maxIteration = GlobalPara::output::maxIteration;
	const double T = GlobalPara::time::T;//目标物理时间。用于控制是否终止计算
	int nFiles = 0;//输出文件个数，用于显示
	int iCurrentAutosaveFile = 0;// 自动保存文件的指标，构成循环队列
	double residual_vector[residual_vector_size]{ 1,1,1,1 };
	double t = previous_t;//当前物理时间。由续算时的起始时间决定

	std::string solveInfo;
	std::string outputPathWithSlash = FilePathManager::getInstance()->getOutputDirectory();
	std::string basicFileName = GlobalPara::basic::filename;
	clock_t start_t;
	COORD p1;
	GPU::GPUSolver2 gpuSolver2;
	HistWriter histWriter(outputPathWithSlash + basicFileName + "_hist.dat");
	FieldWriter tecplotWriter;
	FieldWriter continueWriter;

	// 初始化
	start_t = clock();
	histWriter.writeHistFileHead();
	gpuSolver2.initialze();
	p1 = ConsolePrinter::getCursorPosition();
	ConsolePrinter::drawProgressBar(int(t / T));

	// 迭代
	for (int istep = istep_previous; istep <= maxIteration && t <= T; istep++) {
		// 更新文件名
		char szBuffer[20];
		sprintf_s(szBuffer, _countof(szBuffer), "%04d", istep);//若istep小于4位数，则补0
		std::string str_istep_withBracket = "[" + std::string(szBuffer) + "]";
		std::string tecplotFilePath = outputPathWithSlash + basicFileName + str_istep_withBracket + ".dat";
		std::string continueFilePath = outputPathWithSlash + "pause_" + basicFileName + str_istep_withBracket + ".dat";
		std::string continueFilePath_nan = outputPathWithSlash + "nan_" + basicFileName + str_istep_withBracket + ".dat";

		gpuSolver2.iteration(t, T);

		// --- 输出 ---
		std::cout << ".";// 输出“.”表示正在计算，没有卡顿
		bool b_writeContinue = false;
		bool b_writeTecplot = false;
		bool b_pause = false;//暂停信号，为true则终止
		MySignal promptSignal = _NoSignal;// 提示词分类
		// 定期输出进度
		if (istep % GlobalPara::output::step_per_print == 1) {
			ConsolePrinter::clearDisplay(p1, ConsolePrinter::getCursorPosition());
			ConsolePrinter::setCursorPosition(p1);
			ConsolePrinter::drawProgressBar(int(t / T * 100 + 2));//+x是为了防止停在99%

			double calculateTime = (double)(clock() - start_t) / CLOCKS_PER_SEC;
			double calculateSpeed = ((istep - istep_previous) / calculateTime);
			solveInfo = ConsolePrinter::assemblySolveInfo(calculateTime, istep, maxIteration, calculateSpeed, nFiles, t, T, residual_vector);
			std::cout << solveInfo;
		}
		// 发散则立刻跳出循环
		if (residual_vector[0] > big_value) {
			// 这里直接break，因此不需要用到后面的b_writeContinue、b_writeTecplot
			// 况且writeContinueFile的参数不同，有nan
			ConsolePrinter::printInfo(ConsolePrinter::InfoType::type_nan_detected);
			//writeTecplotFile(tecplotFilePath, t);// 别删
			// 根据单元U更新节点ruvp
			gpuSolver2.updateOutputNodeField();
			tecplotWriter.setFilePath(tecplotFilePath);
			tecplotWriter.writeTecplotFile_GPU(t, "title", gpuSolver2.node_host, gpuSolver2.element_host, gpuSolver2.outputNodeField);// 完成于2024-02-28
			// 写nan文件
			continueWriter.writeContinueFile_GPU(
				istep, t, continueFilePath_nan,
				gpuSolver2.node_host, gpuSolver2.element_host, gpuSolver2.elementField_host, &(pFVM2D->boundaryManager)
			);// 完成于2024-02-28
			break;
		}
		// 定期输出流场
		if (istep % GlobalPara::output::step_per_output == 1) {
			//writeTecplotFile(tecplotFilePath, t);
			b_writeTecplot = true;
			nFiles++;
			//.autosave 自动保存文件 续算文件
			//最多保存global::autosaveFileNum个autosave文件，再多则覆写之前的
			std::string autosaveFilePath = outputPathWithSlash + "autosave" + std::to_string(iCurrentAutosaveFile) + "_" + GlobalPara::basic::filename + ".dat";
			continueWriter.writeContinueFile_GPU(
				istep, t, autosaveFilePath,
				gpuSolver2.node_host, gpuSolver2.element_host, gpuSolver2.elementField_host, &(pFVM2D->boundaryManager)
			);// 完成于2024-02-28
			iCurrentAutosaveFile++;
			if (iCurrentAutosaveFile >= GlobalPara::output::autosaveFileNum)iCurrentAutosaveFile = 0;

		}
		// 定期输出残差
		if (istep % GlobalPara::output::step_per_output_hist == 1) {
			//计算残差，结果保存在residual_vector中
			ResidualCalculator::cal_residual_GPU(gpuSolver2.element_U_old, gpuSolver2.elementField_host, ResidualCalculator::NORM_INF, residual_vector);
			histWriter.writeHistFileData(istep, residual_vector, residual_vector_size);
		}

		// --- 终止 ---
		// 按esc终止
		if (_kbhit()) {
			if (_getch() == 27) {
				b_writeContinue = true;
				b_pause = true;
				promptSignal = _EscPressed;
			}
		}
		// 达到规定迭代时间或迭代步数，终止
		else if (t >= T || istep >= maxIteration) {
			b_writeContinue = true;
			b_writeTecplot = true;
			b_pause = true;
			promptSignal = _TimeReached;
		}
		// 残差足够小，终止
		else if (residual_vector[0] <= Constant::epsilon) {
			b_writeContinue = true;
			b_writeTecplot = true;
			b_pause = true;
			promptSignal = _StableReached;
		}

		// 写文件操作
		if (b_writeContinue || b_writeTecplot) {
			// 根据单元U更新节点ruvp
			gpuSolver2.updateOutputNodeField();
			if (b_writeContinue) {
				continueWriter.writeContinueFile_GPU(
					istep, t, continueFilePath,
					gpuSolver2.node_host, gpuSolver2.element_host, gpuSolver2.elementField_host, &(pFVM2D->boundaryManager)
				);
			}
			if (b_writeTecplot) {
				tecplotWriter.setFilePath(tecplotFilePath);
				tecplotWriter.writeTecplotFile_GPU(t, "title", gpuSolver2.node_host, gpuSolver2.element_host, gpuSolver2.outputNodeField);
			}
		}
		// 提示词 放在写文件之后，防止卡顿
		printSignalInfo(promptSignal);
		// 终止操作
		if (b_pause) {
			LogWriter::writeLog(solveInfo);
			break;
		}
	}

	// 释放资源
	gpuSolver2.finalize();
}

void GPU::GPUMain::printSignalInfo(MySignal signal) {
	switch (signal) {
	case _NoSignal:
		break;
	case _EscPressed:
		std::cout << "[ESC]: Computation Interrupted\n";
		break;
	case _TimeReached:
		std::cout << "Computation finished\n";
		break;
	case _StableReached:
		std::cout << "Computation finished as the field is already stable\n";
		break;
	default:
		break;
	}
}
