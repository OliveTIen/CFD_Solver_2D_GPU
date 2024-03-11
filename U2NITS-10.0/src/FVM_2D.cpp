#include "FVM_2D.h"
#include "output/ConsolePrinter.h"
#include "output/ResidualCalculator.h"
#include "output/LogWriter.h"
#include "global/StringProcessor.h"
#include "global/SystemInfo.h"
#include "global/FilePathManager.h"
#include "output/HistWriter.h"
#include "input/InpMeshReader.h"
#include "global/VectorProcessor.h"
#include "input/SU2MeshReader.h"
#include "solvers/GPUSolver2.h"
#include "space/FluxGPU.h"
#include "time/CalculateDt.h"
#include "output/FieldWriter.h"

FVM_2D* FVM_2D::pFVM2D = nullptr;

FVM_2D* FVM_2D::getInstance() {
	if (pFVM2D == nullptr) {
		pFVM2D = new FVM_2D();
	}
	return pFVM2D;
}

void FVM_2D::run() {
	// run_CPU

	//// 初始化
	// 尝试读取续算文件。若读取失败或者_continue==false则从头开始
	bool startFromZero = false;
	if (GlobalPara::basic::_continue) {
		int flag_readContinue = readContinueFile();
		if (flag_readContinue == -1) {
			std::string error_msg = "Warning: Fail to read previous mesh(pause_*.dat). ";
			error_msg += "Will try to start from zero again.\n";
			LogWriter::writeLogAndCout(error_msg, LogWriter::Error);
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
				LogWriter::writeLogAndCout(error_msg, LogWriter::Error);
				return;//退出
			}
			setInitialCondition();
		}
		else if (type == "su2") {
			std::string dir = FilePathManager::getInstance()->getInputDirectory();
			int flag_readMesh = SU2MeshReader::readFile(dir + GlobalPara::basic::filename + ".su2", true);
			if (flag_readMesh == -1) {
				std::string error_msg = "Error: Fail to read mesh. Program will exit.\n";
				LogWriter::writeLogAndCout(error_msg, LogWriter::Error);
				return;//退出
			}
			setInitialCondition();
		}
		else {
			std::string error_msg = "Invalid mesh file type: " + type + ". Program will exit.\n";
			LogWriter::writeLogAndCout(error_msg, LogWriter::Error);
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
		solve_CPU2("_Euler.out", "_Euler.info");
	}
	else {
		std::cout << "Error: Invalid equation type";
	}

}

void FVM_2D::setInitialCondition() {
	switch (GlobalPara::initialCondition::type) {
	case 1:
	{
		//1.常数
		using namespace GlobalPara::boundaryCondition::_2D;
		for (int ie = 0; ie < elements.size(); ie++) {
			Math_2D::ruvp_2_U(inf::ruvp, elements[ie].U, GlobalPara::constant::gamma);
			//if(elements[ie].ID==173)system()
		}
		std::string str;
		str += "InitialCondition:\nUniform flow, ruvp\t" + StringProcessor::doubleArray_2_string(inf::ruvp, 4) + "\n";
		LogWriter::writeLog(str);
	}
		break;

	case 2:
	{
		//2.均匀流和等熵涡的叠加

		using namespace GlobalPara::boundaryCondition::_2D;
		std::string str;
		str += "InitialCondition:\nUniform flow + isentropicVortex, ruvp of Uniform flow:" + StringProcessor::doubleArray_2_string(inf::ruvp, 4) + "\n";
		LogWriter::writeLogAndCout(str, LogWriter::Info);
		isentropicVortex_2(5, 5, 5, inf::ruvp);
		//for (int ie = 0; ie < elements.size(); ie++) {
		//	//均匀流+等熵涡 Isentropic vortex
		//	double deltau, deltav, deltaT;
		//	double x = elements[ie].x; double y = elements[ie].y;
		//	double xc = 5; double yc = 5;
		//	isentropicVortex(x, y, xc, yc, 1, deltau, deltav, deltaT);
		//	double ruvp[4]{};
		//	ruvp[1] = inf::ruvp[1] + deltau;
		//	ruvp[2] = inf::ruvp[2] + deltav;
		//	double rho0 = inf::ruvp[0];
		//	double p0 = inf::ruvp[3];
		//	double T0 = p0 / rho0 / GlobalPara::constant::gamma;//p=rho*R*T
		//	double deltap = deltaT * rho0 / (1 - Constant::R * rho0 * T0 / pow(p0, GlobalPara::constant::gamma));
		//	double deltarho = deltap * rho0 / pow(p0, GlobalPara::constant::gamma);
		//	ruvp[0] = rho0 + deltarho;
		//	ruvp[3] = p0 + deltap;
		//	Math_2D::ruvp_2_U(ruvp, elements[ie].U, GlobalPara::constant::gamma);
		//}

	}
		break;

	}

}

//void FVM_2D::solve_CPU(std::string suffix_out, std::string suffix_info) {
//	// 输入参数没有被用到
//
//	Solver_2D solver;
//	FilePathManager* filePathManager = FilePathManager::getInstance();
//	HistWriter histWriter(filePathManager->outputFolder_path + "\\" + GlobalPara::basic::filename + "_hist.dat");
//	const int residual_vector_size = 4;
//	const double big_value = 1e8;
//	double residual_vector[residual_vector_size]{1,1,1,1};
//	histWriter.writeHistFileHead();
//
//	//输出格式(tecplot)，在文件夹中输出，每一帧存成一个文件
//	std::vector<Element_2D> elements_old;
//	std::vector<double> u_nodes;// 节点u, size=nodes.size()。输出流场时，用作缓存。
//	std::vector<double> v_nodes;
//	std::vector<double> rho_nodes;
//	std::vector<double> p_nodes;
//
//	//变量
//	const int maxIteration = GlobalPara::output::maxIteration;
//	double t = GlobalPara::time::t_previous;//当前物理时间。
//	int istep_previous = GlobalPara::time::istep_previous;
//	const double T = GlobalPara::time::T;//目标物理时间。用于控制是否终止计算
//	int nFiles = 0;//输出文件个数，用于显示
//	int iCurrentAutosaveFile = 0;// 自动保存文件的指标，构成循环队列
//	bool signal_pause = 0;//暂停信号，用于控制
//	clock_t start_t;//计时变量，用于计算求解时间(基于CPU周期)
//	start_t = clock();
//	double calculateTime = 0.0;
//	double calculateSpeed = 0.0;
//	std::string solveInfo;
//	//等熵涡误差计算文件头
//	if (GlobalPara::initialCondition::type == 2)writeFileHeader_isentropicVortex();
//	COORD p1 = ConsolePrinter::getCursorPosition();
//	ConsolePrinter::drawProgressBar(int(t / T));	//绘制进度条
//	for (int istep = istep_previous; istep <= maxIteration && t <= T; istep++) {
//		elements_old = elements;
//		//计算时间
//		double dt = caldt(t, T);
//		t += dt;
//		//计算通量、时间推进
//		solver.evolve(dt);
//		//根据单元U更新节点ruvp
//		calculateNodeValue(rho_nodes, u_nodes, v_nodes, p_nodes);
//		//输出进度表示正在计算，没有卡顿
//		std::cout << ".";
//		//定期输出进度
//		if (istep % GlobalPara::output::step_per_print == 1) {
//			ConsolePrinter::clearDisplay(p1, ConsolePrinter::getCursorPosition());
//			ConsolePrinter::setCursorPosition(p1);
//			ConsolePrinter::drawProgressBar(int(t / T * 100 + 2));//+x是为了防止停在99%
//
//			calculateTime = (double)(clock() - start_t) / CLOCKS_PER_SEC;//
//			calculateSpeed = ((istep - istep_previous)/calculateTime);
//			solveInfo = ConsolePrinter::assemblySolveInfo(calculateTime, istep, maxIteration, calculateSpeed, nFiles, t, T, residual_vector);
//			std::cout << solveInfo;
//		}
//		//检查非法值
//		if (isNan()|| residual_vector[0] > big_value) {
//			std::cout << "\nWarning: \"NaN\" detected. ";
//			std::cout << "Possible Reason: \n"
//				<< "  1. Last time, this program terminated abnormally, leading to broken autosave files.\n"
//				<< "  2. Invalid boundary condition.\n"
//				<< "  3. When you continue to compute, you use a different boundary condition.\n";
//			std::cout << "Computation stopped.\n";
//			char szBuffer[20];
//			sprintf_s(szBuffer, _countof(szBuffer), "%04d", istep);
//			writeTecplotFile(filePathManager->getExePath_withSlash() + "output\\" + GlobalPara::basic::filename + "[" + szBuffer + "].dat",
//				t, rho_nodes, u_nodes, v_nodes, p_nodes);
//			writeContinueFile(filePathManager->getExePath_withSlash() + "output\\nan_" + GlobalPara::basic::filename + "[" + szBuffer + "].dat", t, istep);
//			break;
//		}
//		//定期输出流场
//		if (istep % GlobalPara::output::step_per_output_field == 1) {
//			//.dat 流场显示文件 tecplot格式
//			char szBuffer[20];
//			sprintf_s(szBuffer, _countof(szBuffer), "%04d", istep);//若istep小于4位数，则补0
//			writeTecplotFile(filePathManager->getExePath_withSlash() + "output\\" + GlobalPara::basic::filename + "[" + szBuffer + "].dat",
//				t, rho_nodes, u_nodes, v_nodes, p_nodes);
//			nFiles++;
//			//.autosave 自动保存文件 续算文件
//			updateAutoSaveFile(t, istep, iCurrentAutosaveFile);
//
//		}
//		//定期输出残差
//		if (istep % GlobalPara::output::step_per_output_hist == 1) {
//			//计算残差，结果保存在residual_vector中
//			ResidualCalculator::cal_residual(elements_old, elements, ResidualCalculator::NORM_INF, residual_vector);
//			histWriter.writeHistFileData(istep, residual_vector, residual_vector_size);
//
//			//等熵涡的误差文件 
//			if (GlobalPara::initialCondition::type == 2) {
//				//求解时间
//				//end_t = clock();
//				double time_used = (double)(clock() - start_t) / CLOCKS_PER_SEC;
//				//输出误差文件
//				ResidualCalculator::cal_error_isentropicVortex(0, 0, 10, 10, 5, t, istep, time_used, GlobalPara::boundaryCondition::_2D::inf::ruvp);
//			}
//		}
//		//按esc终止
//		if (_kbhit()) {
//			char ch = _getch();
//			if (ch == 27) {
//				char szBuffer[20];
//				sprintf_s(szBuffer, _countof(szBuffer), "%04d", istep);
//				writeContinueFile(filePathManager->getExePath_withSlash() + "output\\pause_" + GlobalPara::basic::filename + "[" + szBuffer + "].dat",t,istep);
//				std::cout << "[ESC]: Computation Interrupted\n";
//				signal_pause = true;
//			}
//		}
//		//达到规定迭代时间，终止
//		else if (t >= T || istep >= maxIteration) {
//			char szBuffer[20];
//			sprintf_s(szBuffer, _countof(szBuffer), "%04d", istep);
//			writeContinueFile(filePathManager->getExePath_withSlash() + "output\\pause_" + GlobalPara::basic::filename + "[" + szBuffer + "].dat", t, istep);
//			writeTecplotFile(filePathManager->getExePath_withSlash() + "output\\" + GlobalPara::basic::filename + "[" + szBuffer + "].dat",
//				t, rho_nodes, u_nodes, v_nodes, p_nodes);
//			std::cout << "Computation finished\n";
//			signal_pause = true;
//		}
//		//残差足够小，终止
//		else if (residual_vector[0] <= Constant::epsilon) {
//			char szBuffer[20];
//			sprintf_s(szBuffer, _countof(szBuffer), "%04d", istep);
//			writeContinueFile(filePathManager->getExePath_withSlash() + "output\\pause_" + GlobalPara::basic::filename + "[" + szBuffer + "].dat", t, istep);
//			writeTecplotFile(filePathManager->getExePath_withSlash() + "output\\" + GlobalPara::basic::filename + "[" + szBuffer + "].dat",
//				t, rho_nodes, u_nodes, v_nodes, p_nodes);
//			std::cout << "Computation finished as the field is already stable\n";
//#ifdef _WIN32
//			std::cout << "Please close the pop-up window.\n";
//			MessageBox(NULL, "Computation finished", "U2NITS", MB_OK);
//#endif // _WIN32
//			signal_pause = true;
//		}
//		//终止
//		if (signal_pause) {
//			LogWriter::writeLog(solveInfo);
//			break;
//		}
//	}
//
//}

void FVM_2D::solve_CPU2(std::string suffix_out, std::string suffix_info) {
	/*
	
	将GPU代码复制到此处后，应修改的地方包括：
	1.将GPU::GPUSolver2 gpuSolver2; 改为 Solver_2D cpuSolver2;
	2.删除初始化部分try
	3.计算部分，将gpuSolver2.iteration()改为cpuSolver2.evolve(dt)
	4.释放资源部分，删除gpuSolver2.finalize()
	具体可参照Git记录

	*/

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
	enum MySignal {
		_NoSignal,
		_EscPressed,
		_TimeReached,
		_StableReached
	};
	std::vector<Element_2D> elements_old;
	std::vector<double> u_nodes;// 节点u, size=nodes.size()。输出流场时，用作缓存。
	std::vector<double> v_nodes;
	std::vector<double> rho_nodes;
	std::vector<double> p_nodes;
	std::string solveInfo;
	std::string outputPathWithSlash = FilePathManager::getInstance()->getOutputDirectory();
	std::string basicFileName = GlobalPara::basic::filename;
	clock_t start_t;
	COORD p1;
	Solver_2D cpuSolver2;
	HistWriter histWriter(outputPathWithSlash + basicFileName + "_hist.dat");

	// 初始化
	start_t = clock();
	histWriter.writeHistFileHead();
	p1 = ConsolePrinter::getCursorPosition();
	ConsolePrinter::drawProgressBar(int(t / T));


	// 迭代
	for (int istep = istep_previous; istep <= maxIteration && t <= T; istep++) {
		// --- 更新 ---
		// 更新文件名
		char szBuffer[20];
		sprintf_s(szBuffer, _countof(szBuffer), "%04d", istep);//若istep小于4位数，则补0
		std::string str_istep_withBracket = "[" + std::string(szBuffer) + "]";
		std::string tecplotFilePath = outputPathWithSlash + basicFileName + str_istep_withBracket + ".dat";
		std::string continueFilePath = outputPathWithSlash + "pause_" + basicFileName + str_istep_withBracket + ".dat";
		std::string continueFilePath_nan = outputPathWithSlash + "nan_" + basicFileName + str_istep_withBracket + ".dat";
		// 更新旧结果
		elements_old = this->elements;
		// 更新时间和时间步长
		double dt = this->caldt(t, T);
		t += dt;

		// --- 计算 ---
		// GPU计算
		cpuSolver2.evolve(dt);

		// --- 输出 ---
		// 输出进度表示正在计算，没有卡顿
		std::cout << ".";
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
		// 检查非法值 非法则立刻跳出循环
		if (this->isNan() || residual_vector[0] > big_value) {
			// 这里直接break，因此不需要用到后面的b_writeContinue、b_writeTecplot
			// 况且writeContinueFile的参数不同，有nan
			ConsolePrinter::printInfo(ConsolePrinter::InfoType::type_nan_detected);
			//writeTecplotFile(tecplotFilePath, t);// 别删
			// 根据单元U更新节点ruvp
			this->calculateNodeValue(rho_nodes,u_nodes,v_nodes,p_nodes);
			FieldWriter::writeTecplotFile(t, tecplotFilePath, "title", nodes, elements, rho_nodes, u_nodes, v_nodes, p_nodes);
			//this->writeContinueFile(continueFilePath_nan, t, istep);// 别删
			FieldWriter::writeContinueFile(
				istep, t, continueFilePath_nan, nodes, elements, &(this->boundaryManager)
			);
			break;
		}
		// 定期输出流场
		if (istep % GlobalPara::output::step_per_output_field == 1) {
			//writeTecplotFile(tecplotFilePath, t);
			b_writeTecplot = true;
			nFiles++;
			//.autosave 自动保存文件 续算文件
			this->updateAutoSaveFile(t, istep, iCurrentAutosaveFile);

		}
		// 定期输出残差
		if (istep % GlobalPara::output::step_per_output_hist == 1) {
			//计算残差，结果保存在residual_vector中
			ResidualCalculator::cal_residual(elements_old, elements, ResidualCalculator::NORM_INF, residual_vector);
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
		// 达到规定迭代时间，终止
		else if (t >= T || istep >= maxIteration) {
			b_writeContinue = true;
			b_writeTecplot = true;
			b_pause = true;
			promptSignal = _TimeReached;
		}
		// 残差足够小，终止
		else if (residual_vector[0] <= GlobalPara::constant::epsilon) {
			b_writeContinue = true;
			b_writeTecplot = true;
			b_pause = true;
			promptSignal = _StableReached;
		}

		// 写文件操作
		if (b_writeContinue || b_writeTecplot) {
			// 根据单元U更新节点ruvp
			this->calculateNodeValue(rho_nodes, u_nodes, v_nodes, p_nodes);
			if (b_writeContinue) {
				//this->writeContinueFile(continueFilePath, t, istep);
				FieldWriter::writeContinueFile(
					istep, t, continueFilePath, nodes, elements, &(this->boundaryManager)
				);
			}
			if (b_writeTecplot) {
				//writeTecplotFile(tecplotFilePath, t);
				FieldWriter::writeTecplotFile(t, tecplotFilePath, "title", nodes, elements, rho_nodes, u_nodes, v_nodes, p_nodes);
			}
		}

		// 提示词 放在写文件之后，防止卡顿
		switch (promptSignal) {
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
		// 终止操作
		if (b_pause) {
			LogWriter::writeLog(solveInfo);
			break;
		}
	}


}

//void FVM_2D::run_GPU() {
//	//// 初始化
//	// 尝试读取续算文件。若读取失败或者_continue==false则从头开始
//	bool startFromZero = false;
//	if (GlobalPara::basic::_continue) {
//		int flag_readContinue = readContinueFile();
//		if (flag_readContinue == -1) {
//			std::string error_msg = "Info: Fail to read previous mesh(pause_*.dat). ";
//			error_msg += "Will try to start from zero again.\n";
//			LogWriter::writeLogAndCout(error_msg, LogWriter::Info);
//			startFromZero = true;
//			// 防止后面histWriter不写文件头
//			GlobalPara::basic::_continue = false;
//		}
//	}
//	else {
//		startFromZero = true;
//	}
//	// 从头开始。读取网格、初始化
//	if (startFromZero) {
//		const std::string& type = GlobalPara::basic::meshFileType;
//		if (type == "inp") {
//			std::string dir = FilePathManager::getInstance()->getExePath_withSlash() + "input\\";
//			int flag_readMesh = InpMeshReader::readGmeshFile(dir + GlobalPara::basic::filename + ".inp");
//			if (flag_readMesh == -1) {
//				std::string error_msg = "Error: Fail to read mesh. Program will exit.\n";
//				LogWriter::writeLogAndCout(error_msg);
//				return;//退出
//			}
//			setInitialCondition();
//		}
//		else if (type == "su2") {
//			std::string dir = FilePathManager::getInstance()->getExePath_withSlash() + "input\\";
//			int flag_readMesh = SU2MeshReader::readFile(dir + GlobalPara::basic::filename + ".su2", true);
//			if (flag_readMesh == -1) {
//				std::string error_msg = "Error: Fail to read mesh. Program will exit.\n";
//				LogWriter::writeLogAndCout(error_msg);
//				return;//退出
//			}
//			setInitialCondition();
//		}
//		else {
//			std::string error_msg = "Invalid mesh file type: " + type + ". Program will exit.\n";
//			LogWriter::writeLogAndCout(error_msg);
//			return;
//		}
//
//	}
//
//	// 日志记录边界参数
//	std::string str;
//	str += "BoundaryCondition:\n";
//	str += "inlet::ruvp\t" + StringProcessor::doubleArray_2_string(GlobalPara::boundaryCondition::_2D::inlet::ruvp, 4)
//		+ "\noutlet::ruvp\t" + StringProcessor::doubleArray_2_string(GlobalPara::boundaryCondition::_2D::outlet::ruvp, 4)
//		+ "\ninf::ruvp\t" + StringProcessor::doubleArray_2_string(GlobalPara::boundaryCondition::_2D::inf::ruvp, 4)
//		+ "\n";
//	LogWriter::writeLog(str);
//
//	//// 求解
//	if (GlobalPara::physicsModel::equation == _EQ_euler) {
//		solve_GPU();
//	}
//	else {
//		std::cout << "Error: Invalid equation type";
//	}
//
//}
//
//void FVM_2D::solve_GPU() {
//	/*
//	[开发中]以下为GPU代码
//	TODO:
//	*/
//
//	const int residual_vector_size = 4;
//	const int istep_previous = GlobalPara::time::istep_previous;
//	const double previous_t = GlobalPara::time::t_previous;
//	const double big_value = 1e8;
//	const int maxIteration = GlobalPara::output::maxIteration;
//	const double T = GlobalPara::time::T;//目标物理时间。用于控制是否终止计算
//	int nFiles = 0;//输出文件个数，用于显示
//	int iCurrentAutosaveFile = 0;// 自动保存文件的指标，构成循环队列
//	double residual_vector[residual_vector_size]{ 1,1,1,1 };
//	double t = previous_t;//当前物理时间。由续算时的起始时间决定
//	enum MySignal {
//		_NoSignal,
//		_EscPressed,
//		_TimeReached,
//		_StableReached
//	};
//	std::string solveInfo;
//	std::string outputPathWithSlash = FilePathManager::getInstance()->outputFolder_path + "\\";
//	std::string basicFileName = GlobalPara::basic::filename;
//	clock_t start_t;
//	COORD p1;
//	GPU::GPUSolver2 gpuSolver2;
//	HistWriter histWriter(outputPathWithSlash + basicFileName + "_hist.dat");
//	FieldWriter tecplotWriter;
//	FieldWriter continueWriter;
//
//	// 初始化
//	start_t = clock();
//	histWriter.writeHistFileHead();
//	gpuSolver2.initialze();
//	p1 = ConsolePrinter::getCursorPosition();
//	ConsolePrinter::drawProgressBar(int(t / T));
//
//
//	// 迭代
//	for (int istep = istep_previous; istep <= maxIteration && t <= T; istep++) {
//		// --- 更新 ---
//		// 更新文件名
//		char szBuffer[20];
//		sprintf_s(szBuffer, _countof(szBuffer), "%04d", istep);//若istep小于4位数，则补0
//		std::string str_istep_withBracket = "[" + std::string(szBuffer) + "]";
//		std::string tecplotFilePath = outputPathWithSlash + basicFileName + str_istep_withBracket + ".dat";
//		std::string continueFilePath = outputPathWithSlash + "pause_" + basicFileName + str_istep_withBracket + ".dat";
//		std::string continueFilePath_nan = outputPathWithSlash + "nan_" + basicFileName + str_istep_withBracket + ".dat";
//		// 仅host端的更新 用U更新ruvp、U_old
//		gpuSolver2.update_ruvp_Uold(this);
//
//		// 更新时间和时间步长
//		// 目前默认是无粘的，因此ifViscous和ifConstViscous都为false
//		double dt = gpuSolver2.calculateDt(
//			t, T, GlobalPara::constant::gamma, Constant::Re, Constant::Pr, GlobalPara::time::CFL, Constant::R,
//			false, false
//		);
//		t += dt;
//
//		//std::cout << "gpuSolver2.iterationDevice(dt)\n";
//		// --- 计算 ---
//		// GPU计算
//		gpuSolver2.iterationDevice(dt);
//		gpuSolver2.device_to_host();
//		gpuSolver2.iterationHost(dt);
//		gpuSolver2.host_to_device();
//
//		// --- 输出 ---
//		// 输出进度表示正在计算，没有卡顿
//		std::cout << ".";
//		bool b_writeContinue = false;
//		bool b_writeTecplot = false;
//		bool b_pause = false;//暂停信号，为true则终止
//		MySignal promptSignal = _NoSignal;// 提示词分类
//		// 定期输出进度
//		if (istep % GlobalPara::output::step_per_print == 1) {
//			ConsolePrinter::clearDisplay(p1, ConsolePrinter::getCursorPosition());
//			ConsolePrinter::setCursorPosition(p1);
//			ConsolePrinter::drawProgressBar(int(t / T * 100 + 2));//+x是为了防止停在99%
//
//			double calculateTime = (double)(clock() - start_t) / CLOCKS_PER_SEC;
//			double calculateSpeed = ((istep - istep_previous) / calculateTime);
//			solveInfo = ConsolePrinter::assemblySolveInfo(calculateTime, istep, maxIteration, calculateSpeed, nFiles, t, T, residual_vector);
//			std::cout << solveInfo;
//		}
//		// 发散则立刻跳出循环
//		if (residual_vector[0] > big_value) {
//			// 这里直接break，因此不需要用到后面的b_writeContinue、b_writeTecplot
//			// 况且writeContinueFile的参数不同，有nan
//			ConsolePrinter::printInfo(ConsolePrinter::InfoType::type_nan_detected);
//			//writeTecplotFile(tecplotFilePath, t);// 别删
//			// 根据单元U更新节点ruvp
//			gpuSolver2.updateOutputNodeField();
//			tecplotWriter.setFilePath(tecplotFilePath);
//			tecplotWriter.writeTecplotFile_GPU(t, "title", gpuSolver2.node_host,gpuSolver2.element_host,gpuSolver2.outputNodeField);// 完成于2024-02-28
//			// 写nan文件
//			continueWriter.writeContinueFile_GPU(
//				istep, t, continueFilePath_nan, 
//				gpuSolver2.node_host, gpuSolver2.element_host, gpuSolver2.elementField_host, &(this->boundaryManager)
//			);// 完成于2024-02-28
//			break;
//		}
//		// 定期输出流场
//		if (istep % GlobalPara::output::step_per_output_field == 1) {
//			//writeTecplotFile(tecplotFilePath, t);
//			b_writeTecplot = true;
//			nFiles++;
//			//.autosave 自动保存文件 续算文件
//			//最多保存global::autosaveFileNum个autosave文件，再多则覆写之前的
//			std::string autosaveFilePath = outputPathWithSlash + "autosave" + std::to_string(iCurrentAutosaveFile) + "_" + GlobalPara::basic::filename + ".dat";
//			continueWriter.writeContinueFile_GPU(
//				istep, t, autosaveFilePath,
//				gpuSolver2.node_host, gpuSolver2.element_host, gpuSolver2.elementField_host, &(this->boundaryManager)
//			);// 完成于2024-02-28
//			iCurrentAutosaveFile++;
//			if (iCurrentAutosaveFile >= GlobalPara::output::autosaveFileNum)iCurrentAutosaveFile = 0;
//
//		}
//		// 定期输出残差
//		if (istep % GlobalPara::output::step_per_output_hist == 1) {
//			//计算残差，结果保存在residual_vector中
//			ResidualCalculator::cal_residual_GPU(gpuSolver2.element_U_old, gpuSolver2.elementField_host, ResidualCalculator::NORM_INF, residual_vector);
//			histWriter.writeHistFileData(istep, residual_vector, residual_vector_size);
//		}
//
//		// --- 终止 ---
//		// 按esc终止
//		if (_kbhit()) {
//			if (_getch() == 27) {
//				b_writeContinue = true;
//				b_pause = true;
//				promptSignal = _EscPressed;
//			}
//		}
//		// 达到规定迭代时间或迭代步数，终止
//		else if (t >= T || istep >= maxIteration) {
//			b_writeContinue = true;
//			b_writeTecplot = true;
//			b_pause = true;
//			promptSignal = _TimeReached;
//		}
//		// 残差足够小，终止
//		else if (residual_vector[0] <= Constant::epsilon) {
//			b_writeContinue = true;
//			b_writeTecplot = true;
//			b_pause = true;
//			promptSignal = _StableReached;
//		}
//
//		// 写文件操作
//		if (b_writeContinue || b_writeTecplot) {
//			// 根据单元U更新节点ruvp
//			gpuSolver2.updateOutputNodeField();
//			if (b_writeContinue) {
//				continueWriter.writeContinueFile_GPU(
//					istep, t, continueFilePath,
//					gpuSolver2.node_host, gpuSolver2.element_host, gpuSolver2.elementField_host, &(this->boundaryManager)
//				);
//			}
//			if (b_writeTecplot) {
//				tecplotWriter.setFilePath(tecplotFilePath);
//				tecplotWriter.writeTecplotFile_GPU(t, "title", gpuSolver2.node_host, gpuSolver2.element_host, gpuSolver2.outputNodeField);
//			}
//		}
//
//		// 提示词 放在写文件之后，防止卡顿
//		switch (promptSignal) {
//		case _NoSignal:
//			break;
//		case _EscPressed:
//			std::cout << "[ESC]: Computation Interrupted\n";
//			break;
//		case _TimeReached:
//			std::cout << "Computation finished\n";
//			break;
//		case _StableReached:
//			std::cout << "Computation finished as the field is already stable\n";
//			break;
//		default:
//			break;
//		}
//		// 终止操作
//		if (b_pause) {
//			LogWriter::writeLog(solveInfo);
//			break;
//		}
//	}
//
//	// 释放资源
//	gpuSolver2.finalize();
//}

double FVM_2D::caldt(double t, double T) {
	// 根据单元LambdaC、单元面积、CFL数计算dt_i，取最小值，返回
	// 输入：elements、CFL、当前时间t、目标时间T
	// 输出：dt
	double dt = 1e10;// a big value
	// 以下为新代码
	for (int ie = 0; ie < elements.size(); ie++) {
		double LambdaC_i = elements[ie].calLambdaFlux(this);
		double Omega_i = elements[ie].calArea(this);
		double dt_i = GlobalPara::time::CFL * Omega_i / LambdaC_i;
		dt = (std::min)(dt, dt_i);
	}
	// 如果t+dt超出设定T，则限制dt
	dt = (t + dt > T) ? (T - t) : dt;//if (t + dt > T)dt = T - t;
	return dt;
}

//double FVM_2D::caldt_GPU(double t, double T, void* gpuSolver2) {
//	REAL Re = 1e8;
//	REAL Pr = 0.73;
//	// GlobalPara::physicsModel::equation;
//	// 目前默认是无粘的，因此ifViscous和ifConstViscous都为false
//	GPU::GPUSolver2* g = (GPU::GPUSolver2*)gpuSolver2;
//	return GPU::calculateDt(
//		t, T, GlobalPara::constant::gamma, Re, Pr, GlobalPara::time::CFL, Constant::R,
//		false, false,
//		(*g)
//	);
//}

void FVM_2D::writeTecplotFile(std::string f_name, double t_current, 
	std::vector<double>& rho_nodes,
	std::vector<double>& u_nodes,
	std::vector<double>& v_nodes,
	std::vector<double>& p_nodes
	) {
	// 根据单元U更新节点ruvp
	//calculateNodeValue();

	std::ofstream f(f_name);
	if (!f.is_open())std::cout << "Error: Fail to open file " << f_name << std::endl;

	//header
	f << R"(TITLE = "Example: 2D Finite Element Data")" << "\n";
	f << R"(VARIABLES = "X", "Y")";
	if (GlobalPara::output::output_var_ruvp[0])		f << R"(, "rho")";
	if (GlobalPara::output::output_var_ruvp[1])		f << R"(, "u")";
	if (GlobalPara::output::output_var_ruvp[2])		f << R"(, "v")";
	if (GlobalPara::output::output_var_ruvp[3])		f << R"(, "p")";
	f << "\n";

	f << "ZONE NODES="<< nodes.size() 
		<<", ELEMENTS = "<< elements.size() 
		<<", DATAPACKING = POINT, ZONETYPE = FEQUADRILATERAL"
		<<", SOLUTIONTIME="<<std::scientific<< t_current<<std::fixed
		<<std::endl;
	
	//nodes
	for (int i = 0; i < nodes.size(); i++) {
		f << nodes[i].x << ' ' << nodes[i].y;
		if (GlobalPara::output::output_var_ruvp[0])		f << ' ' << rho_nodes[i];
		if (GlobalPara::output::output_var_ruvp[1])		f << ' ' << u_nodes[i];
		if (GlobalPara::output::output_var_ruvp[2])		f << ' ' << v_nodes[i];
		if (GlobalPara::output::output_var_ruvp[3])		f << ' ' << p_nodes[i];
		f << std::endl;
	}

	//elements
	for (int ie = 0; ie < elements.size(); ie++) {
		f << elements[ie].nodes[0] << " ";
		f << elements[ie].nodes[1] << " ";
		f << elements[ie].nodes[2] << " ";
		f << elements[ie].nodes[2] << std::endl;
	}

	f.close();
}

void FVM_2D::writeContinueFile(std::string f_name, double t, int istep) {
	//std::cout << "wirtefile:" << f_name << std::endl;
	int flag_OUT = -1;
	std::ofstream outfile(f_name);
	if (!outfile.is_open())std::cout << "Error: Fail to open file " << f_name << std::endl;

	outfile << "t, istep" << std::endl;
	outfile << t << " " << istep << std::endl;

	outfile << "nodes: ID, x, y" << std::endl;
	for (int in = 0; in < nodes.size(); in++) {
		//std::vector<double> uu_i = nodes[in].calNodeValue();
		outfile << nodes[in].ID << " ";
		outfile << nodes[in].x << " ";
		outfile << nodes[in].y << std::endl;
	}

	outfile << "elements: ID, nodes, U" << std::endl;
	for (int ie = 0; ie < elements.size(); ie++) {
		outfile << elements[ie].ID << " ";
		outfile << elements[ie].nodes[0] << " ";
		outfile << elements[ie].nodes[1] << " ";
		outfile << elements[ie].nodes[2] << " ";
		outfile << elements[ie].U[0] << " ";
		outfile << elements[ie].U[1] << " ";
		outfile << elements[ie].U[2] << " ";
		outfile << elements[ie].U[3] << std::endl;
	}

	outfile << "boundaries: ID, name; edgeIDs" << std::endl;
	for (int ib = 0; ib < boundaryManager.boundaries.size(); ib++) {
		outfile << boundaryManager.boundaries[ib].ID << " ";
		outfile << boundaryManager.boundaries[ib].name << std::endl;
		//一行太长会导致读取时buffer长度不够，从而导致难以预料的结果(死循环、读取失败等)
		int nEdge = (int)boundaryManager.boundaries[ib].pEdges.size();//控制输出' '还是'\n'
		bool check = 0;//控制
		for (int iEdge = 0; iEdge < nEdge;iEdge++) {
			outfile << boundaryManager.boundaries[ib].pEdges[iEdge]->ID;
			if (iEdge % 10 == 9){
				outfile << std::endl; 
				check = 1;
			}
			else {
				outfile << ' ';
				check = 0;
			}
		}
		if(!check)outfile << std::endl;//若上次输出的' '
	}
	outfile.close();

}

int FVM_2D::readContinueFile() {
	std::string m_path = FilePathManager::getInstance()->getOutputDirectory();
	std::vector<std::string> files = FilePathManager::getInstance()->getFiles(m_path);//files为字符串向量，存储了路径下所有文件名
	int index_maxNstep = -1;
	//搜寻是否有pause_filename[xxx].dat文件，若有，则找xxx最大的
	int maxNstep = 0;
	int tmp_index = 0;//用于指示"[", "]"
	std::string tmp_str;
	for (int i = 0; i < files.size(); i++) {
		//搜寻以pause_filename开头的文件名
		std::string str_match = "pause_" + GlobalPara::basic::filename;
		if (files[i].substr(0, int(str_match.size())) == str_match) {
			//剔除"pause_filename["
			int pre_length = int(str_match.size() + 1);//"pause_filename["的长度
			std::string str = files[i].substr(pre_length);
			//剔除"].dat"
			int post_length = 5;//"].dat"的长度
			for (int i = 0; i < post_length; i++) {
				str.pop_back();
			}
			//将xxx存入num
			int num = std::stoi(str);
			if (num > maxNstep) {
				index_maxNstep = i;
				maxNstep = num;
			}
		}
	}
	//若无，则返回-1
	if (index_maxNstep == -1)return -1;
	std::ifstream infile(m_path + files[index_maxNstep]);
	if (!infile) {
		return -1;
	}

	std::cout << "Continue from: "<< files[index_maxNstep] << std::endl;

	int state = -1;//1-Node, 2-Element
	int maxnodeID = 1;
	int maxedgeID = 1;
	int maxelementID = 1;
	int iline_set = 0;//临时变量，用于记录读取set时经过了多少行
	std::vector<std::vector<int>> edges_of_all_sets;//临时变量，存储各set的edgeID
	const int bufferLength = 600;//!!!小心长度不够
	char buffer[bufferLength];
	std::string tLine;
	std::vector<std::string> tWords;
	std::vector<std::vector<std::string>> tWordsMatrix;
	std::vector<int> edge_ID_long;
	int nNullRow = 0;
	bool loop = 1;
	while (loop) {
		//get words and set state
		infile.getline(buffer, bufferLength);
		if (infile.eof())loop=0;
		tLine = buffer;
		tWords = StringProcessor::splitString(tLine);
		if (nNullRow >= 10)loop = 0;
		if (tWords.size() == 0) {
			nNullRow++;
			continue;
		}//空行，强迫开始下一次循环，防止tWords[0]内存错误
		else {
			if (tWords[0] == "t")state = 0;
			if (tWords[0] == "nodes:")state = 1;
			if (tWords[0] == "edges:")state = 2;
			if (tWords[0] == "elements:")state = 3;
			if (tWords[0] == "boundaries:")state = 4;
		}

		//translate
		if (state == 0 && tWords[0].substr(0, 1) != "t") {//t_solve, istep_solve
			GlobalPara::time::t_previous = std::stod(tWords[0]);
			GlobalPara::time::istep_previous = (int)std::stod(tWords[1]);
		}
		else if (state == 1 && tWords[0].substr(0, 1) != "n") {//nodes: ID, x, y
			Node_2D node;
			node.ID = (int)std::stod(tWords[0]);
			node.x = std::stod(tWords[1]);
			node.y = std::stod(tWords[2]);
			nodes.push_back(node);
			maxnodeID = (std::max)(maxnodeID, node.ID);
		}
		else if (state == 3 && tWords[0].substr(0, 1) != "e") {//elements: ID, nodes, U
			Element_2D ele;
			ele.ID = (int)std::stod(tWords[0]);
			ele.nodes[0] = (int)std::stod(tWords[1]);
			ele.nodes[1] = (int)std::stod(tWords[2]);
			ele.nodes[2] = (int)std::stod(tWords[3]);
			ele.U[0] = std::stod(tWords[4]);
			ele.U[1] = std::stod(tWords[5]);
			ele.U[2] = std::stod(tWords[6]);
			ele.U[3] = std::stod(tWords[7]);
			elements.push_back(ele);
			maxelementID = (std::max)(maxelementID, ele.ID);
		}
		else if (state == 4 && tWords[0].substr(0, 1) != "b") {//boundaries: ID, name, edgeIDs
			iline_set++;
			if (!isdigit(tWords[1][0])) {
				/*
				"1 inf" 
				1.vBoundarySets新增条目、该条目的部分初始化
				2.edges_of_all_sets新增条目
				*/
				//set的ID,name初始化
				VirtualBoundarySet_2D vb;
				vb.ID = (int)std::stod(tWords[0]);
				vb.name = tWords[1];
				boundaryManager.boundaries.push_back(vb);
				//edges_of_all_sets新增条目
				edges_of_all_sets.push_back(std::vector<int>());
			}
			else {
				// 1.edges_of_all_sets的初始化
				// 引用最后一个BoundarySet
				std::vector<int>& edges_of_current_set = edges_of_all_sets[edges_of_all_sets.size() - 1];
				std::vector<int> intVector = StringProcessor::stringVector_2_intVector(tWords);
				VectorProcessor::appendToVector(edges_of_current_set, intVector);
			}
		}
	}

	infile.close();

	iniPNodeTable(maxnodeID);
	iniEdges();
	iniPEdgeTable();
	iniPElementTable(maxelementID);
	iniElement_xy_pEdges();
	iniNode_neighborElements();
	iniEdges_lengths();


	// 初始化boundaryManager.boundaries的pEdges，但不初始化pEdge所指的edge的setID等信息，该部分工作留给后面函数
	// 前置条件：需要初始化pEdgeTable
	if (this->pEdgeTable.size() == 0) {
		LogWriter::writeLogAndCout("Error: uninitialized pEdgeTable. (BoundaryManager_2D::iniBoundarySetPEdges)\n");
		throw "uninitialized pEdgeTable";
	}
	else {
		//初始化boundaryManager.boundaries的每个set的pEdges
		for (int iset = 0; iset < edges_of_all_sets.size(); iset++) {
			std::vector<int>& edgeIDs = edges_of_all_sets[iset];//第iset条set的edgeIDs
			for (int iw = 0; iw < edgeIDs.size(); iw++) {//第iset条set的第iw个edge
				//在f->pEdgeTable中，根据edgeID查询pEdge，完成初始化
				int edge_ID = edgeIDs[iw];
				Edge_2D* pE = this->pEdgeTable[edge_ID];
				this->boundaryManager.boundaries[iset].pEdges.push_back(pE);
			}
		}
	}

	// 设置edges中的setID，设置boundaryManager.boundaries的type
	// 前置条件：有edges向量，有boundaries向量，且boundaries有name、pEdges、ID
	{
		int ret = getInstance()->boundaryManager.iniBoundaryEdgeSetID_and_iniBoundaryType(getInstance());
		if (ret != 0)return ret;
	}

	return 0;
}

void FVM_2D::updateAutoSaveFile(double t, int istep, int& iCurrentAutosaveFile) {
	//最多保存global::autosaveFileNum个autosave文件，再多则覆写之前的
	writeContinueFile(FilePathManager::getInstance()->getOutputDirectory() + "autosave" + std::to_string(iCurrentAutosaveFile) + "_" + GlobalPara::basic::filename + ".dat", t, istep);
	iCurrentAutosaveFile++;
	if (iCurrentAutosaveFile >= GlobalPara::output::autosaveFileNum)iCurrentAutosaveFile = 0;
}

void FVM_2D::calculateNodeValue(
	std::vector<double>& rho_nodes,
	std::vector<double>& u_nodes,
	std::vector<double>& v_nodes,
	std::vector<double>& p_nodes) {

	if (GlobalPara::output::output_var_ruvp[0]) 		rho_nodes.resize(nodes.size());
	if (GlobalPara::output::output_var_ruvp[1]) 		u_nodes.resize(nodes.size());
	if (GlobalPara::output::output_var_ruvp[2]) 		v_nodes.resize(nodes.size());
	if (GlobalPara::output::output_var_ruvp[3]) 		p_nodes.resize(nodes.size());

	for (int i = 0; i < nodes.size(); i++) {
		// 计算节点函数值。取邻居单元的分布函数的平均值
		double U[4];
		nodes[i].calNodeValue(U);
		double ruvp[4];
		Math_2D::U_2_ruvp(U, ruvp, GlobalPara::constant::gamma);
		if (GlobalPara::output::output_var_ruvp[0]) 		rho_nodes[i] = ruvp[0];
		if (GlobalPara::output::output_var_ruvp[1]) 		u_nodes[i] = ruvp[1];
		if (GlobalPara::output::output_var_ruvp[2]) 		v_nodes[i] = ruvp[2];
		if (GlobalPara::output::output_var_ruvp[3]) 		p_nodes[i] = ruvp[3];
	}

}
/*
void FVM_2D::calculateNodeValue_GPU(void* pSolver) {

	GPU::GPUSolver2* gpuSolver = (GPU::GPUSolver2*) pSolver;
	auto& node_ruvp = gpuSolver->outputNodeField.ruvp;
	int num_node = gpuSolver->outputNodeField.num_node;
	int num_element = gpuSolver->element_host.num_element;
	const int nNodePerElement = 4;
	const int nValuePerNode = 4;
	int* node_neighborElement_num;// 记录每个节点的邻居单元数量
	// 申请资源 初始化
	node_neighborElement_num = new int[num_node]();// 小括号自动初始化为0
	for (int i = 0; i < 4; i++) {
		memset(node_ruvp[i], 0, num_node * sizeof(REAL));
	}

	// 将单元值加到其子节点值上
	for (int iElement = 0; iElement < num_element; iElement++) {
		// 统计每个节点的邻居单元数量，并修改节点值
		for (int jNode = 0; jNode < nNodePerElement; jNode++) {
			// 获取单元节点ID
			int GPUID_of_node = gpuSolver->element_host.nodes[jNode][iElement];// GPUID of node 0
			// nodeID[3]=-1时，跳过该节点
			if (GPUID_of_node < 0 || GPUID_of_node >= num_node)continue;// 跳过本次循环，但不跳出循环体
			// 该ID对应的邻居单元数量+1
			node_neighborElement_num[GPUID_of_node]++;
			// 该ID对应的节点所有值加上单元值
			for (int kValue = 0; kValue < nValuePerNode; kValue++) {
				node_ruvp[kValue][GPUID_of_node] += gpuSolver->element_vruvp[kValue][iElement];
			}
		}
	}

	// 节点值除以邻居单元数，得到平均值，作为节点ruvp值
	for (int iNode = 0; iNode < num_node; iNode++) {
		// 为了避免除以0，分母等于0则跳过
		if (node_neighborElement_num[iNode] == 0)continue;
		// node_ruvp除以邻居单元数，得到平均值
		for (int kValue = 0; kValue < nValuePerNode; kValue++) {
			node_ruvp[kValue][iNode] /= node_neighborElement_num[iNode];
		}
	}

	// 释放资源
	//for (int i = 0; i < 4; i++) {
	//	delete[] node_ruvp[i];
	//}
	delete[] node_neighborElement_num;
}
*/
void FVM_2D::iniPNodeTable(int maxnodeID) {
	pNodeTable.resize(maxnodeID + 1);
	for (int in = 0; in < nodes.size(); in++) {
		pNodeTable[nodes[in].ID] = &(nodes[in]);
	}
}

void FVM_2D::iniEdges() {
	for (int ie = 0; ie < elements.size(); ie++) {
		int n0, n1, n2;
		n0 = elements[ie].nodes[0];
		n1 = elements[ie].nodes[1];
		n2 = elements[ie].nodes[2];
		iniEdges_registerSingle(n0, n1, &(elements[ie]));
		iniEdges_registerSingle(n1, n2, &(elements[ie]));
		iniEdges_registerSingle(n2, n0, &(elements[ie]));
	}
}

void FVM_2D::iniPEdgeTable() {
	pEdgeTable.resize(edges.size() + 1);
	for (int i = 0; i < edges.size(); i++) {
		pEdgeTable[edges[i].ID] = &(edges[i]);
	}
}

void FVM_2D::iniNode_neighborElements() {
	// 前置条件：elements, elements.nodes, pNodeTable
	for (int ie = 0; ie < elements.size(); ie++) {
		for (int jn = 0; jn < 3; jn++) {
			Node_2D* pN = getNodeByID(elements[ie].nodes[jn]);
			pN->neighborElements.push_back(&(elements[ie]));
		}
	}
}

Edge_2D* FVM_2D::getEdgeByNodeIDs(int n0, int n1) {
	// 在edges中根据nodeIDs查询edge
	//n0, n1表示节点ID
	if(edges.size()==0)return nullptr;
	for (int i = 0; i < edges.size(); i++) {
		if (edges[i].nodes[0] == n0 && edges[i].nodes[1] == n1)return &(edges[i]);
		if (edges[i].nodes[0] == n1 && edges[i].nodes[1] == n0)return &(edges[i]);
	}
	return nullptr;
}

void FVM_2D::iniEdges_registerSingle(int n0, int n1, Element_2D* pE) {
	Edge_2D* pEdge = getEdgeByNodeIDs(n0, n1);
	if (pEdge == nullptr) {
		Edge_2D tmp_edge;
		tmp_edge.ID = (int)edges.size() + 1;
		tmp_edge.nodes[0] = n0;
		tmp_edge.nodes[1] = n1;
		tmp_edge.pElement_L = pE;
		edges.push_back(tmp_edge);
	}
	else {
		/*
		隐含假设：一个edge如果已经属于某个element A，则A是该edge的左element
		一个edge最多同时属于2个element
		*/ 
		pEdge->pElement_R = pE;
	}
}

void FVM_2D::iniPElementTable(int maxelementID) {
	pElementTable.resize(maxelementID + 1);
	for (int ie = 0; ie < elements.size(); ie++) {
		pElementTable[elements[ie].ID] = &(elements[ie]);
	}
}

void FVM_2D::iniElement_xy_pEdges() {
	/*
	初始化单元的x,y,pEdges，便于以后查找
	前置条件：
	  需要有elements数组
	  需要已知element的nodeIDs
	  node的坐标已初始化
	  有edges数组，且edge的nodeIDs已初始化
	隐患：
	  三角形边界，pEdges只有3个。对于四边形，pEdges的node需要重新确定
	*/

	for (int ie = 0; ie < elements.size(); ie++) {
		//已经calxy，以后不必担心。但是必须要先读取node再读取element

		Element_2D* pElement = &(elements[ie]);
		//初始化单元中心坐标
		double tmp_x = 0;
		double tmp_y = 0;
		for (int i = 0; i < 3; i++) {
			tmp_x += this->getNodeByID(pElement->nodes[i])->x;
			tmp_y += this->getNodeByID(pElement->nodes[i])->y;
		}
		pElement->x = tmp_x / 3.0;
		pElement->y = tmp_y / 3.0;


		pElement->pEdges[0] = this->getEdgeByNodeIDs(pElement->nodes[0], pElement->nodes[1]);
		pElement->pEdges[1] = this->getEdgeByNodeIDs(pElement->nodes[1], pElement->nodes[2]);
		pElement->pEdges[2] = this->getEdgeByNodeIDs(pElement->nodes[2], pElement->nodes[0]);


	}
	hasInitElementXY = true;//已经初始化单元中心坐标。
}

void FVM_2D::iniEdges_lengths() {
	if (hasInitElementXY == false) {
		LogWriter::writeLogAndCout("Error: hasInitElementXY == false, in FVM_2D::iniEdges_lengths()\n");
	}
	for (int i = 0; i < edges.size(); i++) {
		edges[i].length = (float)edges[i].getLength();
		edges[i].refLength = (float)edges[i].getRefLength();
	}
	hasInitEdgeLengths = true;
}

Node_2D* FVM_2D::getNodeByID(int ID) {
	// 防止越界
	if (ID < 0 || ID >= pNodeTable.size()) {
		std::string error_msg = "Error: ID out of range in FVM_2D::getNodeByID(), ID=" + std::to_string(ID) + "\n";
		LogWriter::writeLogAndCout(error_msg);
		throw error_msg;
	}
	return pNodeTable[ID];
}

void FVM_2D::isentropicVortex(double x, double y, double xc, double yc, double chi, double& deltau, double& deltav, double& deltaT) {
	double xbar = x - xc;
	double ybar = y - yc;
	double r2 = xbar * xbar + ybar * ybar;
	const double gamma = GlobalPara::constant::gamma;
	const double PI = GlobalPara::constant::PI;
	deltau = chi / 2.0 / PI * exp(0.5 * (1 - r2)) * (-ybar);
	deltav = chi / 2.0 / PI * exp(0.5 * (1 - r2)) * xbar;
	deltaT = -(gamma - 1) / chi / chi * 8 * gamma * PI * PI * exp(1 - r2);
}

void FVM_2D::isentropicVortex_2(double xc, double yc, double chi, const double* ruvp0) {
	//ruvp0：均匀流参数
	for (int ie = 0; ie < elements.size(); ie++) {
		Element_2D& e = elements[ie];
		double rho, u, v, p;
		double xbar, ybar, r2, du, dv, dT;
		const double PI = GlobalPara::constant::PI;
		const double gamma = GlobalPara::constant::gamma;
		xbar = e.x - xc;
		ybar = e.y - yc;
		r2 = xbar * xbar + ybar * ybar;
		du = chi / 2. / PI * exp(0.5 * (1. - r2)) * (-ybar);
		dv = chi / 2. / PI * exp(0.5 * (1. - r2)) * xbar;
		u = ruvp0[1] + du;
		v = ruvp0[2] + dv;
		dT = -(gamma - 1.) * chi * chi / (8. * gamma * PI * PI) * exp(1. - r2);
		rho = pow(ruvp0[3] + dT, 1. / (gamma - 1.));
		p = rho * (ruvp0[3] + dT);
		double ruvp[4]{ rho,u,v,p };
		Math_2D::ruvp_2_U(ruvp, e.U, gamma);
	}
	//lambda表达式 [ capture ] ( params ) opt -> ret { body; }; http://c.biancheng.net/view/3741.html
	//例如 auto f = [](int a) -> int { return a + 1; };
	//auto fun_WA = [](const double* xy) {
}


void FVM_2D::writeFileHeader_isentropicVortex() {
	std::ofstream outfile(FilePathManager::getInstance()->getOutputDirectory() + "error_isentropicVortex_" + GlobalPara::basic::filename + ".txt", std::ios::app);//追加模式
	outfile << SystemInfo::getCurrentDateTime() << "\n";
	outfile
		<< "istep" << "\t"
		<< "cpu_time[s]" << "\t"
		<< "error_L1" << "\t"
		<< "error_L2" << "\t"
		<< "error_max"<< "\n";
	outfile.close();
}

bool FVM_2D::isStable(std::vector<Element_2D> old) {
	double sum = 0;
	double res;
	for (int ie = 0; ie < elements.size(); ie++) {
		//sum += abs(elements[ie].calculateu() - old[ie].calculateu());
		res = elements[ie].U[1] - old[ie].U[1];
		sum += abs(res);
	}
	//std::cout << sum << std::endl;
	if (sum <= 1e-12)return 1;
	else return 0;
}

bool FVM_2D::isNan() {
	// 检查每个元素的守恒量，输出所有异常值
	bool is_nan = 0;
	for (int ie = 0; ie < elements.size(); ie++) {
		
		std::string str;
		if (isnan(elements[ie].U[0])) {// 这里调用的是系统的isnan函数
			is_nan = 1;
			const double& x = elements[ie].x;
			const double& y = elements[ie].y;
			str = 
				"Warning: rho ==\"NaN\", in element (x=" + std::to_string(x) + ", y=" + std::to_string(y) 
				+ ", U[0,1,2,3]=" + std::to_string(elements[ie].U[0]) + ", " + std::to_string(elements[ie].U[1])
				+ ", " + std::to_string(elements[ie].U[2]) + ", " + std::to_string(elements[ie].U[3]) + "\n";
			//break;
		}
		else if(elements[ie].U[0] < 0) {
			const double& x = elements[ie].x;
			const double& y = elements[ie].y;
			str =
				"Warning: rho < 0, in element (x=" + std::to_string(x) + ", y=" + std::to_string(y)
				+ ", U[0,1,2,3]=" + std::to_string(elements[ie].U[0]) + ", " + std::to_string(elements[ie].U[1])
				+ ", " + std::to_string(elements[ie].U[2]) + ", " + std::to_string(elements[ie].U[3]) + "\n";
		}
		LogWriter::writeLogAndCout(str);
	}
	return is_nan;
}

//bool FVM_2D::isNan_GPU(void* gpuSolver2) {
//	/*
//	GPU运算时不会给出NaN值，因此这里不需要检查
//	*/
//	//bool is_nan = 0;
//	//GPU::GPUSolver2* pSolver = (GPU::GPUSolver2*)gpuSolver2;
//	//for (int ie = 0; ie < pSolver->element_host.num_element; ie++) {
//	//	std::string str;
//	//	
//	//}
//	GPU::GPUSolver2* pSolver = (GPU::GPUSolver2*)gpuSolver2;
//	double value = pSolver->elementField_host.U[0][0];
//	if (value < 0) {
//		std::cout << "Warning: rho < 0, in GPU solver, value=" << value << std::endl;
//	}
//
//	//return is_nan;
//
//	return false;
//}
