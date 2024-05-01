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

//void FVM_2D::setInitialCondition() {
//	switch (GlobalPara::initialCondition::type) {
//	case 1:
//	{
//		//1.常数
//		using namespace GlobalPara::boundaryCondition::_2D;
//		for (int ie = 0; ie < elements.size(); ie++) {
//			Math_2D::ruvp_2_U(inf::ruvp, elements[ie].U, GlobalPara::constant::gamma);
//			//if(elements[ie].ID==173)system()
//		}
//		std::string str;
//		str += "InitialCondition = Uniform flow, ruvp\t" + StringProcessor::doubleArray_2_string(inf::ruvp, 4) + "\n";
//		LogWriter::log(str);
//	}
//		break;
//
//	case 2:
//	{
//		//2.均匀流和等熵涡的叠加
//
//		using namespace GlobalPara::boundaryCondition::_2D;
//		std::string str;
//		str += "InitialCondition = Uniform flow + isentropicVortex, ruvp of Uniform flow:" + StringProcessor::doubleArray_2_string(inf::ruvp, 4) + "\n";
//		LogWriter::logAndPrint(str, LogWriter::Info);
//		isentropicVortex_2(5, 5, 5, inf::ruvp);
//	}
//		break;
//
//	case 3:
//	{
//		using namespace GlobalPara::boundaryCondition::_2D;
//		for (int ie = 0; ie < elements.size(); ie++) {
//			if (elements[ie].x < 0) {
//				Math_2D::ruvp_2_U(inlet::ruvp, elements[ie].U, GlobalPara::constant::gamma);
//			}
//			else {
//				Math_2D::ruvp_2_U(outlet::ruvp, elements[ie].U, GlobalPara::constant::gamma);
//			}
//		}
//		std::stringstream ss;
//		ss << "InitialCondition = Shock Tube.\n";
//		ss << "ruvp_inlet = " << StringProcessor::doubleArray_2_string(inlet::ruvp, 4) << ",";
//		ss << "ruvp_outlet = " << StringProcessor::doubleArray_2_string(outlet::ruvp, 4) << "\n";
//		LogWriter::logAndPrint(ss.str(), LogWriter::Info, LogWriter::Info);
//
//	}
//
//	case 1001:
//	{
//		// sod激波管(会自动设置boundaryCondition)
//		// 参照论文shock_tube，见D:\tgl\Local\CFD\shock_tube_code
//		// 事实上计算远场边界时也是用的ruvp，如果一开始use_ruvp=false
//		// 则会用Ma等计算ruvp，参加TomlFileManager::initialize_ruvp()
//		// 综上，use_ruvp没必要修改，因为后面也用不上
//		std::stringstream ss;
//		ss << "InitialCondition = SOD Shock Tube. Boundary condition is automatically set.\n";
//		LogWriter::logAndPrint(ss.str(), LogWriter::Info, LogWriter::Info);
//		real ruvpL[4]{ 1,0,0,1 };
//		real ruvpR[4]{ 0.125,0,0,0.1 };
//		using namespace GlobalPara::boundaryCondition::_2D;
//		// 修改边界条件
//		for (int i = 0; i < 4; i++) {
//			inlet::ruvp[i] = ruvpL[i];
//			outlet::ruvp[i] = ruvpR[i];
//		}
//		// 初场
//
//	}
//
//	default:
//	{
//		std::stringstream ss;
//		ss << "Error: Invalid initialCondition type " << GlobalPara::initialCondition::type << ".\n";
//		LogWriter::logAndPrint(ss.str(), LogWriter::Error, LogWriter::Error);
//		exit(GlobalPara::initialCondition::type);
//	}
//	}
//
//}

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
//			LogWriter::log(solveInfo);
//			break;
//		}
//	}
//
//}

//void FVM_2D::solve_CPU2(std::string suffix_out, std::string suffix_info) {
//	/*
//	
//	将GPU代码复制到此处后，应修改的地方包括：
//	1.将GPU::GPUSolver2 gpuSolver2; 改为 Solver_2D cpuSolver2;
//	2.删除初始化部分try
//	3.计算部分，将gpuSolver2.iteration()改为cpuSolver2.evolve(dt)
//	4.释放资源部分，删除gpuSolver2.finalize()
//	具体可参照Git记录
//
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
//	std::vector<Element_2D> elements_old;
//	std::vector<double> u_nodes;// 节点u, size=nodes.size()。输出流场时，用作缓存。
//	std::vector<double> v_nodes;
//	std::vector<double> rho_nodes;
//	std::vector<double> p_nodes;
//	std::string solveInfo;
//	std::string outputPathWithSlash = FilePathManager::getInstance()->getOutputDirectory();
//	std::string basicFileName = GlobalPara::basic::filename;
//	clock_t start_t;
//	COORD p1;
//	Solver_2D cpuSolver2;
//	HistWriter histWriter(outputPathWithSlash + basicFileName + "_hist.dat");
//
//	// 初始化
//	start_t = clock();
//	histWriter.writeHistFileHead();
//	p1 = ConsolePrinter::getCursorPosition();
//	//ConsolePrinter::drawProgressBar(int(t / T));
//
//
//	// 迭代
//	const int offset = 0;
//	for (int istep = istep_previous; istep <= maxIteration && t <= T; istep++) {
//		// --- 更新 ---
//		// 更新文件名
//		char szBuffer[20];
//		sprintf_s(szBuffer, _countof(szBuffer), "%04d", istep);//若istep小于4位数，则补0
//		std::string str_istep_withBracket = "[" + std::string(szBuffer) + "]";
//		std::string tecplotFilePath = outputPathWithSlash + basicFileName + str_istep_withBracket + ".dat";
//		std::string continueFilePath = outputPathWithSlash + "pause_" + basicFileName + str_istep_withBracket + ".dat";
//		std::string continueFilePath_nan = outputPathWithSlash + "nan_" + basicFileName + str_istep_withBracket + ".dat";
//		// 更新旧结果
//		elements_old = this->elements;
//		// 更新时间和时间步长
//		double dt = this->caldt(t, T);
//		t += dt;
//
//		// --- 计算 ---
//		// GPU计算
//		cpuSolver2.evolve(dt);
//
//		// --- 输出 ---
//		// 输出进度表示正在计算，没有卡顿
//		std::cout << ".";
//		bool b_writeContinue = false;
//		bool b_writeTecplot = false;
//		bool b_pause = false;//暂停信号，为true则终止
//		MySignal promptSignal = _NoSignal;// 提示词分类
//		// 定期输出进度
//		if (istep % GlobalPara::output::step_per_print == offset) {
//			ConsolePrinter::clearDisplay(p1, ConsolePrinter::getCursorPosition());
//			ConsolePrinter::setCursorPosition(p1);
//			ConsolePrinter::drawProgressBar(int(double(istep) / double(maxIteration) * 100.0 + 2));;//+x是为了防止停在99%
//
//			double calculateTime = (double)(clock() - start_t) / CLOCKS_PER_SEC;
//			double calculateSpeed = ((istep - istep_previous) / calculateTime);
//			solveInfo = ConsolePrinter::assemblySolveInfo(calculateTime, istep, maxIteration, calculateSpeed, nFiles, t, T, residual_vector);
//			std::cout << solveInfo;
//		}
//		// 检查非法值 非法则立刻跳出循环
//		if (this->isNan() || residual_vector[0] > big_value) {
//			// 这里直接break，因此不需要用到后面的b_writeContinue、b_writeTecplot
//			// 况且writeContinueFile的参数不同，有nan
//			ConsolePrinter::printInfo(ConsolePrinter::InfoType::type_nan_detected);
//			//writeTecplotFile(tecplotFilePath, t);// 别删
//			// 根据单元U更新节点ruvp
//			this->calculateNodeValue(rho_nodes,u_nodes,v_nodes,p_nodes);
//			FieldWriter::writeTecplotFile(t, tecplotFilePath, "title", nodes, elements, rho_nodes, u_nodes, v_nodes, p_nodes);
//			//this->writeContinueFile(continueFilePath_nan, t, istep);// 别删
//			FieldWriter::writeContinueFile(
//				istep, t, continueFilePath_nan, nodes, elements, &(this->boundaryManager)
//			);
//			break;
//		}
//		// 定期输出流场
//		if (istep % GlobalPara::output::step_per_output_field == offset) {
//			//writeTecplotFile(tecplotFilePath, t);
//			b_writeTecplot = true;
//			nFiles++;
//			//.autosave 自动保存文件 续算文件
//			this->updateAutoSaveFile(t, istep, iCurrentAutosaveFile);
//
//		}
//		// 定期输出残差
//		if (istep % GlobalPara::output::step_per_output_hist == offset) {
//			//计算残差，结果保存在residual_vector中
//			ResidualCalculator::cal_residual(elements_old, elements, ResidualCalculator::NORM_INF, residual_vector);
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
//		// 达到规定迭代时间，终止
//		else if (t >= T || istep >= maxIteration) {
//			b_writeContinue = true;
//			b_writeTecplot = true;
//			b_pause = true;
//			promptSignal = _TimeReached;
//		}
//		// 残差足够小，终止
//		else if (residual_vector[0] <= GlobalPara::constant::epsilon) {
//			b_writeContinue = true;
//			b_writeTecplot = true;
//			b_pause = true;
//			promptSignal = _StableReached;
//		}
//
//		// 写文件操作
//		if (b_writeContinue || b_writeTecplot) {
//			// 根据单元U更新节点ruvp
//			this->calculateNodeValue(rho_nodes, u_nodes, v_nodes, p_nodes);
//			if (b_writeContinue) {
//				//this->writeContinueFile(continueFilePath, t, istep);
//				FieldWriter::writeContinueFile(
//					istep, t, continueFilePath, nodes, elements, &(this->boundaryManager)
//				);
//			}
//			if (b_writeTecplot) {
//				//writeTecplotFile(tecplotFilePath, t);
//				FieldWriter::writeTecplotFile(t, tecplotFilePath, "title", nodes, elements, rho_nodes, u_nodes, v_nodes, p_nodes);
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
//			LogWriter::log(solveInfo);
//			break;
//		}
//	}
//
//
//}

//void FVM_2D::run_GPU() {
//	//// 初始化
//	// 尝试读取续算文件。若读取失败或者_continue==false则从头开始
//	bool startFromZero = false;
//	if (GlobalPara::basic::_continue) {
//		int flag_readContinue = readContinueFile();
//		if (flag_readContinue == -1) {
//			std::string error_msg = "Info: Fail to read previous mesh(pause_*.dat). ";
//			error_msg += "Will try to start from zero again.\n";
//			LogWriter::logAndPrint(error_msg, LogWriter::Info);
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
//				LogWriter::logAndPrint(error_msg);
//				return;//退出
//			}
//			setInitialCondition();
//		}
//		else if (type == "su2") {
//			std::string dir = FilePathManager::getInstance()->getExePath_withSlash() + "input\\";
//			int flag_readMesh = SU2MeshReader::readFile(dir + GlobalPara::basic::filename + ".su2", true);
//			if (flag_readMesh == -1) {
//				std::string error_msg = "Error: Fail to read mesh. Program will exit.\n";
//				LogWriter::logAndPrint(error_msg);
//				return;//退出
//			}
//			setInitialCondition();
//		}
//		else {
//			std::string error_msg = "Invalid mesh file type: " + type + ". Program will exit.\n";
//			LogWriter::logAndPrint(error_msg);
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
//	LogWriter::log(str);
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
//	gpuSolver2.initialize();
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
//			LogWriter::log(solveInfo);
//			break;
//		}
//	}
//
//	// 释放资源
//	gpuSolver2.finalize();
//}

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

void FVM_2D::updateAutoSaveFile(double t, int istep, int& iCurrentAutosaveFile) {
	//最多保存global::autosaveFileNum个autosave文件，再多则覆写之前的
	writeContinueFile(FilePathManager::getInstance()->getOutputDirectory() + "autosave" + std::to_string(iCurrentAutosaveFile) + "_" + GlobalPara::basic::filename + ".dat", t, istep);
	iCurrentAutosaveFile++;
	if (iCurrentAutosaveFile >= GlobalPara::output::autosaveFileNum)iCurrentAutosaveFile = 0;
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
	if (nodes.size() == 0) {
		LogWriter::logError("null nodes exception, @iniPNodeTable\n");
		exit(-1);
	}
	pNodeTable.resize(maxnodeID + 1);
	for (int in = 0; in < nodes.size(); in++) {
		pNodeTable[nodes[in].ID] = &(nodes[in]);
	}
}

void FVM_2D::iniEdges() {
	// 确保单元、节点信息已读取。这是组装edge的前提条件
	if (elements.size() == 0) {
		LogWriter::logError("null elements exception, @FVM_2D::iniEdges\n");
		exit(-1);
	}
	if (nodes.size() == 0) {
		LogWriter::logError("null nodes exception, @FVM_2D::iniEdges\n");
		exit(-1);
	}
	if (pNodeTable.size() == 0) {
		LogWriter::logError("null pNodeTable exception, @FVM_2D::iniEdges\n");
		exit(-1);
	}
	// 叉乘计算单元面积，确保节点顺序是逆时针。
	for (int ie = 0; ie < elements.size(); ie++) {
		// 叉乘计算单元面积，若面积为负，应交换某两个节点顺序，使得节点顺序是逆时针
		Element_2D& element_i = elements[ie];
		double xn[3]{}, yn[3]{};
		for (int i = 0; i < 3; i++) {
			xn[i] = this->getNodeByID(element_i.nodes[i])->x;
			yn[i] = this->getNodeByID(element_i.nodes[i])->y;
		}
		double area = 0.5 * (xn[0] * (yn[1] - yn[2]) + xn[1] * (yn[2] - yn[0]) + xn[2] * (yn[0] - yn[1]));
		if (area < 0) {
			element_i.area = -area;
			int n1 = element_i.nodes[1];
			int n2 = element_i.nodes[2];
			element_i.nodes[1] = n2;
			element_i.nodes[2] = n1;
		}
		else {
			element_i.area = area;
		}
	}
	// 组装edge
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
	if (edges.size() == 0) {
		LogWriter::logError("null edges exception, @iniPEdgeTable\n");
		exit(-1);
	}
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
	// 依赖于elements, edges,
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
	if (elements.size() == 0) {
		LogWriter::logError("null elements exception, @iniPElementTable\n");
		exit(-1);
	}
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
	不足：
	  该方式时间复杂度为O(n^2)，首先要遍历单元，然后getEdgeByNodeIDs要遍历所有edge
	  究竟该如何并行？
	*/
	if (elements.size() == 0) {
		LogWriter::logError("null elements exception, @iniElement_xy_pEdges\n");
		exit(-1);
	}
	for (int ie = 0; ie < elements.size(); ie++) {
		//已经calxy，以后不必担心。但是必须要先读取node再读取element

		Element_2D* pElement = &(elements[ie]);
		//初始化单元中心坐标
		double sum_x = 0;
		double sum_y = 0;
		for (int i = 0; i < 3; i++) {
			sum_x += this->getNodeByID(pElement->nodes[i])->x;
			sum_y += this->getNodeByID(pElement->nodes[i])->y;
		}
		pElement->x = sum_x / 3.0;
		pElement->y = sum_y / 3.0;


		pElement->pEdges[0] = this->getEdgeByNodeIDs(pElement->nodes[0], pElement->nodes[1]);
		pElement->pEdges[1] = this->getEdgeByNodeIDs(pElement->nodes[1], pElement->nodes[2]);
		pElement->pEdges[2] = this->getEdgeByNodeIDs(pElement->nodes[2], pElement->nodes[0]);


	}
	hasInitElementXY = true;//已经初始化单元中心坐标。
}

void FVM_2D::iniElement_xy_pEdges_parallel() {

	
}

void FVM_2D::iniEdges_lengths() {
	if (hasInitElementXY == false) {
		LogWriter::logAndPrintError("hasInitElementXY == false, in FVM_2D::iniEdges_lengths()\n");
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
		std::string error_msg = "ID out of range in FVM_2D::getNodeByID(), ID=" + std::to_string(ID) + "\n";
		LogWriter::logAndPrintError(error_msg);
		throw error_msg;
	}
	return pNodeTable[ID];
}

//void FVM_2D::isentropicVortex(double x, double y, double xc, double yc, double chi, double& deltau, double& deltav, double& deltaT) {
//	double xbar = x - xc;
//	double ybar = y - yc;
//	double r2 = xbar * xbar + ybar * ybar;
//	const double gamma = GlobalPara::constant::gamma;
//	const double PI = GlobalPara::constant::PI;
//	deltau = chi / 2.0 / PI * exp(0.5 * (1 - r2)) * (-ybar);
//	deltav = chi / 2.0 / PI * exp(0.5 * (1 - r2)) * xbar;
//	deltaT = -(gamma - 1) / chi / chi * 8 * gamma * PI * PI * exp(1 - r2);
//}
//
//void FVM_2D::isentropicVortex_2(double xc, double yc, double chi, const double* ruvp0) {
//	//ruvp0：均匀流参数
//	for (int ie = 0; ie < elements.size(); ie++) {
//		Element_2D& e = elements[ie];
//		double rho, u, v, p;
//		double xbar, ybar, r2, du, dv, dT;
//		const double PI = GlobalPara::constant::PI;
//		const double gamma = GlobalPara::constant::gamma;
//		xbar = e.x - xc;
//		ybar = e.y - yc;
//		r2 = xbar * xbar + ybar * ybar;
//		du = chi / 2. / PI * exp(0.5 * (1. - r2)) * (-ybar);
//		dv = chi / 2. / PI * exp(0.5 * (1. - r2)) * xbar;
//		u = ruvp0[1] + du;
//		v = ruvp0[2] + dv;
//		dT = -(gamma - 1.) * chi * chi / (8. * gamma * PI * PI) * exp(1. - r2);
//		rho = pow(ruvp0[3] + dT, 1. / (gamma - 1.));
//		p = rho * (ruvp0[3] + dT);
//		double ruvp[4]{ rho,u,v,p };
//		Math_2D::ruvp_2_U(ruvp, e.U, gamma);
//	}
//	//lambda表达式 [ capture ] ( params ) opt -> ret { body; }; http://c.biancheng.net/view/3741.html
//	//例如 auto f = [](int a) -> int { return a + 1; };
//	//auto fun_WA = [](const double* xy) {
//}
//
//void FVM_2D::writeFileHeader_isentropicVortex() {
//	std::ofstream outfile(FilePathManager::getInstance()->getOutputDirectory() + "error_isentropicVortex_" + GlobalPara::basic::filename + ".txt", std::ios::app);//追加模式
//	outfile << SystemInfo::getCurrentDateTime() << "\n";
//	outfile
//		<< "istep" << "\t"
//		<< "cpu_time[s]" << "\t"
//		<< "error_L1" << "\t"
//		<< "error_L2" << "\t"
//		<< "error_max"<< "\n";
//	outfile.close();
//}

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
				"rho ==\"NaN\", in element (x=" + std::to_string(x) + ", y=" + std::to_string(y) 
				+ ", U[0,1,2,3]=" + std::to_string(elements[ie].U[0]) + ", " + std::to_string(elements[ie].U[1])
				+ ", " + std::to_string(elements[ie].U[2]) + ", " + std::to_string(elements[ie].U[3]) + "\n";
			//break;
		}
		else if(elements[ie].U[0] < 0) {
			const double& x = elements[ie].x;
			const double& y = elements[ie].y;
			str =
				"rho < 0, in element (x=" + std::to_string(x) + ", y=" + std::to_string(y)
				+ ", U[0,1,2,3]=" + std::to_string(elements[ie].U[0]) + ", " + std::to_string(elements[ie].U[1])
				+ ", " + std::to_string(elements[ie].U[2]) + ", " + std::to_string(elements[ie].U[3]) + "\n";
		}
		LogWriter::logAndPrint(str,LogWriter::Warning);
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
