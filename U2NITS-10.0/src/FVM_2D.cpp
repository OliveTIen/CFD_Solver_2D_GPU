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
//		//1.����
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
//		//2.�������͵����еĵ���
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
//		// sod������(���Զ�����boundaryCondition)
//		// ��������shock_tube����D:\tgl\Local\CFD\shock_tube_code
//		// ��ʵ�ϼ���Զ���߽�ʱҲ���õ�ruvp�����һ��ʼuse_ruvp=false
//		// �����Ma�ȼ���ruvp���μ�TomlFileManager::initialize_ruvp()
//		// ���ϣ�use_ruvpû��Ҫ�޸ģ���Ϊ����Ҳ�ò���
//		std::stringstream ss;
//		ss << "InitialCondition = SOD Shock Tube. Boundary condition is automatically set.\n";
//		LogWriter::logAndPrint(ss.str(), LogWriter::Info, LogWriter::Info);
//		real ruvpL[4]{ 1,0,0,1 };
//		real ruvpR[4]{ 0.125,0,0,0.1 };
//		using namespace GlobalPara::boundaryCondition::_2D;
//		// �޸ı߽�����
//		for (int i = 0; i < 4; i++) {
//			inlet::ruvp[i] = ruvpL[i];
//			outlet::ruvp[i] = ruvpR[i];
//		}
//		// ����
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
//	// �������û�б��õ�
//
//	Solver_2D solver;
//	FilePathManager* filePathManager = FilePathManager::getInstance();
//	HistWriter histWriter(filePathManager->outputFolder_path + "\\" + GlobalPara::basic::filename + "_hist.dat");
//	const int residual_vector_size = 4;
//	const double big_value = 1e8;
//	double residual_vector[residual_vector_size]{1,1,1,1};
//	histWriter.writeHistFileHead();
//
//	//�����ʽ(tecplot)�����ļ����������ÿһ֡���һ���ļ�
//	std::vector<Element_2D> elements_old;
//	std::vector<double> u_nodes;// �ڵ�u, size=nodes.size()���������ʱ���������档
//	std::vector<double> v_nodes;
//	std::vector<double> rho_nodes;
//	std::vector<double> p_nodes;
//
//	//����
//	const int maxIteration = GlobalPara::output::maxIteration;
//	double t = GlobalPara::time::t_previous;//��ǰ����ʱ�䡣
//	int istep_previous = GlobalPara::time::istep_previous;
//	const double T = GlobalPara::time::T;//Ŀ������ʱ�䡣���ڿ����Ƿ���ֹ����
//	int nFiles = 0;//����ļ�������������ʾ
//	int iCurrentAutosaveFile = 0;// �Զ������ļ���ָ�꣬����ѭ������
//	bool signal_pause = 0;//��ͣ�źţ����ڿ���
//	clock_t start_t;//��ʱ���������ڼ������ʱ��(����CPU����)
//	start_t = clock();
//	double calculateTime = 0.0;
//	double calculateSpeed = 0.0;
//	std::string solveInfo;
//	//�������������ļ�ͷ
//	if (GlobalPara::initialCondition::type == 2)writeFileHeader_isentropicVortex();
//	COORD p1 = ConsolePrinter::getCursorPosition();
//	ConsolePrinter::drawProgressBar(int(t / T));	//���ƽ�����
//	for (int istep = istep_previous; istep <= maxIteration && t <= T; istep++) {
//		elements_old = elements;
//		//����ʱ��
//		double dt = caldt(t, T);
//		t += dt;
//		//����ͨ����ʱ���ƽ�
//		solver.evolve(dt);
//		//���ݵ�ԪU���½ڵ�ruvp
//		calculateNodeValue(rho_nodes, u_nodes, v_nodes, p_nodes);
//		//������ȱ�ʾ���ڼ��㣬û�п���
//		std::cout << ".";
//		//�����������
//		if (istep % GlobalPara::output::step_per_print == 1) {
//			ConsolePrinter::clearDisplay(p1, ConsolePrinter::getCursorPosition());
//			ConsolePrinter::setCursorPosition(p1);
//			ConsolePrinter::drawProgressBar(int(t / T * 100 + 2));//+x��Ϊ�˷�ֹͣ��99%
//
//			calculateTime = (double)(clock() - start_t) / CLOCKS_PER_SEC;//
//			calculateSpeed = ((istep - istep_previous)/calculateTime);
//			solveInfo = ConsolePrinter::assemblySolveInfo(calculateTime, istep, maxIteration, calculateSpeed, nFiles, t, T, residual_vector);
//			std::cout << solveInfo;
//		}
//		//���Ƿ�ֵ
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
//		//�����������
//		if (istep % GlobalPara::output::step_per_output_field == 1) {
//			//.dat ������ʾ�ļ� tecplot��ʽ
//			char szBuffer[20];
//			sprintf_s(szBuffer, _countof(szBuffer), "%04d", istep);//��istepС��4λ������0
//			writeTecplotFile(filePathManager->getExePath_withSlash() + "output\\" + GlobalPara::basic::filename + "[" + szBuffer + "].dat",
//				t, rho_nodes, u_nodes, v_nodes, p_nodes);
//			nFiles++;
//			//.autosave �Զ������ļ� �����ļ�
//			updateAutoSaveFile(t, istep, iCurrentAutosaveFile);
//
//		}
//		//��������в�
//		if (istep % GlobalPara::output::step_per_output_hist == 1) {
//			//����в���������residual_vector��
//			ResidualCalculator::cal_residual(elements_old, elements, ResidualCalculator::NORM_INF, residual_vector);
//			histWriter.writeHistFileData(istep, residual_vector, residual_vector_size);
//
//			//�����е�����ļ� 
//			if (GlobalPara::initialCondition::type == 2) {
//				//���ʱ��
//				//end_t = clock();
//				double time_used = (double)(clock() - start_t) / CLOCKS_PER_SEC;
//				//�������ļ�
//				ResidualCalculator::cal_error_isentropicVortex(0, 0, 10, 10, 5, t, istep, time_used, GlobalPara::boundaryCondition::_2D::inf::ruvp);
//			}
//		}
//		//��esc��ֹ
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
//		//�ﵽ�涨����ʱ�䣬��ֹ
//		else if (t >= T || istep >= maxIteration) {
//			char szBuffer[20];
//			sprintf_s(szBuffer, _countof(szBuffer), "%04d", istep);
//			writeContinueFile(filePathManager->getExePath_withSlash() + "output\\pause_" + GlobalPara::basic::filename + "[" + szBuffer + "].dat", t, istep);
//			writeTecplotFile(filePathManager->getExePath_withSlash() + "output\\" + GlobalPara::basic::filename + "[" + szBuffer + "].dat",
//				t, rho_nodes, u_nodes, v_nodes, p_nodes);
//			std::cout << "Computation finished\n";
//			signal_pause = true;
//		}
//		//�в��㹻С����ֹ
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
//		//��ֹ
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
//	��GPU���븴�Ƶ��˴���Ӧ�޸ĵĵط�������
//	1.��GPU::GPUSolver2 gpuSolver2; ��Ϊ Solver_2D cpuSolver2;
//	2.ɾ����ʼ������try
//	3.���㲿�֣���gpuSolver2.iteration()��ΪcpuSolver2.evolve(dt)
//	4.�ͷ���Դ���֣�ɾ��gpuSolver2.finalize()
//	����ɲ���Git��¼
//
//	*/
//
//	const int residual_vector_size = 4;
//	const int istep_previous = GlobalPara::time::istep_previous;
//	const double previous_t = GlobalPara::time::t_previous;
//	const double big_value = 1e8;
//	const int maxIteration = GlobalPara::output::maxIteration;
//	const double T = GlobalPara::time::T;//Ŀ������ʱ�䡣���ڿ����Ƿ���ֹ����
//	int nFiles = 0;//����ļ�������������ʾ
//	int iCurrentAutosaveFile = 0;// �Զ������ļ���ָ�꣬����ѭ������
//	double residual_vector[residual_vector_size]{ 1,1,1,1 };
//	double t = previous_t;//��ǰ����ʱ�䡣������ʱ����ʼʱ�����
//	enum MySignal {
//		_NoSignal,
//		_EscPressed,
//		_TimeReached,
//		_StableReached
//	};
//	std::vector<Element_2D> elements_old;
//	std::vector<double> u_nodes;// �ڵ�u, size=nodes.size()���������ʱ���������档
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
//	// ��ʼ��
//	start_t = clock();
//	histWriter.writeHistFileHead();
//	p1 = ConsolePrinter::getCursorPosition();
//	//ConsolePrinter::drawProgressBar(int(t / T));
//
//
//	// ����
//	const int offset = 0;
//	for (int istep = istep_previous; istep <= maxIteration && t <= T; istep++) {
//		// --- ���� ---
//		// �����ļ���
//		char szBuffer[20];
//		sprintf_s(szBuffer, _countof(szBuffer), "%04d", istep);//��istepС��4λ������0
//		std::string str_istep_withBracket = "[" + std::string(szBuffer) + "]";
//		std::string tecplotFilePath = outputPathWithSlash + basicFileName + str_istep_withBracket + ".dat";
//		std::string continueFilePath = outputPathWithSlash + "pause_" + basicFileName + str_istep_withBracket + ".dat";
//		std::string continueFilePath_nan = outputPathWithSlash + "nan_" + basicFileName + str_istep_withBracket + ".dat";
//		// ���¾ɽ��
//		elements_old = this->elements;
//		// ����ʱ���ʱ�䲽��
//		double dt = this->caldt(t, T);
//		t += dt;
//
//		// --- ���� ---
//		// GPU����
//		cpuSolver2.evolve(dt);
//
//		// --- ��� ---
//		// ������ȱ�ʾ���ڼ��㣬û�п���
//		std::cout << ".";
//		bool b_writeContinue = false;
//		bool b_writeTecplot = false;
//		bool b_pause = false;//��ͣ�źţ�Ϊtrue����ֹ
//		MySignal promptSignal = _NoSignal;// ��ʾ�ʷ���
//		// �����������
//		if (istep % GlobalPara::output::step_per_print == offset) {
//			ConsolePrinter::clearDisplay(p1, ConsolePrinter::getCursorPosition());
//			ConsolePrinter::setCursorPosition(p1);
//			ConsolePrinter::drawProgressBar(int(double(istep) / double(maxIteration) * 100.0 + 2));;//+x��Ϊ�˷�ֹͣ��99%
//
//			double calculateTime = (double)(clock() - start_t) / CLOCKS_PER_SEC;
//			double calculateSpeed = ((istep - istep_previous) / calculateTime);
//			solveInfo = ConsolePrinter::assemblySolveInfo(calculateTime, istep, maxIteration, calculateSpeed, nFiles, t, T, residual_vector);
//			std::cout << solveInfo;
//		}
//		// ���Ƿ�ֵ �Ƿ�����������ѭ��
//		if (this->isNan() || residual_vector[0] > big_value) {
//			// ����ֱ��break����˲���Ҫ�õ������b_writeContinue��b_writeTecplot
//			// ����writeContinueFile�Ĳ�����ͬ����nan
//			ConsolePrinter::printInfo(ConsolePrinter::InfoType::type_nan_detected);
//			//writeTecplotFile(tecplotFilePath, t);// ��ɾ
//			// ���ݵ�ԪU���½ڵ�ruvp
//			this->calculateNodeValue(rho_nodes,u_nodes,v_nodes,p_nodes);
//			FieldWriter::writeTecplotFile(t, tecplotFilePath, "title", nodes, elements, rho_nodes, u_nodes, v_nodes, p_nodes);
//			//this->writeContinueFile(continueFilePath_nan, t, istep);// ��ɾ
//			FieldWriter::writeContinueFile(
//				istep, t, continueFilePath_nan, nodes, elements, &(this->boundaryManager)
//			);
//			break;
//		}
//		// �����������
//		if (istep % GlobalPara::output::step_per_output_field == offset) {
//			//writeTecplotFile(tecplotFilePath, t);
//			b_writeTecplot = true;
//			nFiles++;
//			//.autosave �Զ������ļ� �����ļ�
//			this->updateAutoSaveFile(t, istep, iCurrentAutosaveFile);
//
//		}
//		// ��������в�
//		if (istep % GlobalPara::output::step_per_output_hist == offset) {
//			//����в���������residual_vector��
//			ResidualCalculator::cal_residual(elements_old, elements, ResidualCalculator::NORM_INF, residual_vector);
//			histWriter.writeHistFileData(istep, residual_vector, residual_vector_size);
//		}
//
//		// --- ��ֹ ---
//		// ��esc��ֹ
//		if (_kbhit()) {
//			if (_getch() == 27) {
//				b_writeContinue = true;
//				b_pause = true;
//				promptSignal = _EscPressed;
//			}
//		}
//		// �ﵽ�涨����ʱ�䣬��ֹ
//		else if (t >= T || istep >= maxIteration) {
//			b_writeContinue = true;
//			b_writeTecplot = true;
//			b_pause = true;
//			promptSignal = _TimeReached;
//		}
//		// �в��㹻С����ֹ
//		else if (residual_vector[0] <= GlobalPara::constant::epsilon) {
//			b_writeContinue = true;
//			b_writeTecplot = true;
//			b_pause = true;
//			promptSignal = _StableReached;
//		}
//
//		// д�ļ�����
//		if (b_writeContinue || b_writeTecplot) {
//			// ���ݵ�ԪU���½ڵ�ruvp
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
//		// ��ʾ�� ����д�ļ�֮�󣬷�ֹ����
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
//		// ��ֹ����
//		if (b_pause) {
//			LogWriter::log(solveInfo);
//			break;
//		}
//	}
//
//
//}

//void FVM_2D::run_GPU() {
//	//// ��ʼ��
//	// ���Զ�ȡ�����ļ�������ȡʧ�ܻ���_continue==false���ͷ��ʼ
//	bool startFromZero = false;
//	if (GlobalPara::basic::_continue) {
//		int flag_readContinue = readContinueFile();
//		if (flag_readContinue == -1) {
//			std::string error_msg = "Info: Fail to read previous mesh(pause_*.dat). ";
//			error_msg += "Will try to start from zero again.\n";
//			LogWriter::logAndPrint(error_msg, LogWriter::Info);
//			startFromZero = true;
//			// ��ֹ����histWriter��д�ļ�ͷ
//			GlobalPara::basic::_continue = false;
//		}
//	}
//	else {
//		startFromZero = true;
//	}
//	// ��ͷ��ʼ����ȡ���񡢳�ʼ��
//	if (startFromZero) {
//		const std::string& type = GlobalPara::basic::meshFileType;
//		if (type == "inp") {
//			std::string dir = FilePathManager::getInstance()->getExePath_withSlash() + "input\\";
//			int flag_readMesh = InpMeshReader::readGmeshFile(dir + GlobalPara::basic::filename + ".inp");
//			if (flag_readMesh == -1) {
//				std::string error_msg = "Error: Fail to read mesh. Program will exit.\n";
//				LogWriter::logAndPrint(error_msg);
//				return;//�˳�
//			}
//			setInitialCondition();
//		}
//		else if (type == "su2") {
//			std::string dir = FilePathManager::getInstance()->getExePath_withSlash() + "input\\";
//			int flag_readMesh = SU2MeshReader::readFile(dir + GlobalPara::basic::filename + ".su2", true);
//			if (flag_readMesh == -1) {
//				std::string error_msg = "Error: Fail to read mesh. Program will exit.\n";
//				LogWriter::logAndPrint(error_msg);
//				return;//�˳�
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
//	// ��־��¼�߽����
//	std::string str;
//	str += "BoundaryCondition:\n";
//	str += "inlet::ruvp\t" + StringProcessor::doubleArray_2_string(GlobalPara::boundaryCondition::_2D::inlet::ruvp, 4)
//		+ "\noutlet::ruvp\t" + StringProcessor::doubleArray_2_string(GlobalPara::boundaryCondition::_2D::outlet::ruvp, 4)
//		+ "\ninf::ruvp\t" + StringProcessor::doubleArray_2_string(GlobalPara::boundaryCondition::_2D::inf::ruvp, 4)
//		+ "\n";
//	LogWriter::log(str);
//
//	//// ���
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
//	[������]����ΪGPU����
//	TODO:
//	*/
//
//	const int residual_vector_size = 4;
//	const int istep_previous = GlobalPara::time::istep_previous;
//	const double previous_t = GlobalPara::time::t_previous;
//	const double big_value = 1e8;
//	const int maxIteration = GlobalPara::output::maxIteration;
//	const double T = GlobalPara::time::T;//Ŀ������ʱ�䡣���ڿ����Ƿ���ֹ����
//	int nFiles = 0;//����ļ�������������ʾ
//	int iCurrentAutosaveFile = 0;// �Զ������ļ���ָ�꣬����ѭ������
//	double residual_vector[residual_vector_size]{ 1,1,1,1 };
//	double t = previous_t;//��ǰ����ʱ�䡣������ʱ����ʼʱ�����
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
//	// ��ʼ��
//	start_t = clock();
//	histWriter.writeHistFileHead();
//	gpuSolver2.initialize();
//	p1 = ConsolePrinter::getCursorPosition();
//	ConsolePrinter::drawProgressBar(int(t / T));
//
//
//	// ����
//	for (int istep = istep_previous; istep <= maxIteration && t <= T; istep++) {
//		// --- ���� ---
//		// �����ļ���
//		char szBuffer[20];
//		sprintf_s(szBuffer, _countof(szBuffer), "%04d", istep);//��istepС��4λ������0
//		std::string str_istep_withBracket = "[" + std::string(szBuffer) + "]";
//		std::string tecplotFilePath = outputPathWithSlash + basicFileName + str_istep_withBracket + ".dat";
//		std::string continueFilePath = outputPathWithSlash + "pause_" + basicFileName + str_istep_withBracket + ".dat";
//		std::string continueFilePath_nan = outputPathWithSlash + "nan_" + basicFileName + str_istep_withBracket + ".dat";
//		// ��host�˵ĸ��� ��U����ruvp��U_old
//		gpuSolver2.update_ruvp_Uold(this);
//
//		// ����ʱ���ʱ�䲽��
//		// ĿǰĬ������ճ�ģ����ifViscous��ifConstViscous��Ϊfalse
//		double dt = gpuSolver2.calculateDt(
//			t, T, GlobalPara::constant::gamma, Constant::Re, Constant::Pr, GlobalPara::time::CFL, Constant::R,
//			false, false
//		);
//		t += dt;
//
//		//std::cout << "gpuSolver2.iterationDevice(dt)\n";
//		// --- ���� ---
//		// GPU����
//		gpuSolver2.iterationDevice(dt);
//		gpuSolver2.device_to_host();
//		gpuSolver2.iterationHost(dt);
//		gpuSolver2.host_to_device();
//
//		// --- ��� ---
//		// ������ȱ�ʾ���ڼ��㣬û�п���
//		std::cout << ".";
//		bool b_writeContinue = false;
//		bool b_writeTecplot = false;
//		bool b_pause = false;//��ͣ�źţ�Ϊtrue����ֹ
//		MySignal promptSignal = _NoSignal;// ��ʾ�ʷ���
//		// �����������
//		if (istep % GlobalPara::output::step_per_print == 1) {
//			ConsolePrinter::clearDisplay(p1, ConsolePrinter::getCursorPosition());
//			ConsolePrinter::setCursorPosition(p1);
//			ConsolePrinter::drawProgressBar(int(t / T * 100 + 2));//+x��Ϊ�˷�ֹͣ��99%
//
//			double calculateTime = (double)(clock() - start_t) / CLOCKS_PER_SEC;
//			double calculateSpeed = ((istep - istep_previous) / calculateTime);
//			solveInfo = ConsolePrinter::assemblySolveInfo(calculateTime, istep, maxIteration, calculateSpeed, nFiles, t, T, residual_vector);
//			std::cout << solveInfo;
//		}
//		// ��ɢ����������ѭ��
//		if (residual_vector[0] > big_value) {
//			// ����ֱ��break����˲���Ҫ�õ������b_writeContinue��b_writeTecplot
//			// ����writeContinueFile�Ĳ�����ͬ����nan
//			ConsolePrinter::printInfo(ConsolePrinter::InfoType::type_nan_detected);
//			//writeTecplotFile(tecplotFilePath, t);// ��ɾ
//			// ���ݵ�ԪU���½ڵ�ruvp
//			gpuSolver2.updateOutputNodeField();
//			tecplotWriter.setFilePath(tecplotFilePath);
//			tecplotWriter.writeTecplotFile_GPU(t, "title", gpuSolver2.node_host,gpuSolver2.element_host,gpuSolver2.outputNodeField);// �����2024-02-28
//			// дnan�ļ�
//			continueWriter.writeContinueFile_GPU(
//				istep, t, continueFilePath_nan, 
//				gpuSolver2.node_host, gpuSolver2.element_host, gpuSolver2.elementField_host, &(this->boundaryManager)
//			);// �����2024-02-28
//			break;
//		}
//		// �����������
//		if (istep % GlobalPara::output::step_per_output_field == 1) {
//			//writeTecplotFile(tecplotFilePath, t);
//			b_writeTecplot = true;
//			nFiles++;
//			//.autosave �Զ������ļ� �����ļ�
//			//��ౣ��global::autosaveFileNum��autosave�ļ����ٶ���д֮ǰ��
//			std::string autosaveFilePath = outputPathWithSlash + "autosave" + std::to_string(iCurrentAutosaveFile) + "_" + GlobalPara::basic::filename + ".dat";
//			continueWriter.writeContinueFile_GPU(
//				istep, t, autosaveFilePath,
//				gpuSolver2.node_host, gpuSolver2.element_host, gpuSolver2.elementField_host, &(this->boundaryManager)
//			);// �����2024-02-28
//			iCurrentAutosaveFile++;
//			if (iCurrentAutosaveFile >= GlobalPara::output::autosaveFileNum)iCurrentAutosaveFile = 0;
//
//		}
//		// ��������в�
//		if (istep % GlobalPara::output::step_per_output_hist == 1) {
//			//����в���������residual_vector��
//			ResidualCalculator::cal_residual_GPU(gpuSolver2.element_U_old, gpuSolver2.elementField_host, ResidualCalculator::NORM_INF, residual_vector);
//			histWriter.writeHistFileData(istep, residual_vector, residual_vector_size);
//		}
//
//		// --- ��ֹ ---
//		// ��esc��ֹ
//		if (_kbhit()) {
//			if (_getch() == 27) {
//				b_writeContinue = true;
//				b_pause = true;
//				promptSignal = _EscPressed;
//			}
//		}
//		// �ﵽ�涨����ʱ��������������ֹ
//		else if (t >= T || istep >= maxIteration) {
//			b_writeContinue = true;
//			b_writeTecplot = true;
//			b_pause = true;
//			promptSignal = _TimeReached;
//		}
//		// �в��㹻С����ֹ
//		else if (residual_vector[0] <= Constant::epsilon) {
//			b_writeContinue = true;
//			b_writeTecplot = true;
//			b_pause = true;
//			promptSignal = _StableReached;
//		}
//
//		// д�ļ�����
//		if (b_writeContinue || b_writeTecplot) {
//			// ���ݵ�ԪU���½ڵ�ruvp
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
//		// ��ʾ�� ����д�ļ�֮�󣬷�ֹ����
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
//		// ��ֹ����
//		if (b_pause) {
//			LogWriter::log(solveInfo);
//			break;
//		}
//	}
//
//	// �ͷ���Դ
//	gpuSolver2.finalize();
//}

//double FVM_2D::caldt_GPU(double t, double T, void* gpuSolver2) {
//	REAL Re = 1e8;
//	REAL Pr = 0.73;
//	// GlobalPara::physicsModel::equation;
//	// ĿǰĬ������ճ�ģ����ifViscous��ifConstViscous��Ϊfalse
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
		//һ��̫���ᵼ�¶�ȡʱbuffer���Ȳ������Ӷ���������Ԥ�ϵĽ��(��ѭ������ȡʧ�ܵ�)
		int nEdge = (int)boundaryManager.boundaries[ib].pEdges.size();//�������' '����'\n'
		bool check = 0;//����
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
		if(!check)outfile << std::endl;//���ϴ������' '
	}
	outfile.close();

}

void FVM_2D::updateAutoSaveFile(double t, int istep, int& iCurrentAutosaveFile) {
	//��ౣ��global::autosaveFileNum��autosave�ļ����ٶ���д֮ǰ��
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
	int* node_neighborElement_num;// ��¼ÿ���ڵ���ھӵ�Ԫ����
	// ������Դ ��ʼ��
	node_neighborElement_num = new int[num_node]();// С�����Զ���ʼ��Ϊ0
	for (int i = 0; i < 4; i++) {
		memset(node_ruvp[i], 0, num_node * sizeof(REAL));
	}

	// ����Ԫֵ�ӵ����ӽڵ�ֵ��
	for (int iElement = 0; iElement < num_element; iElement++) {
		// ͳ��ÿ���ڵ���ھӵ�Ԫ���������޸Ľڵ�ֵ
		for (int jNode = 0; jNode < nNodePerElement; jNode++) {
			// ��ȡ��Ԫ�ڵ�ID
			int GPUID_of_node = gpuSolver->element_host.nodes[jNode][iElement];// GPUID of node 0
			// nodeID[3]=-1ʱ�������ýڵ�
			if (GPUID_of_node < 0 || GPUID_of_node >= num_node)continue;// ��������ѭ������������ѭ����
			// ��ID��Ӧ���ھӵ�Ԫ����+1
			node_neighborElement_num[GPUID_of_node]++;
			// ��ID��Ӧ�Ľڵ�����ֵ���ϵ�Ԫֵ
			for (int kValue = 0; kValue < nValuePerNode; kValue++) {
				node_ruvp[kValue][GPUID_of_node] += gpuSolver->element_vruvp[kValue][iElement];
			}
		}
	}

	// �ڵ�ֵ�����ھӵ�Ԫ�����õ�ƽ��ֵ����Ϊ�ڵ�ruvpֵ
	for (int iNode = 0; iNode < num_node; iNode++) {
		// Ϊ�˱������0����ĸ����0������
		if (node_neighborElement_num[iNode] == 0)continue;
		// node_ruvp�����ھӵ�Ԫ�����õ�ƽ��ֵ
		for (int kValue = 0; kValue < nValuePerNode; kValue++) {
			node_ruvp[kValue][iNode] /= node_neighborElement_num[iNode];
		}
	}

	// �ͷ���Դ
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
	// ȷ����Ԫ���ڵ���Ϣ�Ѷ�ȡ��������װedge��ǰ������
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
	// ��˼��㵥Ԫ�����ȷ���ڵ�˳������ʱ�롣
	for (int ie = 0; ie < elements.size(); ie++) {
		// ��˼��㵥Ԫ����������Ϊ����Ӧ����ĳ�����ڵ�˳��ʹ�ýڵ�˳������ʱ��
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
	// ��װedge
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
	// ǰ��������elements, elements.nodes, pNodeTable
	for (int ie = 0; ie < elements.size(); ie++) {
		for (int jn = 0; jn < 3; jn++) {
			Node_2D* pN = getNodeByID(elements[ie].nodes[jn]);
			pN->neighborElements.push_back(&(elements[ie]));
		}
	}
}

Edge_2D* FVM_2D::getEdgeByNodeIDs(int n0, int n1) {
	// ��edges�и���nodeIDs��ѯedge
	//n0, n1��ʾ�ڵ�ID
	if(edges.size()==0)return nullptr;
	for (int i = 0; i < edges.size(); i++) {
		if (edges[i].nodes[0] == n0 && edges[i].nodes[1] == n1)return &(edges[i]);
		if (edges[i].nodes[0] == n1 && edges[i].nodes[1] == n0)return &(edges[i]);
	}
	return nullptr;
}

void FVM_2D::iniEdges_registerSingle(int n0, int n1, Element_2D* pE) {
	// ������elements, edges,
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
		�������裺һ��edge����Ѿ�����ĳ��element A����A�Ǹ�edge����element
		һ��edge���ͬʱ����2��element
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
	��ʼ����Ԫ��x,y,pEdges�������Ժ����
	ǰ��������
	  ��Ҫ��elements����
	  ��Ҫ��֪element��nodeIDs
	  node�������ѳ�ʼ��
	  ��edges���飬��edge��nodeIDs�ѳ�ʼ��
	������
	  �����α߽磬pEdgesֻ��3���������ı��Σ�pEdges��node��Ҫ����ȷ��
	���㣺
	  �÷�ʽʱ�临�Ӷ�ΪO(n^2)������Ҫ������Ԫ��Ȼ��getEdgeByNodeIDsҪ��������edge
	  ��������β��У�
	*/
	if (elements.size() == 0) {
		LogWriter::logError("null elements exception, @iniElement_xy_pEdges\n");
		exit(-1);
	}
	for (int ie = 0; ie < elements.size(); ie++) {
		//�Ѿ�calxy���Ժ󲻱ص��ġ����Ǳ���Ҫ�ȶ�ȡnode�ٶ�ȡelement

		Element_2D* pElement = &(elements[ie]);
		//��ʼ����Ԫ��������
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
	hasInitElementXY = true;//�Ѿ���ʼ����Ԫ�������ꡣ
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
	// ��ֹԽ��
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
//	//ruvp0������������
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
//	//lambda���ʽ [ capture ] ( params ) opt -> ret { body; }; http://c.biancheng.net/view/3741.html
//	//���� auto f = [](int a) -> int { return a + 1; };
//	//auto fun_WA = [](const double* xy) {
//}
//
//void FVM_2D::writeFileHeader_isentropicVortex() {
//	std::ofstream outfile(FilePathManager::getInstance()->getOutputDirectory() + "error_isentropicVortex_" + GlobalPara::basic::filename + ".txt", std::ios::app);//׷��ģʽ
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
	// ���ÿ��Ԫ�ص��غ�������������쳣ֵ
	bool is_nan = 0;
	for (int ie = 0; ie < elements.size(); ie++) {
		
		std::string str;
		if (isnan(elements[ie].U[0])) {// ������õ���ϵͳ��isnan����
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
//	GPU����ʱ�������NaNֵ��������ﲻ��Ҫ���
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
