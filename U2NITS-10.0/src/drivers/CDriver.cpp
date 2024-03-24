#include <ctime>
#include "CDriver.h"
#include "../FVM_2D.h"
#include "../global/GlobalPara.h"
#include "../global/FilePathManager.h"
#include "../global/StringProcessor.h"
#include "../global/SystemInfo.h"
#include "../solvers/GPUSolver2.h"
#include "../input/CInput.h"
#include "../input/SU2MeshReader.h"
#include "../input/InpMeshReader.h"
#include "../input/TomlFileManager.h"
#include "../output/ConsolePrinter.h"
#include "../output/FieldWriter.h"
#include "../output/HistWriter.h"
#include "../output/ResidualCalculator.h"
#include "../output/LogWriter.h"


void U2NITS::CDriver::run_current() {
	CInput input;

	input.readConfig();
	input.readField();

	
	//// ���
	solve_current();


	// ֹͣ
#ifdef _WIN32
	system("pause");
#elif defined __linux__
	getchar();
#endif
}

void U2NITS::CDriver::solve_current() {
	/*
	[������]����ΪGPU����
	TODO:
	*/
	const int residual_vector_size = 4;
	const int istep_previous = GlobalPara::time::istep_previous;
	const double previous_t = GlobalPara::time::t_previous;
	const double big_value = 1e8;
	const int maxIteration = GlobalPara::output::maxIteration;
	const double T = GlobalPara::time::T;//Ŀ������ʱ�䡣���ڿ����Ƿ���ֹ����
	int nFiles = 0;//����ļ�������������ʾ
	int iCurrentAutosaveFile = 0;// �Զ������ļ���ָ�꣬����ѭ������
	double residual_vector[residual_vector_size]{ 1,1,1,1 };
	double t = previous_t;//��ǰ����ʱ�䡣������ʱ����ʼʱ�����

	std::string solveInfo;
	std::string outputPathWithSlash = FilePathManager::getInstance()->getOutputDirectory();
	std::string basicFileName = GlobalPara::basic::filename;
	clock_t start_t;
	COORD p1;
	GPU::GPUSolver2 gpuSolver2;
	HistWriter histWriter(outputPathWithSlash + basicFileName + "_hist.dat");
	bool useGPU = GlobalPara::basic::useGPU;
	if (useGPU) {
		LogWriter::writeLogAndCout("Using GPU to solve 2D PDE.\n");
	}
	else {
		LogWriter::writeLogAndCout("Using CPU to solve 2D PDE.\n");
	}

	// ��ʼ��
	start_t = clock();
	histWriter.writeHistFileHead();
	gpuSolver2.initialze();// ����host�ڴ�ҲҪ��ʼ������˲���ʡ��
	p1 = ConsolePrinter::getCursorPosition();
	//ConsolePrinter::drawProgressBar(int(t / T));

	// ����
	const int offset = 0;
	for (int istep = istep_previous; istep <= maxIteration && t <= T; istep++) {
		// �����ļ���
		char szBuffer[20];
		sprintf_s(szBuffer, _countof(szBuffer), "%04d", istep);//��istepС��4λ������0
		std::string str_istep_withBracket = "[" + std::string(szBuffer) + "]";
		std::string tecplotFilePath = outputPathWithSlash + basicFileName + str_istep_withBracket + ".dat";
		std::string continueFilePath = outputPathWithSlash + "pause_" + basicFileName + str_istep_withBracket + ".dat";
		std::string continueFilePath_nan = outputPathWithSlash + "nan_" + basicFileName + str_istep_withBracket + ".dat";

		if (useGPU) {
			//gpuSolver2.iterationGPU(t, T);
			gpuSolver2.iterationHalfGPU(t, T);
		}
		else {
			gpuSolver2.iterationTotalCPU(t, T);
		}

		// --- ��� ---
		std::cout << ".";// �����.����ʾ���ڼ��㣬û�п���
		bool b_writeContinue = false;
		bool b_writeTecplot = false;
		bool b_pause = false;//��ͣ�źţ�Ϊtrue����ֹ
		MySignal promptSignal = _NoSignal;// ��ʾ�ʷ���
		// �����������
		if (istep % GlobalPara::output::step_per_print == offset) {
			ConsolePrinter::clearDisplay(p1, ConsolePrinter::getCursorPosition());
			ConsolePrinter::setCursorPosition(p1);
			ConsolePrinter::drawProgressBar(int(double(istep) / double(maxIteration) * 100.0 + 2));//+x��Ϊ�˷�ֹͣ��99%

			double calculateTime = (double)(clock() - start_t) / CLOCKS_PER_SEC;
			double calculateSpeed = ((istep - istep_previous) / calculateTime);
			solveInfo = ConsolePrinter::assemblySolveInfo(calculateTime, istep, maxIteration, calculateSpeed, nFiles, t, T, residual_vector);
			std::cout << solveInfo;
		}
		// ��ɢ����������ѭ��
		if (residual_vector[0] > big_value) {
			// ����ֱ��break����˲���Ҫ�õ������b_writeContinue��b_writeTecplot
			// ����writeContinueFile�Ĳ�����ͬ����nan
			ConsolePrinter::printInfo(ConsolePrinter::InfoType::type_nan_detected);
			//writeTecplotFile(tecplotFilePath, t);// ��ɾ
			// ���ݵ�ԪU���½ڵ�ruvp
			gpuSolver2.updateOutputNodeField();
			FieldWriter::writeTecplotFile_GPU(t, tecplotFilePath, "title", gpuSolver2.node_host, gpuSolver2.element_host, gpuSolver2.outputNodeField);// �����2024-02-28
			// дnan�ļ�
			FieldWriter::writeContinueFile_GPU(
				istep, t, continueFilePath_nan,
				gpuSolver2.node_host, gpuSolver2.element_host, gpuSolver2.elementField_host
			);// �����2024-02-28
			break;
		}
		// �����������
		if (istep % GlobalPara::output::step_per_output_field == offset) {
			//writeTecplotFile(tecplotFilePath, t);
			b_writeTecplot = true;
			nFiles++;
			//.autosave �Զ������ļ� �����ļ�
			//��ౣ��global::autosaveFileNum��autosave�ļ����ٶ���д֮ǰ��
			std::string autosaveFilePath = outputPathWithSlash + "autosave" + std::to_string(iCurrentAutosaveFile) + "_" + GlobalPara::basic::filename + ".dat";
			iCurrentAutosaveFile++;
			if (iCurrentAutosaveFile >= GlobalPara::output::autosaveFileNum)iCurrentAutosaveFile = 0;

		}
		// ��������в�
		if (istep % GlobalPara::output::step_per_output_hist == offset) {
			//����в���������residual_vector��
			ResidualCalculator::cal_residual_GPU(gpuSolver2.element_U_old, gpuSolver2.elementField_host, ResidualCalculator::NORM_INF, residual_vector);
			histWriter.writeHistFileData(istep, residual_vector, residual_vector_size);
		}

		// --- ��ֹ ---
		// ��esc��ֹ
		if (_kbhit()) {
			if (_getch() == 27) {
				b_writeContinue = true;
				b_pause = true;
				promptSignal = _EscPressed;
			}
		}
		// �ﵽ�涨����ʱ��������������ֹ
		else if (t >= T || istep >= maxIteration) {
			b_writeContinue = true;
			b_writeTecplot = true;
			b_pause = true;
			promptSignal = _TimeReached;
		}
		// �в��㹻С����ֹ
		else if (residual_vector[0] <= GlobalPara::constant::epsilon) {
			b_writeContinue = true;
			b_writeTecplot = true;
			b_pause = true;
			promptSignal = _StableReached;
		}

		// д�ļ�����
		if (b_writeContinue || b_writeTecplot) {
			// ���ݵ�ԪU���½ڵ�ruvp
			gpuSolver2.updateOutputNodeField();
			if (b_writeContinue) {
				FieldWriter::writeContinueFile_GPU(
					istep, t, continueFilePath,
					gpuSolver2.node_host, gpuSolver2.element_host, gpuSolver2.elementField_host
				);
			}
			if (b_writeTecplot) {
				FieldWriter::writeTecplotFile_GPU(t, tecplotFilePath, "title", gpuSolver2.node_host, gpuSolver2.element_host, gpuSolver2.outputNodeField);
			}
		}
		// ��ʾ�� ����д�ļ�֮�󣬷�ֹ����
		printSignalInfo(promptSignal);
		// ��ֹ����
		if (b_pause) {
			LogWriter::writeLog(solveInfo);
			break;
		}
	}

	// �ͷ���Դ
	gpuSolver2.finalize();


}

void U2NITS::CDriver::run_2() {
	/*
���裺���Ż���
��ȡ���ã����Ʋ����� CInput.readConfig()
��ʼ����������ȡ�����ļ�/��ȡ�����ļ�+��ʼ�� CInput.readField()
�����ʼ��Ϣ��������Ϣ���߽�����ȣ�COutput.printScreen()
����в��ļ�ͷ COutput.writeHistFileHead() ���жϷŵ�������ȥ
��ʱ����ʼ CTimer.start()
��ʼ����Դ��GPU�ڴ棩CSolver.initialize()
������ʼ
	���㣨�����Ͳв CSolver.iterate() �в�������ʱ����
	���������Ϣ�������ļ�������һ�ι��λ�á�������ȣ�COutput.update()����(COutput.updateData COutput.printScreen())
	����ļ����������вCOutput.writeFile() ����COutput.writeField() COutput.writeResidual()(����CSolver.updateResidual() )
	�ж��Ƿ�����ѭ������ֹ��this->checkStop() ��Ϣ���У�
	���������Ϣ����ʾ�ʣ�COutput.updateScreen()
�ͷ���Դ CSolver.finalize()

ΪʲôҪ���̣߳�
����IO�����ȽϺ�ʱ����ʱCPU���У��������Ч��
һ��IO������ʱ�൱��2~3�ε���������������ʱ��IO�ͼ���ʱ�䶼������
���ն��߳�֪ʶ
����Ϣ���ݵķ�ʽ�����Ϣ

	*/
}

void U2NITS::CDriver::printSignalInfo(MySignal signal) {
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
