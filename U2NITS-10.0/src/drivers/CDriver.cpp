#include "CDriver.h"
#include "../FVM_2D.h"
#include "../solvers/GPUSolver2.h"
#include "../input/CInput.h"
#include "../output/COutput.h"


void U2NITS::CDriver::run_current() {

/*
[������]����ΪGPU����

*/

	CInput input;
	COutput out;
	GPU::GPUSolver2 solver;

	out.console.printWelcome_2();
	input.readConfig();
	input.readField_old();

	solver.initializeHost();// �����ڴ棬��FVM_2D�����ݳ�ʼ��
	if (GlobalPara::basic::useGPU) {
		solver.initializeDevice();
	}

	const int istep_previous = GlobalPara::time::istep_previous;
	const int maxIteration = GlobalPara::output::maxIteration;
	double t = GlobalPara::time::t_previous;// ��ǰ����ʱ��
	const double T = GlobalPara::time::T;//Ŀ������ʱ�䡣���ڿ����Ƿ���ֹ����
	out.hist.setFilePath(out.getDir() + GlobalPara::basic::filename + "_hist.dat");



	// ����
	out.timer.start();// ��ʱ
	out.console.setClearStartPosition();
	for (int istep = istep_previous; istep <= maxIteration && t <= T; istep++) {
		// �����ļ���
		out.updateFileName(istep);

		if (GlobalPara::basic::useGPU) {
			// solver.iterationGPU(t, T);
			solver.iterationHalfGPU(t, T);
		}
		else {
			solver.iterationTotalCPU(t, T);
		}
		std::cout << ".";// �����.����ʾ���ڼ��㣬û�п���

		SignalPack sp = emitSignalPack(istep, maxIteration, t, T, solver.residualVector[0]);

		// ��ӡ����
		if (sp.b_print) {
			out.console.clear();
			out.console.drawProgressBar(int(double(istep) / double(maxIteration) * 100.0 + 2));//+x��Ϊ�˷�ֹͣ��99%

			double calculateTime = out.timer.getIime();
			double calculateSpeed = ((istep - istep_previous) / calculateTime);
			std::string solveInfo = out.console.assemblySolveInfo(calculateTime, istep, maxIteration, calculateSpeed, out.field.getNumTecplotFileWritten(), t, T, solver.residualVector);
			out.log.print(solveInfo);
			if (sp.pauseSignal != _NoSignal) {
				out.log.log(solveInfo);
			}
		}
		// д�ļ�����
		if (sp.b_writeContinue || sp.b_writeTecplot) {
			solver.updateOutputNodeField();// ���ݵ�ԪU���½ڵ�ruvp
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
			//����в���������residual_vector��
			solver.updateResidual();
			out.hist.writeHistFile_2(istep, solver.residualVector, solver.residualVectorSize);
		}
		// ��ֹ����
		if (sp.pauseSignal != _NoSignal) {
			out.log.logAndPrintSignalInfo(sp.pauseSignal);// ��ʾ�� ����д�ļ�֮�󣬷�ֹ����
			break;// ���ɷŽ����������У���ΪҪ��ֹѭ��
		}
	}

	// �ͷ���Դ
	solver.freeHostMemory();
	if (GlobalPara::basic::useGPU) {
		solver.freeDeviceMemory();
	}

	// ֹͣ
#ifdef _WIN32
	system("pause");
#elif defined __linux__
	getchar();
#endif

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
	����ļ����������вCOutput.writeFile() ����COutput.writeField() COutput.writeHist()(����CSolver.updateResidual() )
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

U2NITS::CDriver::SignalPack U2NITS::CDriver::emitSignalPack(int istep, int maxIteration, real t, real T, real residualRho) {
	SignalPack s;

	// �����ź�
	const int offset = 0;
	const double big_value = 1e8;

	//PauseSignal pauseSignal = _NoSignal;// ��ͣ�ź� ��ʾ�ʷ���
	if (istep % GlobalPara::output::step_per_print == offset) {// �����������
		s.b_print = true;
	}
	if (istep % GlobalPara::output::step_per_output_field == offset) {// �����������
		s.b_writeTecplot = true;
	}
	if (istep % GlobalPara::output::step_per_output_hist == offset) {// ��������в�
		s.b_writeHist = true;
	}
	if (residualRho > big_value) {// �в������ֹ
		s.b_writeContinue = true;
		s.b_writeTecplot = true;
		s.b_nanDetected = true;
		s.pauseSignal = _NanDetected;
	}
	else if (_kbhit()) {// ��esc��ֹ
		if (_getch() == 27) {
			s.b_writeContinue = true;
			s.pauseSignal = _EscPressed;
		}
	}
	else if (t >= T || istep >= maxIteration) {// �ﵽ�涨����ʱ��������������ֹ
		s.b_writeContinue = true;
		s.b_writeTecplot = true;
		s.pauseSignal = _TimeReached;
	}
	else if (residualRho <= GlobalPara::constant::epsilon) {// �в��㹻С����ֹ
		s.b_writeContinue = true;
		s.b_writeTecplot = true;
		s.pauseSignal = _StableReached;
	}

	return s;
}

void U2NITS::CDriver::onSignalPack(const SignalPack& sp) {

}
