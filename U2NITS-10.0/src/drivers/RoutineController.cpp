#include "RoutineController.h"
#include "../output/LogWriter.h"
#include <sstream>

RoutineController* RoutineController::p_instance = nullptr;
RoutineController* RoutineController::getInstance() {
	if (p_instance == nullptr) {
		p_instance = new RoutineController();
	}
	return p_instance;
}

void RoutineController::setStrategy(int s) {
	m_strategy = s;
	//��������ȷ�Ժ������CInput���
	increase.initialize(1000, 1.1);
	decrease.initialize(500, 0.9);
}


void RoutineController::applyStrategy(myfloat residual_rho, myfloat& CFL, int step) {
	m_current_time = step;
	m_current_CFL = CFL;
	m_current_residual_rho = residual_rho;
	switch (m_strategy) {
	case 0:
		// �޲���
		return;
		break;

	case 1:
		// ���������ݲв����CFL����ʱ�䲽��
		strategy_dynamic_CFL(CFL);
		break;

	default:
		return;
		break;
	}
	m_num_of_strategy_calls++;
}

void RoutineController::strategy_dynamic_CFL(myfloat& CFL) {
	/*
	�޸�CFL��
	���޸�Ӧ����ȫ�ֱ������Լ�GPU��ȫ�ֱ�������������CFLֻ�ڼ���dtʱ�õ�����GPU�������dt��host��ɣ�������贴������
	*/
	// ��һ�ε��ã���¼
	if (m_num_of_strategy_calls == 0) {// m_tick��applyStrategy�и���
		m_record_res_rho = m_current_residual_rho;
		LogWriter::logAndPrint("use strategy_dynamic_CFL.\n");
		return;
	}

	bool apply_success = false;
	// ����һ�����������СCFL
	if (m_current_time > decrease_start_time && m_current_residual_rho > m_record_res_rho * 10.0) {
		if (tryApplyAction(decrease, CFL, 0.75)) {
			apply_success = true;
		}

	}
	// ��С�ڼ�¼��һ�룬������CFL����С̽����ֵ
	if (m_current_time > increase_start_time && m_current_residual_rho < m_record_res_rho * 0.5) {
		if (tryApplyAction(increase, CFL, 1.1)) {
			m_record_res_rho *= 0.9;
			apply_success = true;
		}
	}
	// ��һֱ���䣬��������CFL�������ı�m_res_rho
	if (m_current_time > increase_start_time) {
		if (tryApplyAction(increase, CFL, 1.05)) {
			apply_success = true;
		}
	}
	// ���Ʒ�Χ
	if (CFL < tolerance_min_CFL) {
		CFL = tolerance_min_CFL;
	}
	if (CFL > tolerance_max_CFL) {
		CFL = tolerance_max_CFL;
	}

	if (apply_success) {
		std::stringstream ss;
		ss << "step = " << m_current_time << "\tresidual_rho = " << m_current_residual_rho << "\tCFL = " << CFL << "\n";
		LogWriter::log(ss.str());
	}
}