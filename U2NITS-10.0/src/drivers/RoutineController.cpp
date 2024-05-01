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
	//检查参数正确性和输出由CInput完成
	increase.initialize(1000, 1.1);
	decrease.initialize(500, 0.9);
}


void RoutineController::applyStrategy(myfloat residual_rho, myfloat& CFL, int step) {
	m_current_time = step;
	m_current_CFL = CFL;
	m_current_residual_rho = residual_rho;
	switch (m_strategy) {
	case 0:
		// 无策略
		return;
		break;

	case 1:
		// 定常，根据残差调整CFL数和时间步长
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
	修改CFL数
	将修改应用至全局变量，以及GPU的全局变量副本。由于CFL只在计算dt时用到，而GPU程序计算dt在host完成，因此无需创建副本
	*/
	// 第一次调用，记录
	if (m_num_of_strategy_calls == 0) {// m_tick在applyStrategy中更新
		m_record_res_rho = m_current_residual_rho;
		LogWriter::logAndPrint("use strategy_dynamic_CFL.\n");
		return;
	}

	bool apply_success = false;
	// 若大一个量级，则减小CFL
	if (m_current_time > decrease_start_time && m_current_residual_rho > m_record_res_rho * 10.0) {
		if (tryApplyAction(decrease, CFL, 0.75)) {
			apply_success = true;
		}

	}
	// 若小于记录的一半，则增大CFL，减小探测阈值
	if (m_current_time > increase_start_time && m_current_residual_rho < m_record_res_rho * 0.5) {
		if (tryApplyAction(increase, CFL, 1.1)) {
			m_record_res_rho *= 0.9;
			apply_success = true;
		}
	}
	// 若一直不变，则尝试增大CFL，但不改变m_res_rho
	if (m_current_time > increase_start_time) {
		if (tryApplyAction(increase, CFL, 1.05)) {
			apply_success = true;
		}
	}
	// 限制范围
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