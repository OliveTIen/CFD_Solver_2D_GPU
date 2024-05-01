#ifndef ROUTINE_CONTROLLER_H
#define ROUTINE_CONTROLLER_H
#include "../gpu/datatype/DefineType.h"

class RoutineController {
	class Action {
	public:
		myfloat m_value = 1.0;
	private:
		int m_cold_time = 1000;// 冷却时间
		int m_last_used_time = 0;
	public:
		Action(){}
		Action(int _cold_time, myfloat _value) {
			initialize(_cold_time, _value);
		}
		// 初始化CD和值
		void initialize(int _cold_time, myfloat _value) {
			m_cold_time = _cold_time;
			m_value = _value;
		}
		// 更新计时器
		void update_time(int _time) {
			m_last_used_time = _time;
		}
		// 修改v。未来可考虑接受lambda函数作为参数。此处只是对输入参数做乘法
		void operate(myfloat& _value, myfloat para) const { _value *= para; }
		// 冷却结束，可使用
		bool is_CD_ready(int _time) const {
			return (_time - m_last_used_time) >= m_cold_time;
		}
		// 若CD好了，则使用技能，并重置CD，返回true。否则返回false
		bool try_operate(myfloat& _value, int _time, myfloat para) {
			if (is_CD_ready(_time)) {
				operate(_value, para);
				update_time(_time);
				return true;
			}
			return false;
		}
		myfloat get_value() { return m_value; }
	};


public:
	myfloat tolerance_max_CFL = 500;// 允许最大CFL
	myfloat tolerance_min_CFL = 0.5;// 允许最小CFL
	int increase_start_time = 5000;// 允许增大CFL的起始时间步
	int decrease_start_time = 1000;// 允许减小CFL的起始时间步

private:
	static RoutineController* p_instance;// 单例，类指针
	int m_strategy = 0;// 时间推进策略

	int m_num_of_strategy_calls = 0;// 策略已调用次数
	int m_current_time = 0;// 当前时间步
	myfloat m_current_CFL = 0.8;
	myfloat m_current_residual_rho = 1e4;
	myfloat m_record_res_rho = 10;// 残差

	Action increase;
	Action decrease;

public:
	static RoutineController* getInstance();
	// 设置时间推进策略。只被TomlFileManager调用
	void setStrategy(int s);
	int getStrategy() { return m_strategy; }
	// 应用时间推进策略，动态调整CFL数
	void applyStrategy(myfloat residual_rho, myfloat& CFL, int step);

private:
	RoutineController() {};
	// 策略：动态调整CFL
	void strategy_dynamic_CFL(myfloat& CFL);
	// 尝试用action修改参数值。若已冷却，则成功，更新action计时器，返回true
	bool tryApplyAction(Action& a, myfloat& value) {return a.try_operate(value, m_current_time, a.get_value());}
	bool tryApplyAction(Action& a, myfloat& value, myfloat para) {return a.try_operate(value, m_current_time, para);}
};


#endif // !ROUTINE_CONTROLLER_H
