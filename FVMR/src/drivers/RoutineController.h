#ifndef ROUTINE_CONTROLLER_H
#define ROUTINE_CONTROLLER_H
#include "../gpu/datatype/DefineType.h"

class RoutineController {
	class Action {
	public:
		myfloat m_value = 1.0;
	private:
		int m_cold_time = 1000;
		int m_last_used_time = 0;
	public:
		Action(){}
		Action(int _cold_time, myfloat _value) {
			initialize(_cold_time, _value);
		}
		void initialize(int _cold_time, myfloat _value) {
			m_cold_time = _cold_time;
			m_value = _value;
		}
		void update_time(int _time) {
			m_last_used_time = _time;
		}
		// modify value �޸�v��δ���ɿ��ǽ���lambda������Ϊ�������˴�ֻ�Ƕ�����������˷�
		void operate(myfloat& _value, myfloat para) const { _value *= para; }
		// is cold time ready ��ȴ��������ʹ��
		bool is_CD_ready(int _time) const {
			return (_time - m_last_used_time) >= m_cold_time;
		}
		// ��CD���ˣ���ʹ�ü��ܣ�������CD������true�����򷵻�false
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
	myfloat tolerance_max_CFL = 500;// �������CFL(���淢��CFL������������˸ò���δʹ��)
	myfloat tolerance_min_CFL = 0.5;// ������СCFL
	int increase_start_time = 5000;// ��������CFL����ʼʱ�䲽
	int decrease_start_time = 1000;// �����СCFL����ʼʱ�䲽

private:
	static RoutineController* p_instance;// ��������ָ��
	int m_strategy = 0;// ʱ���ƽ�����

	int m_num_of_strategy_calls = 0;// �����ѵ��ô���
	int m_current_time = 0;// ��ǰʱ�䲽
	myfloat m_current_CFL = 0.8;
	myfloat m_current_residual_rho = 1e4;
	myfloat m_record_res_rho = 10;// �в�

	Action increase;
	Action decrease;

public:
	static RoutineController* getInstance();
	// ����ʱ���ƽ����ԡ�ֻ��TomlFileManager����
	void setStrategy(int s);
	int getStrategy() { return m_strategy; }
	// Ӧ��ʱ���ƽ����ԣ���̬����CFL��
	void applyStrategy(myfloat residual_rho, myfloat& CFL, int step, bool& apply_success);

private:
	RoutineController() {};
	// ���ԣ���̬����CFL
	void strategy_dynamic_CFL(myfloat& CFL, bool& apply_success);
	// ������action�޸Ĳ���ֵ��������ȴ����ɹ�������action��ʱ��������true
	bool tryApplyAction(Action& a, myfloat& value) {return a.try_operate(value, m_current_time, a.get_value());}
	bool tryApplyAction(Action& a, myfloat& value, myfloat para) {return a.try_operate(value, m_current_time, para);}
};


#endif // !ROUTINE_CONTROLLER_H
