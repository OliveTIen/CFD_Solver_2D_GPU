#ifndef	CBOUNDARY_DOUBLE_SHOCK_REFLECT_H
#define	CBOUNDARY_DOUBLE_SHOCK_REFLECT_H
#include <string>

/*
����˫��շ���ı߽�
�ڼ���ͨ��ʱ������˫��ձ߽磬���ø���ĺ����жϱ��Ƿ�λ�ڼ������Σ���������Ӧ��Զ���߽�
*/
class CBoundaryDoubleShockReflect {
private:
	static CBoundaryDoubleShockReflect* m_pointer;
	static const double sqrt3;// sqrt(3.0)
	double shock_x = 0.0;
	double shock_y = 0.0;
	double shock_angle_degree = 60.0;
	double shock_speed = 20.0;
	double m_physical_t = 0.0;
	double m_physical_dt = 0.0;

	CBoundaryDoubleShockReflect() {};

public:
	// 
	static CBoundaryDoubleShockReflect* getInstance();


	// ���ڲ�֪����μ��㼤���˶��ٶȣ��÷��������á�ֻ�ʺ����ڳ�ʼ������(t=0)
	bool isUpStreamOfShock_forElement(double x, double y);
	// �÷���ֻ�ʺϼ���90�Ƚǵ�����������Լ����ʱ��
	bool isUpStreamOfShock_atBoundary(double x, double y, double t_RK);
	void setVar(double x, double y, double angle) {
		shock_x = x;
		shock_y = y;
		shock_angle_degree = angle;
	}
	void set_t(double _t) { m_physical_t = _t; }
	void set_dt(double _dt) { m_physical_dt = _dt; }
	double get_t() { return m_physical_t; }
	double get_dt() { return m_physical_dt; }
	double get_t_plus_dt() { return m_physical_t + m_physical_dt; }
};

#endif