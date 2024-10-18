#ifndef	CBOUNDARY_DOUBLE_SHOCK_REFLECT_H
#define	CBOUNDARY_DOUBLE_SHOCK_REFLECT_H
#include <string>

/*
更新双马赫反射的边界
在计算通量时，对于双马赫边界，会用该类的函数判断边是否位于激波上游，并设置相应的远场边界
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
	// 获取类实例
	static CBoundaryDoubleShockReflect* getInstance();
	// 由于不知道如何计算激波运动速度，该方法被弃用。只适合用于初始化初场(t=0)
	bool isUpStreamOfShock_forElement(double x, double y);
	// 该方法只适合计算90度角的情况，但可以计算各时间
	bool isUpStreamOfShock_atBoundary(double x, double y);
	// 设置成员变量shock_x,shock_y,shock_angle_degree的值
	void set_shock_x_y_angle(double x, double y, double angle) {
		shock_x = x;
		shock_y = y;
		shock_angle_degree = angle;
	}
	// 设置成员变量t的值
	void set_t(double _t) { m_physical_t = _t; }
	// 设置成员变量dt的值
	void set_dt(double _dt) { m_physical_dt = _dt; }
	double get_t() { return m_physical_t; }
	double get_dt() { return m_physical_dt; }
	double get_t_plus_dt() { return m_physical_t + m_physical_dt; }
	double get_shock_x() { return shock_x; }
	double get_shock_y() { return shock_y; }
	double get_shock_speed() { return shock_speed; }
	double get_sqrt3() { return sqrt3; }
};

#endif