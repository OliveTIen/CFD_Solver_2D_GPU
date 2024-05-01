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
	// 
	static CBoundaryDoubleShockReflect* getInstance();


	// 由于不知道如何计算激波运动速度，该方法被弃用。只适合用于初始化初场(t=0)
	bool isUpStreamOfShock_forElement(double x, double y);
	// 该方法只适合计算90度角的情况，但可以计算各时间
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