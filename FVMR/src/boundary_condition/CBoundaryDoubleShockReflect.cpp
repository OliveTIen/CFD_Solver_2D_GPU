#include "CBoundaryDoubleShockReflect.h"
#include "../math/CommonValue.h"


CBoundaryDoubleShockReflect* CBoundaryDoubleShockReflect::m_pointer = nullptr;

const double CBoundaryDoubleShockReflect::sqrt3 = sqrt(3.0);

CBoundaryDoubleShockReflect* CBoundaryDoubleShockReflect::getInstance() {
	if (m_pointer == nullptr) {
		m_pointer = new CBoundaryDoubleShockReflect();
	}
	return m_pointer;
}


bool CBoundaryDoubleShockReflect::isUpStreamOfShock_forElement(double x, double y) {
	/*
	功能：更新双马赫反射的上边界
	输入：坐标
	输出：是否位于激波上游
	*/
	// 斜激波法向量
	double angle_rad = shock_angle_degree * U2NITS::Math::PI / 180.0;
	double norm_x = -sin(angle_rad);
	double norm_y = cos(angle_rad);
	// 待求点矢径
	double r_x = x - shock_x;
	double r_y = y - shock_y;
	// 矢径与法向量点积大于0，则为上游
	return (r_x * norm_x + r_y * norm_y > 0);

}


bool CBoundaryDoubleShockReflect::isUpStreamOfShock_atBoundary(double x, double y) {
	/*
	https://zhuanlan.zhihu.com/p/630069961
	该方法仅适合边界
	xy边界元坐标，t物理时间，shock_x要求是激波
	本来公式是 x < 1/6 + (1 + 20 * t)/sqrt3，其中1/6是激波与交点的初始横坐标，1是
	域高度
	此处要求y是区域顶部的y坐标
	*/
	double t_RK = get_t_plus_dt();
	double right = shock_x + (y + shock_speed * t_RK) / sqrt3;
	if (x < right) return true;
	else return false;
}
