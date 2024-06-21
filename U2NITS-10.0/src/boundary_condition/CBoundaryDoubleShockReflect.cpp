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
	���ܣ�����˫��շ�����ϱ߽�
	���룺����
	������Ƿ�λ�ڼ�������
	*/
	// б����������
	double angle_rad = shock_angle_degree * U2NITS::Math::PI / 180.0;
	double norm_x = -sin(angle_rad);
	double norm_y = cos(angle_rad);
	// �����ʸ��
	double r_x = x - shock_x;
	double r_y = y - shock_y;
	// ʸ���뷨�����������0����Ϊ����
	return (r_x * norm_x + r_y * norm_y > 0);

}


bool CBoundaryDoubleShockReflect::isUpStreamOfShock_atBoundary(double x, double y) {
	/*
	https://zhuanlan.zhihu.com/p/630069961
	�÷������ʺϱ߽�
	xy�߽�Ԫ���꣬t����ʱ�䣬shock_xҪ���Ǽ���
	������ʽ�� x < 1/6 + (1 + 20 * t)/sqrt3������1/6�Ǽ����뽻��ĳ�ʼ�����꣬1��
	��߶�
	�˴�Ҫ��y�����򶥲���y����
	*/
	double t_RK = get_t_plus_dt();
	double right = shock_x + (y + shock_speed * t_RK) / sqrt3;
	if (x < right) return true;
	else return false;
}
