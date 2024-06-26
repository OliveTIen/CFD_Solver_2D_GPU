#include "drivers/CDriver.h"
using namespace std;

int main() {
	U2NITS::CDriver::getInstance()->run_20240517();
	return 0;
}


//命名规范：
//ini_x初始化变量的值
//x_to_x计算并输出。输出方式可以为return、&
//set_x改变成员变量的值
//get_x输出成员变量的值
//cal_x计算并改变成员变量的值，不输出



