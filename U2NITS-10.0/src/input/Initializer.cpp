#include "Initializer.h"
#include "../global/GlobalPara.h"
#include "../global/StringProcessor.h"
#include "../gpu/datatype/Define.h"
#include "../FVM_2D.h"
#include "../output/LogWriter.h"

void Initializer::setInitialAndBoundaryCondition() {

	auto& elements = FVM_2D::getInstance()->elements;
	switch (GlobalPara::initialCondition::type) {
	case 1:
	{
		//1.����
		printCondition("Uniform flow");
		setInitialUniform();
		break;
	}
	

	case 2:
	{
		//2.�������͵����еĵ���
		printCondition("Uniform flow + isentropicVortex");
		using namespace GlobalPara::boundaryCondition::_2D;
		setInitialIsentropicVortex(5, 5, 5, inf::ruvp);
		break;
	}

	case 3:
	{
		printCondition("Shock Tube");
		setInitialShockTube();
		break;
	}

	case 1001:
	{
		// sod������(���Զ�����boundaryCondition)
		// ��������shock_tube��https://zhuanlan.zhihu.com/p/154508317������D:\tgl\Local\CFD\shock_tube_code
		// ��ʵ�ϼ���Զ���߽�ʱҲ���õ�ruvp�����һ��ʼuse_ruvp=false
		// �����Ma�ȼ���ruvp���μ�TomlFileManager::initialize_ruvp()
		// ���ϣ�use_ruvpû��Ҫ�޸ģ���Ϊ����Ҳ�ò���
		// ��ȷ�����Matlab����
		printCondition("SOD Shock Tube");
		std::cout << "Boundary condition is automatically set.\n";
		using namespace GlobalPara::boundaryCondition::_2D;
		// �޸ı߽�����
		real ruvpL[4]{ 1,0,0,1 };
		real ruvpR[4]{ 0.125,0,0,0.1 };
		for (int i = 0; i < 4; i++) {
			inlet::ruvp[i] = ruvpL[i];
			outlet::ruvp[i] = ruvpR[i];
		}
		// ����
		setInitialShockTube();
		break;
	}
	case 1002:
	{
		printCondition("Lax Shock Tube");
		std::cout << "Boundary condition is automatically set.\n";
		using namespace GlobalPara::boundaryCondition::_2D;
		// �޸ı߽�����
		real ruvpL[4]{ 0.445,0.698,0,3.528 };
		real ruvpR[4]{ 0.5,0,0,0.571 };
		for (int i = 0; i < 4; i++) {
			inlet::ruvp[i] = ruvpL[i];
			outlet::ruvp[i] = ruvpR[i];
		}
		// ����
		setInitialShockTube();
		break;
	}
	case 1003:
	{
		printCondition("Double Mach Reflection");
		// ˫��շ����ʼֵ���� https://zhuanlan.zhihu.com/p/630069961
	}

	default:
	{
		std::stringstream ss;
		ss << "Error: Invalid initialCondition type " << GlobalPara::initialCondition::type << ".\n";
		LogWriter::writeLogAndCout(ss.str(), LogWriter::Error, LogWriter::Error);
		exit(GlobalPara::initialCondition::type);
	}
	}
}

void Initializer::printCondition(std::string initialConditionName) {
	using namespace GlobalPara::boundaryCondition::_2D;
	std::stringstream ss;
	ss << "InitialCondition = "<< initialConditionName <<"\n";
	LogWriter::writeLogAndCout(ss.str(), LogWriter::Info, LogWriter::Info);
	ss.str("");// ���
	ss << "ruvp_inf = " << StringProcessor::doubleArray_2_string(inf::ruvp, 4) << ",";
	ss << "ruvp_inlet = " << StringProcessor::doubleArray_2_string(inlet::ruvp, 4) << ",";
	ss << "ruvp_outlet = " << StringProcessor::doubleArray_2_string(outlet::ruvp, 4) << "\n";
	LogWriter::writeLog(ss.str(), LogWriter::Info);


}

void Initializer::setInitialUniform() {
	auto& elements = FVM_2D::getInstance()->elements;
	using namespace GlobalPara::boundaryCondition::_2D;
	for (int ie = 0; ie < elements.size(); ie++) {
		Math_2D::ruvp_2_U(inf::ruvp, elements[ie].U, GlobalPara::constant::gamma);
	}
}

void Initializer::setInitialIsentropicVortex(double xc, double yc, double chi, const double* ruvp0) {
	//ruvp0������������
	auto& elements = FVM_2D::getInstance()->elements;
	for (int ie = 0; ie < elements.size(); ie++) {
		Element_2D& e = elements[ie];
		double rho, u, v, p;
		double xbar, ybar, r2, du, dv, dT;
		const double PI = GlobalPara::constant::PI;
		const double gamma = GlobalPara::constant::gamma;
		xbar = e.x - xc;
		ybar = e.y - yc;
		r2 = xbar * xbar + ybar * ybar;
		du = chi / 2. / PI * exp(0.5 * (1. - r2)) * (-ybar);
		dv = chi / 2. / PI * exp(0.5 * (1. - r2)) * xbar;
		u = ruvp0[1] + du;
		v = ruvp0[2] + dv;
		dT = -(gamma - 1.) * chi * chi / (8. * gamma * PI * PI) * exp(1. - r2);
		rho = pow(ruvp0[3] + dT, 1. / (gamma - 1.));
		p = rho * (ruvp0[3] + dT);
		double ruvp[4]{ rho,u,v,p };
		Math_2D::ruvp_2_U(ruvp, e.U, gamma);
	}
	//lambda���ʽ [ capture ] ( params ) opt -> ret { body; }; http://c.biancheng.net/view/3741.html
	//���� auto f = [](int a) -> int { return a + 1; };
	//auto fun_WA = [](const double* xy) {
}

void Initializer::setInitialShockTube() {
	auto& elements = FVM_2D::getInstance()->elements;
	using namespace GlobalPara::boundaryCondition::_2D;
	for (int ie = 0; ie < elements.size(); ie++) {
		if (elements[ie].x < 0) {
			Math_2D::ruvp_2_U(inlet::ruvp, elements[ie].U, GlobalPara::constant::gamma);
		}
		else {
			Math_2D::ruvp_2_U(outlet::ruvp, elements[ie].U, GlobalPara::constant::gamma);
		}
	}
}
