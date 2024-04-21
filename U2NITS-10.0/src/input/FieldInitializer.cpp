#include "FieldInitializer.h"
#include "../global/GlobalPara.h"
#include "../global/StringProcessor.h"
#include "../gpu/datatype/DefineType.h"
#include "../FVM_2D.h"
#include "../output/LogWriter.h"
#include "../global/CExit.h"
#include "../math/Common.h"
#include "TomlFileManager.h"
#include "../boundary_condition/CBoundaryDoubleShockReflect.h"

void FieldInitializer::setInitialAndBoundaryCondition() {

	auto& elements = FVM_2D::getInstance()->elements;
	switch (GlobalPara::initialCondition::type) {
	case 1:
	{
		printCondition("Uniform flow");
		setInitialUniform();
		break;
	}
	case 2:
	{
		printCondition("Uniform flow + isentropicVortex");
		setInitialIsentropicVortex(5, 5, 5, GlobalPara::boundaryCondition::_2D::inf::ruvp);
		break;
	}
	case 3:
	{
		printCondition("Shock Tube");
		setInitialShockTube();
		break;
	}
	case 4:
	{
		printCondition("Double Mach Reflection");
		setInitialDoubleShockReflection();
		break;
	}

	default:
	{
		std::stringstream ss;
		ss << "Error: Invalid initialCondition type " << GlobalPara::initialCondition::type << ".\n";
		LogWriter::logAndPrintError(ss.str());
		CExit::saveAndExit(GlobalPara::initialCondition::type);
		break;
	}
	}
}


void FieldInitializer::printCondition(std::string initialConditionName) {
	using namespace GlobalPara::boundaryCondition::_2D;
	std::stringstream ss;
	ss << "InitialCondition = "<< initialConditionName <<"\n";
	LogWriter::logAndPrint(ss.str());

	//ss.str("");// 清空
	//ss << "ruvp_inlet = " << StringProcessor::doubleArray_2_string(inlet::ruvp, 4) << ",";
	//ss << "ruvp_outlet = " << StringProcessor::doubleArray_2_string(outlet::ruvp, 4) << "\n";
	//ss << "ruvp_inf = " << StringProcessor::doubleArray_2_string(inf::ruvp, 4) << ",";
	//LogWriter::log(ss.str());


}

void FieldInitializer::setInitialUniform() {
	// 均匀场
	auto& elements = FVM_2D::getInstance()->elements;
	using namespace GlobalPara::boundaryCondition::_2D;
	for (int ie = 0; ie < elements.size(); ie++) {
		Math_2D::ruvp_2_U(inf::ruvp, elements[ie].U, GlobalPara::constant::gamma);
	}
}

void FieldInitializer::setInitialIsentropicVortex(double xc, double yc, double chi, const double* ruvp0) {
	// 等熵涡。ruvp0：均匀流参数
	auto& elements = FVM_2D::getInstance()->elements;
	for (int ie = 0; ie < elements.size(); ie++) {
		Element_2D& e = elements[ie];
		double rho, u, v, p;
		double xbar, ybar, r2, du, dv, dT;
		const double PI = U2NITS::Math::PI;
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

}

void FieldInitializer::setInitialShockTube() {
	// 激波管。参照论文shock_tube（https://zhuanlan.zhihu.com/p/154508317），见D:\tgl\Local\CFD\shock_tube_code
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

void FieldInitializer::setInitialDoubleShockReflection() {
	// 双马赫反射初始值参照 https://zhuanlan.zhihu.com/p/630069961
	// 激波直线由点斜式(x,y,angle)确定
	double shock_x = 0;
	double shock_y = 0;
	double shock_angle_degree = 60;
	TomlFileManager::getInstance()->getValueOrExit("initialCondition.doubleShockReflection.shock_x", shock_x);
	TomlFileManager::getInstance()->getValueOrExit("initialCondition.doubleShockReflection.shock_y", shock_y);
	TomlFileManager::getInstance()->getValueOrExit("initialCondition.doubleShockReflection.shock_angle_degree", shock_angle_degree);

	if (shock_angle_degree > 90 || shock_angle_degree < 0) {
		std::cout << "shock_angle_degree out of range [0,90]. please try again.\n";
		CExit::pressAnyKeyToExit();
	}

	CBoundaryDoubleShockReflect* pDSR = CBoundaryDoubleShockReflect::getInstance();
	pDSR->setVar(shock_x, shock_y, shock_angle_degree);

	constexpr bool isDebug = false;
	if (isDebug) {

		// 记录边界参数
		auto writeBoundaryCondition = [](double* inlet_ruvp, double* outlet_ruvp, double* inf_ruvp, const int num_ruvp)->void {
			std::string str;
			str += "BoundaryCondition:\n";
			str += "inlet::ruvp\t" + StringProcessor::doubleArray_2_string(inlet_ruvp, num_ruvp)
				+ "\noutlet::ruvp\t" + StringProcessor::doubleArray_2_string(outlet_ruvp, num_ruvp)
				+ "\ninf::ruvp\t" + StringProcessor::doubleArray_2_string(inf_ruvp, num_ruvp)
				+ "\n";
			LogWriter::logAndPrint(str);
		};
		using namespace GlobalPara::boundaryCondition::_2D;
		writeBoundaryCondition(inlet::ruvp, outlet::ruvp, inf::ruvp, 4);

		std::cout << "Debug success.\n";
		CExit::pressAnyKeyToExit();
	}

	auto& elements = FVM_2D::getInstance()->elements;
	using namespace GlobalPara::boundaryCondition::_2D;
	for (int ie = 0; ie < elements.size(); ie++) {
		Element_2D& element = elements[ie];
		CBoundaryDoubleShockReflect::getInstance();
		if (pDSR->isUpStreamOfShock_1(element.x, element.y)) {
			Math_2D::ruvp_2_U(inlet::ruvp, element.U, GlobalPara::constant::gamma);
		}
		else {
			Math_2D::ruvp_2_U(outlet::ruvp, element.U, GlobalPara::constant::gamma);
		}
	}
}
