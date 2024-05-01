#include "FieldInitializer.h"
#include "../global/GlobalPara.h"
#include "../global/StringProcessor.h"
#include "../gpu/datatype/DefineType.h"
#include "../FVM_2D.h"
#include "../output/LogWriter.h"
#include "../global/CExit.h"
#include "../math/PhysicalKernel.h"
#include "TomlFileManager.h"
#include "../boundary_condition/CBoundaryDoubleShockReflect.h"

FieldInitializer* FieldInitializer::p_instance = nullptr;

FieldInitializer* FieldInitializer::getInstance() {
	if (p_instance == nullptr) {
		p_instance = new FieldInitializer();
	}
	return p_instance;
}

enum InitialConditionType {
	type_none,
	type_uniform_flow,         // 均匀场
	type_isentropic_vortex,    // 等熵涡
	type_shock_tube,           // 激波管
	type_double_mach_reflection// 双马赫反射
};

void FieldInitializer::setInitialAndBoundaryCondition() {
	
	auto& elements = FVM_2D::getInstance()->elements;
	switch (m_initial_type) {
	case type_uniform_flow:
	{
		printCondition("Uniform flow");
		setInitialUniform();
		break;
	}
	case type_isentropic_vortex:
	{
		printCondition("Uniform flow + isentropicVortex");
		setInitialIsentropicVortex();
		break;
	}
	case type_shock_tube:
	{
		printCondition("Shock Tube");
		setInitialShockTube();
		break;
	}
	case type_double_mach_reflection:
	{
		printCondition("Double Mach Reflection");
		setInitialDoubleShockReflection();
		break;
	}

	default:
	{
		std::stringstream ss;
		ss << "Error: Invalid initialCondition type " << m_initial_type << ".\n";
		LogWriter::logAndPrintError(ss.str());
		CExit::saveAndExit(-1);
		break;
	}
	}
}

void FieldInitializer::initialize_using_config(void* tomlFileManager) {
	/*
	用TomlFileManager提供的API进行初始化
	对于双马赫反射，还要求读取激波位置、角度参数
	*/
	TomlFileManager* t = (TomlFileManager*)tomlFileManager;
	int initialCondition_type = type_uniform_flow;// 初始化方式
	t->getValue("initialCondition.type", initialCondition_type);
	FieldInitializer::getInstance()->set_initial_type(initialCondition_type);
	if (initialCondition_type == type_double_mach_reflection) {
		double shock_x = 0;
		double shock_y = 0;
		double shock_angle_degree = 60;
		t->getValueOrExit("initialCondition.doubleShockReflection.shock_x", shock_x);
		t->getValueOrExit("initialCondition.doubleShockReflection.shock_y", shock_y);
		t->getValueOrExit("initialCondition.doubleShockReflection.shock_angle_degree", shock_angle_degree);

		if (shock_angle_degree > 90 || shock_angle_degree < 0) {
			std::cout << "shock_angle_degree out of range [0,90]. please try again.\n";
			CExit::pressAnyKeyToExit();
		}

		CBoundaryDoubleShockReflect::getInstance()->setVar(shock_x, shock_y, shock_angle_degree);
	}
}


void FieldInitializer::printCondition(std::string initialConditionName) {
	using namespace GlobalPara::boundaryCondition::_2D;
	std::stringstream ss;
	ss << "InitialCondition = "<< initialConditionName <<"\n";
	LogWriter::logAndPrint(ss.str());
}

void FieldInitializer::setInitialUniform() {
	// 均匀场
	auto& elements = FVM_2D::getInstance()->elements;
	using namespace GlobalPara::boundaryCondition::_2D;
	for (int ie = 0; ie < elements.size(); ie++) {
		Math_2D::ruvp_2_U(inf::ruvp, elements[ie].U, GlobalPara::constant::gamma);
	}
}

void FieldInitializer::setInitialIsentropicVortex() {
	double vortex_x = 0;
	double vortex_y = 0;
	double vortex_strength = 1;
	TomlFileManager::getInstance()->getValueOrExit("initialCondition.isentropicVortex.vortex_x", vortex_x);
	TomlFileManager::getInstance()->getValueOrExit("initialCondition.isentropicVortex.vortex_y", vortex_y);
	TomlFileManager::getInstance()->getValueOrExit("initialCondition.isentropicVortex.vortex_strength", vortex_strength);

	double xc = vortex_x;
	double yc = vortex_y;
	double chi = vortex_strength;
	double chi2 = chi * chi;
	const double* ruvp_inf = GlobalPara::boundaryCondition::_2D::inf::ruvp;
	const double gamma = GlobalPara::constant::gamma;
	const double ga1 = gamma - 1.0;
	constexpr double PI = U2NITS::Math::PI;
	constexpr double two_pi = 2.0 * PI;
	constexpr double pi2 = PI * PI;
	double c_du = chi / two_pi;// du表达式的系数，正数，与涡强度有关
	double c_dT = ga1 * chi2 / (8. * gamma * pi2);// dT表达式的系数，正数，与涡强度有关

	for (Element_2D& e: FVM_2D::getInstance()->elements) {
		double dx = e.x - xc;
		double dy = e.y - yc;
		double r2 = dx * dx + dy * dy;
		double one_minus_r2 = 1.0 - r2;
		double c_distance = exp(0.5 * one_minus_r2);
		double du = c_du * c_distance * (-dy);
		double dv = c_du * c_distance * dx;
		double u = ruvp_inf[1] + du;
		double v = ruvp_inf[2] + dv;
		double dT = -c_dT * exp(one_minus_r2);
		double rho = pow(ruvp_inf[3] + dT, 1. / ga1);
		double p = rho * (ruvp_inf[3] + dT);
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
		if (CBoundaryDoubleShockReflect::getInstance()->isUpStreamOfShock_forElement(element.x, element.y)) {
			Math_2D::ruvp_2_U(inlet::ruvp, element.U, GlobalPara::constant::gamma);
		}
		else {
			Math_2D::ruvp_2_U(outlet::ruvp, element.U, GlobalPara::constant::gamma);
		}
	}
}
