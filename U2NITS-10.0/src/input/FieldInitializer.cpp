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
#include "MaterialManager.h"

FieldInitializer* FieldInitializer::p_instance = nullptr;

FieldInitializer* FieldInitializer::getInstance() {
	if (p_instance == nullptr) {
		p_instance = new FieldInitializer();
	}
	return p_instance;
}

void FieldInitializer::setInitialAndBoundaryCondition() {

	// ���ó�ʼ�ͱ߽�����
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

	case type_tecplot_file:
	{
		printCondition("Tecplot File");
		setInitialTecplotFile();
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

void FieldInitializer::read_initial_type_from_config() {
	TomlFileManager* t = TomlFileManager::getInstance();
	t->getValue("initialCondition.type", m_initial_type);
}
void FieldInitializer::initialize_using_config() {
	/*
	��TomlFileManager�ṩ��API��ʼ��m_initial_type
	����˫��շ��䣬��Ҫ���ȡ����λ�á��ǶȲ���
	*/
	TomlFileManager* t = TomlFileManager::getInstance();
	read_initial_type_from_config();

	/*
	[to do]����Ӧ�ƶ���CBoundaryDoubleShockReflect::initialize_using_config
	*/
	if (m_initial_type == type_double_mach_reflection) {
		myfloat shock_x = 0;
		myfloat shock_y = 0;
		myfloat shock_angle_degree = 60;
		t->getValueOrExit("initialCondition.doubleShockReflection.shock_x", shock_x);
		t->getValueOrExit("initialCondition.doubleShockReflection.shock_y", shock_y);
		t->getValueOrExit("initialCondition.doubleShockReflection.shock_angle_degree", shock_angle_degree);

		if (shock_angle_degree > 90 || shock_angle_degree < 0) {
			std::cout << "shock_angle_degree out of range [0,90]. please try again.\n";
			CExit::pressAnyKeyToExit();
		}

		CBoundaryDoubleShockReflect::getInstance()->set_shock_x_y_angle(shock_x, shock_y, shock_angle_degree);
	}

	b_has_read_config = true;
}

int FieldInitializer::get_initial_type() {
	if (!b_has_read_config) {
		LogWriter::logAndPrintError("!b_has_read_config\n");
		exit(-1);
	}
	return m_initial_type; 
}


void FieldInitializer::printCondition(std::string initialConditionName) {
	using namespace GlobalPara::boundaryCondition::_2D;
	std::stringstream ss;
	ss << "InitialCondition = "<< initialConditionName <<"\n";
	LogWriter::logAndPrint(ss.str());
}

void FieldInitializer::setInitialUniform() {
	// ���ȳ�
	auto& elements = FVM_2D::getInstance()->elements;
	using namespace GlobalPara::boundaryCondition::_2D;
	for (int ie = 0; ie < elements.size(); ie++) {
		Math_2D::ruvp_2_U(inf::ruvp, elements[ie].U, GlobalPara::constant::gamma);
	}
}

void FieldInitializer::setInitialIsentropicVortex() {
	myfloat vortex_x = 0;
	myfloat vortex_y = 0;
	myfloat vortex_strength = 1;
	TomlFileManager::getInstance()->getValueOrExit("initialCondition.isentropicVortex.vortex_x", vortex_x);
	TomlFileManager::getInstance()->getValueOrExit("initialCondition.isentropicVortex.vortex_y", vortex_y);
	TomlFileManager::getInstance()->getValueOrExit("initialCondition.isentropicVortex.vortex_strength", vortex_strength);

	myfloat xc = vortex_x;
	myfloat yc = vortex_y;
	myfloat chi = vortex_strength;
	myfloat chi2 = chi * chi;
	const myfloat* ruvp_inf = GlobalPara::boundaryCondition::_2D::inf::ruvp;
	const myfloat gamma = GlobalPara::constant::gamma;
	const myfloat ga1 = gamma - 1.0;
	constexpr myfloat PI = U2NITS::Math::PI;
	constexpr myfloat two_pi = 2.0 * PI;
	constexpr myfloat pi2 = PI * PI;
	myfloat c_du = chi / two_pi;// du���ʽ��ϵ��������������ǿ���й�
	myfloat c_dT = ga1 * chi2 / (8. * gamma * pi2);// dT���ʽ��ϵ��������������ǿ���й�

	for (Element_2D& e: FVM_2D::getInstance()->elements) {
		myfloat dx = e.x - xc;
		myfloat dy = e.y - yc;
		myfloat r2 = dx * dx + dy * dy;
		myfloat one_minus_r2 = 1.0 - r2;
		myfloat c_distance = exp(0.5 * one_minus_r2);
		myfloat du = c_du * c_distance * (-dy);
		myfloat dv = c_du * c_distance * dx;
		myfloat u = ruvp_inf[1] + du;
		myfloat v = ruvp_inf[2] + dv;
		myfloat dT = -c_dT * exp(one_minus_r2);
		myfloat rho = pow(ruvp_inf[3] + dT, 1. / ga1);
		myfloat p = rho * (ruvp_inf[3] + dT);
		myfloat ruvp[4]{ rho,u,v,p };
		Math_2D::ruvp_2_U(ruvp, e.U, gamma);
	}

}

void FieldInitializer::setInitialShockTube() {
	// �����ܡ���������shock_tube��https://zhuanlan.zhihu.com/p/154508317������D:\tgl\Local\CFD\shock_tube_code
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
	// ˫��շ����ʼֵ���� https://zhuanlan.zhihu.com/p/630069961
	// ����ֱ���ɵ�бʽ(x,y,angle)ȷ��


	constexpr bool isDebug = false;
	if (isDebug) {

		// ��¼�߽����
		auto writeBoundaryCondition = [](myfloat* inlet_ruvp, myfloat* outlet_ruvp, myfloat* inf_ruvp, const int num_ruvp)->void {
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

void FieldInitializer::setInitialTecplotFile() {
	LogWriter::logAndPrintError("unimplemented.\n");
	exit(-1);
	/*
	2024-05-27����
	�������TecplotFileReader
	
	*/
}
