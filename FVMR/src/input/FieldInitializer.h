#ifndef FIELD_INITIALIZER_H
#define FIELD_INITIALIZER_H

#include <string>


enum InitialConditionType {
	type_none,
	type_uniform_flow,          // 均匀场
	type_isentropic_vortex,     // 等熵涡
	type_shock_tube,            // 激波管
	type_double_mach_reflection,// 双马赫反射
	type_tecplot_file           // tecplot文件(点格式)
};

class FieldInitializer {
private:
	static FieldInitializer* p_instance;
	int m_initial_type = 1;// uniform flow
	bool b_has_read_config = false;

	// type_shock_tube所用到的局部变量，用于初始化激波位置
	double m_shock_x = 0.0;// 激波横坐标
	double m_shock_y = 0.0;
	double m_shock_normal_x = 1.0;// 激波法向量x分量，法向量从inlet指向outlet
	double m_shock_normal_y = 0.0;

public:
	static FieldInitializer* getInstance();
	void setInitialAndBoundaryCondition();
	// 设置初始化方式
	void set_initial_type(int _type) { m_initial_type = _type; }
	// 用TomlFileManager提供的API进行初始化。只被tomlFileManager调用
	void initialize_using_config();
	int get_initial_type();
	void read_initial_type_from_config();

private:
	FieldInitializer() {};
	void printCondition(std::string initialConditionName);
	void setInitialUniform();
	void setInitialIsentropicVortex();
	void setInitialShockTube();
	void setInitialDoubleShockReflection();
	void setInitialTecplotFile();
};

#endif // !INITIALIZER_H
