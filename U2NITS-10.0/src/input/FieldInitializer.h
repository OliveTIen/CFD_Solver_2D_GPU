#ifndef FIELD_INITIALIZER_H
#define FIELD_INITIALIZER_H

#include <string>


enum InitialConditionType {
	type_none,
	type_uniform_flow,          // ���ȳ�
	type_isentropic_vortex,     // ������
	type_shock_tube,            // ������
	type_double_mach_reflection,// ˫��շ���
	type_tecplot_file           // tecplot�ļ�(���ʽ)
};

class FieldInitializer {
private:
	static FieldInitializer* p_instance;
	int m_initial_type = 1;// uniform flow
	bool b_has_read_config = false;

public:
	static FieldInitializer* getInstance();
	void setInitialAndBoundaryCondition();
	// ���ó�ʼ����ʽ
	void set_initial_type(int _type) { m_initial_type = _type; }
	// ��TomlFileManager�ṩ��API���г�ʼ����ֻ��tomlFileManager����
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
