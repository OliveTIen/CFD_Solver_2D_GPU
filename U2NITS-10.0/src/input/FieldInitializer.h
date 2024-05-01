#ifndef FIELD_INITIALIZER_H
#define FIELD_INITIALIZER_H

#include <string>
class FieldInitializer {
private:
	static FieldInitializer* p_instance;
	int m_initial_type = 0;

public:
	static FieldInitializer* getInstance();
	void setInitialAndBoundaryCondition();
	// ���ó�ʼ����ʽ
	void set_initial_type(int _type) { m_initial_type = _type; }
	// ��TomlFileManager�ṩ��API���г�ʼ����ֻ��tomlFileManager����
	void initialize_using_config(void* tomlFileManager);
	int get_initial_type() { return m_initial_type; }

private:
	FieldInitializer() {};
	void printCondition(std::string initialConditionName);
	void setInitialUniform();
	void setInitialIsentropicVortex();
	void setInitialShockTube();
	void setInitialDoubleShockReflection();
};

#endif // !INITIALIZER_H
