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
	// 设置初始化方式
	void set_initial_type(int _type) { m_initial_type = _type; }
	// 用TomlFileManager提供的API进行初始化。只被tomlFileManager调用
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
