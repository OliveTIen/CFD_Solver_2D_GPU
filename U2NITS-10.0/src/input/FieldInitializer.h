#ifndef FIELD_INITIALIZER_H
#define FIELD_INITIALIZER_H

#include <string>
class FieldInitializer {
public:
	static void setInitialAndBoundaryCondition();


private:
	static void printCondition(std::string initialConditionName);
	static void setInitialUniform();
	static void setInitialIsentropicVortex(double xc, double yc, double chi, const double* ruvp0);
	static void setInitialShockTube();
	static void setInitialDoubleShockReflection();
};

#endif // !INITIALIZER_H
