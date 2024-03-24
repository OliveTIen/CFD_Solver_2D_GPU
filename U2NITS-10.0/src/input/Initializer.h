#ifndef INITIALIZER_H
#define INITIALIZER_H

#include <string>
class Initializer {
public:
	static void setInitialAndBoundaryCondition();

private:
	static void printCondition(std::string initialConditionName);
	static void setInitialUniform();
	static void setInitialIsentropicVortex(double xc, double yc, double chi, const double* ruvp0);
	static void setInitialShockTube();
};

#endif // !INITIALIZER_H
