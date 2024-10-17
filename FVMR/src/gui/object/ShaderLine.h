#ifndef SHADER_LINE_H
#define SHADER_LINE_H
#include "ShaderBasic.h"
class ShaderLine : public ShaderBasic{
public:
	ShaderLine(const char* vertexPath, const char* fragmentPath) :
		ShaderBasic(vertexPath, fragmentPath) {	}
};
#endif