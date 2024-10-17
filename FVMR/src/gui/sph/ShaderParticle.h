#pragma once
#include "../object/ShaderBasic.h"
class ShaderParticle : public ShaderBasic {
public:
	ShaderParticle(const char* vertexPath, const char* fragmentPath) :
		ShaderBasic(vertexPath, fragmentPath) {
	}
};
