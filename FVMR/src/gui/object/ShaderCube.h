#ifndef SHADER_CUBE_H
#define SHADER_CUBE_H

#include "ShaderBasic.h"


class ShaderCube : public ShaderBasic
{
public:
	ShaderCube(const char* vertexPath, const char* fragmentPath):
		ShaderBasic(vertexPath,fragmentPath){ }

	void loadTexture(unsigned int &textureIndex, const char* texturePath, GLint internalFormat = GL_RGB, GLenum format = GL_RGB, 
		bool texture_filter_linear = true);
};
#endif