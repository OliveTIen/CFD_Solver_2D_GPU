#pragma once
#include <glad/glad.h>
#include <cuda_runtime_api.h>
#include "CudaGrid.h"
namespace GPU {
	void setVerticesColor(GLuint vbo, cudaGraphicsResource* cuda_vbo_resource, CudaGrid* cudaGrid);
}