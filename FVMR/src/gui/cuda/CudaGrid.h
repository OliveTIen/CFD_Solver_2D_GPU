#pragma once 
#include <glm/glm.hpp>
#include <vector>
#include "../object/ShaderBasic.h"
#include "../object/Object.h"
#include "../io/typeDefine.h"
#include "../../gpu/datatype/FieldSoA.h"
#include "StructColorMap.h"

class CudaGrid : public Object {
public:
	const char* vert_line_path = nullptr;
	const char* frag_line_path = nullptr;
	const char* vert_tri_path = nullptr;// triangle
	const char* frag_tri_path = nullptr;

	std::vector<glm::vec3> vertices;
	std::vector<glm::vec3> vertices_copy;// used for reset
	std::vector<glm::uvec2> indices_line;

	std::vector<float> vertices_triangle;// coordinate and color
	std::vector<int> indices_triangle;// 

	std::vector<std::pair<std::string, std::vector<glm::uvec2>>> boundaries;
	
	std::vector<glm::vec3> velocities;
	std::vector<glm::vec3> X_hat;
	std::vector<glm::vec3> G;// point force
	std::vector<int> E;// point indices for each edge
	std::vector<float> L;// original edge length
	glm::vec3 gravity_vector = glm::vec3(0.0f, -9.8f, 0.0f);
	unsigned int slices = 20;
	unsigned int vao_line = -1;
	unsigned int vbo_line = -1;
	unsigned int ibo_line = -1;
	unsigned int vao_triangle = -1;
	unsigned int vbo_triangle = -1;
	unsigned int ibo_triangle = -1;
	glm::vec4 userColor = glm::vec4(0.0f, 0.2f, 0.0f, 0.7f);
	glm::mat4 matrix_model = glm::mat4(1.0f);
	float spring_k = 8.0e3f;
	float mass = 1;
	float velocity_damping = 0.99;
	int substep_num = 32;
	ShaderBasic* shader_line = nullptr;
	ShaderBasic* shader_triangle = nullptr;
	bool enable_render_line = true;
	bool enable_render_blend = true;
	float render_blend_alpha = 0.8;
	bool polygonmode_line = false;
	bool render_using_device_data = true;
	struct cudaGraphicsResource* cuda_vbo_resource = nullptr;  // handles OpenGL-CUDA exchange
	GPU::ColorMap colormap_host;
	GPU::ColorMap colormap_device;
	GPU::DataSoA nodeField_device;// only used in CudaGrid

	CudaGrid(const char* _vert_line_path, const char* _frag_line_path, const char* _vert_tri_path, const char* _frag_tri_path) {
		vert_line_path = _vert_line_path;
		frag_line_path = _frag_line_path;
		vert_tri_path = _vert_tri_path;
		frag_tri_path = _frag_tri_path;
		//init(); // delay init
	}
	virtual ~CudaGrid() override {
		cleanup();
	}
	virtual void init() override;
	virtual void reset() override;
	virtual void update(float dt) override;
	virtual void render(const glm::mat4& projection, const glm::mat4& view) const override;
	virtual void cleanup() override;

private:
	bool m_initialized = false;

	void reset_position();

	void initRenderingData_UnstructuredGrid();
	void initPhysicsPointData_UnstructuredGrid();
	void initPhysicsEdgeData_UnstructuredGrid();

	void Get_Gradient(float dt);

	void Set_Boundary_Condition_Fixed(size_t index);

	static void setTriangleVertexPosition(std::vector<float>& tri, const std::vector<glm::vec3>& coord);
	void setTriangleVertexColor();
	void setTriangleVertexColor_GPU();
	glm::vec3 getColor(const GPU::ColorMap& colormap, float input);
	glm::vec3 float3_to_vec3(const float3& input) {
		return glm::vec3(input.x, input.y, input.z);
	}
	glm::vec3 float4_to_vec3(const float4& input) {// discard input.w
		return glm::vec3(input.x, input.y, input.z);
	}
	float3 vec3_to_float3(const glm::vec3& input) {
		return float3{ input.x, input.y, input.z };
	}
};