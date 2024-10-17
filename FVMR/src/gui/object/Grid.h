#pragma once
#include <glm/glm.hpp>
#include <vector>
#include "ShaderLine.h"
#include "Object.h"

class Grid : public Object {
public:
	bool b_structured = false;
	const char* file_path;
	const char* vert_path;
	const char* frag_path;
	std::vector<glm::vec3> vertices;// used for opengl rendering
	std::vector<glm::vec3> vertices_copy;// used for reset
	std::vector<glm::uvec2> indices;// used for opengl rendering
	std::vector<std::pair<std::string, std::vector<glm::uvec2>>> boundaries;
	unsigned int index_length;// used for opengl rendering
	std::vector<glm::vec3> velocities;
	std::vector<glm::vec3> X_hat;
	std::vector<glm::vec3> G;// point force
	std::vector<int> E;// point indices for each edge
	std::vector<float> L;// original edge length
	glm::vec3 gravity_vector = glm::vec3(0.0f, -9.8f, 0.0f);
	unsigned int slices = 20;
	unsigned int vao;
	unsigned int vbo;
	unsigned int ibo;
	glm::vec4 userColor = glm::vec4(0.0f, 0.0f, 0.0f, 1.0f);
	glm::mat4 matrix_model = glm::mat4(1.0f);
	float spring_k = 8.0e3f;
	float mass = 1;
	float velocity_damping = 0.99;
	int substep_num = 32;
	ShaderLine* shader = nullptr;
	static int instance_count;

	Grid(const char* _file_path, const char* _vert_path, const char* _frag_path) { 
		file_path = _file_path;
		vert_path = _vert_path;
		frag_path = _frag_path;
		init(); 
	}
	virtual ~Grid() override { cleanup(); }

	virtual void init() override;
	virtual void reset() override { reset_position(); }
	virtual void update(float dt) override;
	virtual void render(const glm::mat4& projection, const glm::mat4& view) const override;
	virtual void cleanup() override;

private:
	bool buffer_allocated = false;

	void reset_position() {
		if (b_structured) {
			int n = slices + 1;
			float dL = 1.0f / slices;
			for (int j = 0; j < n; ++j) {
				for (int i = 0; i < n; ++i) {
					float x = i * dL;
					float y = 0;
					float z = j * dL;
					//vertices.push_back(glm::vec3(x, y, z));
					vertices[j * n + i] = glm::vec3(x, y, z);
				}
			}
		}
		else {
			vertices = vertices_copy;
		}
	}

	// functions
	void initRenderingData_StructuredGrid();
	void initPhysicsPointData_StructuredGrid();
	void initPhysicsEdgeData_StructuredGrid();

	void initRenderingData_UnstructuredGrid(std::string filepath);
	void initPhysicsPointData_UnstructuredGrid();
	void initPhysicsEdgeData_UnstructuredGrid();

	void Get_Gradient(float dt);
	//void Set_Boundary_Condition_Fixed_Structured();

	void Swap(int& a, int& b) {
		int temp = a;
		a = b;
		b = temp;
	}

	void Quick_Sort(int* a, int l, int r) {
		int j;
		if (l < r) {
			j = Quick_Sort_Partition(a, l, r);
			Quick_Sort(a, l, j - 1);
			Quick_Sort(a, j + 1, r);
		}
	};

	int Quick_Sort_Partition(int* a, int l, int r) {
		int pivot_0, pivot_1, i, j;
		pivot_0 = a[l * 2 + 0];
		pivot_1 = a[l * 2 + 1];
		i = l;
		j = r + 1;
		while (true) {
			do ++i; while (i <= r && (a[i * 2] < pivot_0 || a[i * 2] == pivot_0 && a[i * 2 + 1] <= pivot_1));
			do --j; while (a[j * 2] > pivot_0 || a[j * 2] == pivot_0 && a[j * 2 + 1] > pivot_1);
			if (i >= j)	break;
			Swap(a[i * 2], a[j * 2]);
			Swap(a[i * 2 + 1], a[j * 2 + 1]);
		}
		Swap(a[l * 2 + 0], a[j * 2 + 0]);
		Swap(a[l * 2 + 1], a[j * 2 + 1]);
		return j;
	}

	void Set_Boundary_Condition_Fixed_Structured(size_t i, size_t j) {
		int node_num = slices + 1;
		vertices[j * node_num + i] = X_hat[j * node_num + i];
		velocities[j * node_num + i] *= 0;
	}

	void Set_Boundary_Condition_Fixed(size_t index) {
		vertices[index] = X_hat[index];
		velocities[index] *= 0;
	}
};