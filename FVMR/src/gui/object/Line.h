#pragma once
#include <glm/glm.hpp>
#include "ShaderLine.h"
#include "Object.h"
class Line :public  Object {
public:
	bool use_user_color = false;
	glm::vec4 userColor = glm::vec4(1.0f, 0.0f, 0.0f, 1.0f);
	glm::mat4 matrix_model = glm::mat4(1.0f);
	
	static float vertices[];
	static unsigned int vertex_number;
	static unsigned int vao;
	static unsigned int vbo;
	static ShaderLine* shader;

	Line() { init(); }
	virtual ~Line() override { cleanup(); }
	virtual void init() override;
	virtual void reset() override {};
	virtual void update(float dt) override {};
	virtual void render(const glm::mat4& projection, const glm::mat4& view) const override;
	virtual void cleanup() override;
	static void init_static_buffer();
	static void delete_static_buffer();

private:
	static int instance_count;
	static bool static_buffer_allocated;
};