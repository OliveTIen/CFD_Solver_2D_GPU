#pragma once
#include <glm/glm.hpp>
#include "ShaderCube.h"
#include "Object.h"

class Cube : public Object {
public:
	enum BlockType {
		_type_air,
		_type_stone,
		_type_grass,
		_type_dirt,
		_type_diamond_ore,
		_type_obsidian
	};


	glm::vec3 position = glm::vec3(0.0f, 0.0f, 0.0f);
	glm::vec3 scale = glm::vec3(1.0f, 1.0f, 1.0f);
	glm::vec3 rotate_axis = glm::vec3(1.0f, 0.3f, 0.5f);
	float rotate_angle = 0.0f;
	float rotate_speed = 0.0f;
	float mix_weight = 0.0f;
	BlockType block_type = _type_stone;

	static float vertices[];
	static unsigned int vertex_number;
	static unsigned int VAO;
	static unsigned int VBO;
	static ShaderCube* shader;
	static unsigned int texture_stone;
	static unsigned int texture_dirt;
	static unsigned int texture_diamond_ore;
	static unsigned int texture_obsidian;
	static int instance_count;
	static bool static_buffer_initialized;

	Cube() { init(); }
	Cube(float x, float y, float z, BlockType _t = BlockType::_type_stone) {
		position = glm::vec3(x, y, z); 
		block_type = _t;
		init();
	}
	virtual ~Cube()override { cleanup(); }
	virtual void init() override;
	virtual void reset() override {};
	virtual void update(float dt) override { rotate_angle += rotate_speed * dt; };
	virtual void render(const glm::mat4& projection, const glm::mat4& view) const override;
	virtual void cleanup() override;

	static void init_static_buffer();
	static void delete_static_buffer();
};