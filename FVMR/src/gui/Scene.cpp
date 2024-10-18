#include <glad/glad.h> // used by glad.c
#include <GLFW/glfw3.h>// used by glfw3.lib
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include "Scene.h"
#include "object/Cube.h"
#include "object/Line.h"
#include "object/Grid.h"
#include "sph/ParticleSystem.h"
#include "cuda/CudaGrid.h"

/*
计划
添加Scene::addObject函数，用于把外部设定好的物体添加进场景。为避免多余复制，需要用智能指针
类似Qt中添加子控件一样，不知道怎么实现。
最好能够动态添加
没必要维护object list，可以直接把object放进scene，无需多层级。

*/

void Scene::init(Control* control) {
	m_control = control;

	objectLists["Cube"] = std::vector<Object*>();
	objectLists["Line"] = std::vector<Object*>();
	//objectLists["Grid"] = std::vector<Object*>();
	objectLists["CudaGrid"] = std::vector<Object*>();
	objectLists["ParticleSystem"] = std::vector<Object*>();

	auto it = objectLists.find("Cube");
	if (it != objectLists.end()) {
		auto& cubes = it->second;
		cubes.push_back(new Cube(2.5f, 1.0f, -3.0f, Cube::_type_diamond_ore));
		cubes.push_back(new Cube(-1.3f, 1.0f, -1.5f, Cube::_type_stone));
		static_cast<Cube*>(cubes[0])->scale = glm::vec3(0.5f, 0.5f, 0.5f);
		static_cast<Cube*>(cubes[0])->rotate_speed = -40.0f;
		static_cast<Cube*>(cubes[1])->rotate_speed = 30.0f;
		const int size = 10;
		for (int i = 0; i < size; i++) {
			for (int j = 0; j < size; j++) {
				float x = i - size / 2;
				float z = j - size / 2;
				float y = -3.0f;
				cubes.push_back(new Cube(x, y, z, Cube::_type_dirt));
			}

		}

	}

	it = objectLists.find("Line");
	if (it != objectLists.end()) {
		// add axes (x-red, y-green, z-blue)
		auto& lines = it->second;
		for (int i = 0; i < 6; i++) {
			lines.push_back(new Line);
		}
		static_cast<Line*>(lines[0])->use_user_color = true;
		static_cast<Line*>(lines[1])->use_user_color = true;
		static_cast<Line*>(lines[2])->use_user_color = true;
		static_cast<Line*>(lines[0])->userColor = glm::vec4(1.0f, 0.0f, 0.0f, 1.0f);
		static_cast<Line*>(lines[1])->userColor = glm::vec4(0.0f, 1.0f, 0.0f, 1.0f);
		static_cast<Line*>(lines[2])->userColor = glm::vec4(0.0f, 0.0f, 1.0f, 1.0f);
		static_cast<Line*>(lines[1])->matrix_model = glm::rotate(static_cast<Line*>(lines[1])->matrix_model, glm::half_pi<float>(), glm::vec3(0.0f, 0.0f, 1.0f));
		static_cast<Line*>(lines[2])->matrix_model = glm::rotate(static_cast<Line*>(lines[2])->matrix_model, glm::half_pi<float>(), glm::vec3(0.0f, -1.0f, 0.0f));
	}

	it = objectLists.find("Grid");
	if (it != objectLists.end()) {
		auto& grids = it->second;
		grids.push_back(new Grid("./resources/shocktube.su2", "./shader/line_black.vert", "./shader/line_black.frag"));
	}

	it = objectLists.find("CudaGrid");
	if (it != objectLists.end()) {
		auto& grids = it->second;
		grids.push_back(new CudaGrid(
			"./shader/line_black.vert", "./shader/line_black.frag", 
			"./shader/cuda_grid.vert", "./shader/cuda_grid.frag"
		));
	}

	it = objectLists.find("ParticleSystem");
	if (it != objectLists.end()) {
		auto& list = it->second;
		list.push_back(new ParticleSystem);
	}
}

void Scene::update() {
	// update object positions and rotations
	const Control& control = *m_control;
	for (auto& list : objectLists) {
		for (auto& object : list.second) {
			if (control.key_R_pressed)object->reset();
			object->update(control.deltaTime);
		}
	}
	// reset control key
	control.key_R_pressed = false;

}

void Scene::render() const {
	const Control& control = *m_control;

	// update transformations
	glm::mat4 projection = control.camera.getProjection();
	glm::mat4 view = control.camera.getView();

	for (auto& list : objectLists) {
		for (auto& object : list.second) {
			object->render(projection, view);
		}
	}
}

void Scene::cleanup() {
	auto it = objectLists.find("Cube");
	if (it != objectLists.end()) {
		Cube::delete_static_buffer();
	}

	it = objectLists.find("Line");
	if (it != objectLists.end()) {
		Line::delete_static_buffer();
	}
}
