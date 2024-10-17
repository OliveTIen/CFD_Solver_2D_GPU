#pragma once
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
class Camera {
public:
	glm::vec3 position = glm::vec3(0.34f, 1.2f, 2.0f);
	glm::vec3 front = glm::vec3(-0.02f, -0.0628f, -0.0778f);
	glm::vec3 up = glm::vec3(0.0f, 1.0f, 0.0f);
	float fov = 45.0f;
	float aspect_ratio = 1.0f;
	float move_speed = 2.5f;
	float rotate_sensitivity = 0.1f;
	float yaw = -90.0f;	// 0: look at +x dir. 90: +z dir. 180: -x dir
	float pitch = -40.0f;// range: [-89,89]. -89 means look down


	Camera(float _aspect_ratio) { aspect_ratio = _aspect_ratio; setFront(yaw, pitch); }
	Camera(){ setFront(yaw, pitch); }
	glm::mat4 getProjection() {
		return glm::perspective(glm::radians(fov), aspect_ratio, 0.1f, 100.0f);
	}
	glm::mat4 getView() {
		return glm::lookAt(position, position + front, up);
	}
	void setFront(float yaw, float pitch);
	void setRotation(float xoffset, float yoffset);
	
};
