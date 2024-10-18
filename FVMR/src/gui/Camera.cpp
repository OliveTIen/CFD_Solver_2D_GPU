#include "Camera.h"
#include <iostream>

void Camera::setFront(float _yaw, float _pitch) {
	front.x = cos(glm::radians(_yaw)) * cos(glm::radians(_pitch));
	front.y = sin(glm::radians(_pitch));
	front.z = sin(glm::radians(_yaw)) * cos(glm::radians(_pitch));
	front = glm::normalize(front);
}

void Camera::setRotation(float xoffset, float yoffset) {
	
	xoffset *= rotate_sensitivity;
	yoffset *= rotate_sensitivity;
	yaw += xoffset;
	pitch += yoffset;

	// limit pitch to avoid screen getting flipped
	if (pitch > 89.0f)
		pitch = 89.0f;
	if (pitch < -89.0f)
		pitch = -89.0f;
	if (yaw > 360.0f)
		yaw -= 360.0f;
	if (yaw < -360.f)
		yaw += 360.0f;

	setFront(yaw, pitch);
}
