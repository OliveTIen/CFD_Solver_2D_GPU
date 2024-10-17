#include "Control.h"
#include <iostream>

unsigned int Control::SCR_WIDTH = 1200;
unsigned int Control::SCR_HEIGHT = 1000;
Camera Control::camera = Camera((float)SCR_WIDTH / (float)SCR_HEIGHT);

glm::vec2 Control::lastMouse = glm::vec2(SCR_WIDTH / 2.0f, SCR_HEIGHT / 2.0f);
int Control::moveMode = moveMode_normal;

bool Control::update_fovy_aspect_on_resize = true;
bool Control::fixed_mouse = false;
bool Control::on_fixedMouseMode_entered = true;
// timing
float Control::deltaTime = 0.0f;	// time between current frame and last frame
float Control::lastTime = 0.0f;

glm::vec4 Control::clear_color = glm::vec4(0.45f, 0.55f, 0.60f, 1.00f);
bool Control::enable_depthtest = true;
const char* Control::glsl_version = "#version 130";
// imgui state
bool Control::show_demo_window = false;
bool Control::show_another_window = false;

bool Control::key_R_pressed = false;

Control::Control() {
}


Control::~Control() {

}

void Control::framebuffer_resize_callback(GLFWwindow* window, int width, int height) {
	glViewport(0, 0, width, height);
	if (update_fovy_aspect_on_resize) {
		camera.aspect_ratio = (float)width / (float)height;
	}

}

void Control::key_callback(GLFWwindow* window, int key, int scancode, int action, int mods) {
	// to avoid to be called every frame
	if (glfwGetKey(window, GLFW_KEY_Q) == GLFW_PRESS) {
		fixed_mouse = !fixed_mouse;
		if (fixed_mouse) {
			glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
			on_fixedMouseMode_entered = true;
		}
		else {
			glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_NORMAL);
		}
	}
}

void Control::mouse_callback(GLFWwindow* window, double xposIn, double yposIn) {
	if (!fixed_mouse)return;



	float xpos = static_cast<float>(xposIn);
	float ypos = static_cast<float>(yposIn);

	static bool isFirstCalled = true;
	if (isFirstCalled) {
		lastMouse.x = xpos;
		lastMouse.y = ypos;
		isFirstCalled = false;
	}
	if (on_fixedMouseMode_entered) {
		lastMouse.x = xpos;
		lastMouse.y = ypos;
		on_fixedMouseMode_entered = false;
	}

	float xoffset = xpos - lastMouse.x;
	float yoffset = lastMouse.y - ypos; // reversed since y-coordinates go from bottom to top
	lastMouse.x = xpos;
	lastMouse.y = ypos;

	camera.setRotation(xoffset, yoffset);
	return;

	float sensitivity = 0.1f; // change this value to your liking
	xoffset *= sensitivity;
	yoffset *= sensitivity;

	// "static" means only initialized once
	static float yaw = -90.0f;	// yaw is initialized to -90.0 degrees since a yaw of 0.0 results in a direction vector pointing to the right so we initially rotate a bit to the left.
	static float pitch = 0.0f;

	yaw += xoffset;
	pitch += yoffset;

	// make sure that when pitch is out of bounds, screen doesn't get flipped
	if (pitch > 89.0f)
		pitch = 89.0f;
	if (pitch < -89.0f)
		pitch = -89.0f;

	camera.setFront(yaw, pitch);
	std::cout << "xoff, yoff: " << xoffset << ", " << yoffset << "\n";
	std::cout << "yaw, pitch: " << yaw << ", " << pitch << "\n";

}

void Control::scroll_callback(GLFWwindow* window, double xoffset, double yoffset) {
	camera.fov -= (float)yoffset;
	if (camera.fov < 1.0f)
		camera.fov = 1.0f;
	if (camera.fov > 45.0f)
		camera.fov = 45.0f;
}

void Control::mouse_button_callback(GLFWwindow* window, int button, int action, int mods) {
	if (action == GLFW_PRESS) switch (button) {
	case GLFW_MOUSE_BUTTON_LEFT:
		std::cout << "Mouse left button clicked!" << std::endl;
		break;
	case GLFW_MOUSE_BUTTON_MIDDLE:
		std::cout << "Mouse middle button clicked!" << std::endl;
		break;
	case GLFW_MOUSE_BUTTON_RIGHT:
		std::cout << "Mouse right button clicked!" << std::endl;
		break;
	default:
		break;
	}
}

void Control::processInput(GLFWwindow* window) {
	//if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
	//	glfwSetWindowShouldClose(window, true);

	// WSAD Move
	float stepLength = static_cast<float>(camera.move_speed * deltaTime);
	glm::vec2 input = getKeyboardInputVector2D(window);
	glm::vec3 camera_right = glm::normalize(glm::cross(camera.front, camera.up));
	glm::vec3 camera_horizontal_front = glm::normalize(glm::cross(glm::vec3(0.0f, 1.0f, 0.0f), camera_right));

	switch (moveMode) {
	case moveMode_normal:
		camera.position += stepLength * camera_horizontal_front * input.y;
		camera.position += stepLength * camera_right * input.x;
		break;
	case moveMode_swimming:
		camera.position += stepLength * camera.front * input.y;
		camera.position += stepLength * camera_right * input.x;
		break;
	}

	// Up Down
	if (glfwGetKey(window, GLFW_KEY_SPACE) == GLFW_PRESS)
		camera.position += stepLength * glm::vec3(0.0f, 1.0f, 0.0f);
	if (glfwGetKey(window, GLFW_KEY_LEFT_CONTROL) == GLFW_PRESS)
		camera.position += stepLength * glm::vec3(0.0f, -1.0f, 0.0f);

	if (glfwGetKey(window, GLFW_KEY_R) == GLFW_PRESS)
		key_R_pressed = true;
}

void Control::setCurMousePos(float _xpos, float _ypos) {


}

glm::vec2 Control::getKeyboardInputVector2D(GLFWwindow* window) {
	int x = 0, y = 0;
	if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS)
		y += 1;
	if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS)
		y -= 1;
	if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS)
		x -= 1;
	if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS)
		x += 1;
	glm::vec2 xy = glm::vec2((float)x, (float)y);
	if (x == 0 && y == 0)return glm::vec2(0.0f, 0.0f);
	else return glm::normalize(xy);
}

void Control::clearScreen() {
	glClearColor(clear_color.x * clear_color.w, clear_color.y * clear_color.w, clear_color.z * clear_color.w, clear_color.w);
	if (enable_depthtest)glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	else glClear(GL_COLOR_BUFFER_BIT);
}
