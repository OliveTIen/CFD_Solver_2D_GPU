#pragma once
#include "Camera.h"
#include <glad/glad.h>
#include <GLFW/glfw3.h>

class Control {
public:
	enum MOVE_MODE {
		moveMode_normal,
		moveMode_swimming
	};

	static unsigned int SCR_WIDTH;
	static unsigned int SCR_HEIGHT;
	static Camera camera;
	// camera move
	static glm::vec2 lastMouse;
	static int moveMode;
	static bool update_fovy_aspect_on_resize;
	static bool fixed_mouse;
	static bool on_fixedMouseMode_entered;
	// timing
	static float deltaTime;	// time between current frame and last frame
	static float lastTime;
	// render
	static bool enable_depthtest;
	static const char* glsl_version;
	static glm::vec4 clear_color;
	// imgui state
	static bool show_demo_window;
	static bool show_another_window;
	// key state
	static bool key_R_pressed;

	Control();
	~Control();
	static void glfw_error_callback(int error, const char* description) {
		fprintf(stderr, "GLFW Error %d: %s\n", error, description);
	}
	static void framebuffer_resize_callback(GLFWwindow* window, int width, int height);
	static void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods);
	static void mouse_callback(GLFWwindow* window, double xposIn, double yposIn);
	static void scroll_callback(GLFWwindow* window, double xoffset, double yoffset);
	static void mouse_button_callback(GLFWwindow* window, int button, int action, int mods);
	static void processInput(GLFWwindow* window);
	static void setCurMousePos(float _xpos, float _ypos);
	static void clearScreen();
	//void pick();
	//void build();
	//int DetectRayPolygon(glm::vec3 startPoint, glm::vec3 inAABBPoint, glm::vec3 testDir, glm::vec3 border);
	static glm::vec2 getKeyboardInputVector2D(GLFWwindow* window);


private:



};

