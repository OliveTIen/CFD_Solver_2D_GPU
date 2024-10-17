#pragma once
// imgui
#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"
//
#include "Control.h"
#include "Scene.h"

class WinGui {
public:
	void init();
	void render();
	void cleanup();
	bool windowShouldClose();

private:
	GLFWwindow* m_window = nullptr;
	ImGuiIO* m_io = nullptr;
	Control* m_control = nullptr;
	Scene* m_scene = nullptr;

	void panel_ParticleSystem();
	void panel_FVMR();
	void update_imgui();

	void button_ColorPicker(ImVec4& color);
	void buttons_ColorTest();// 用于测试color相关控件
};