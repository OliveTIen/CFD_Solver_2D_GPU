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
	void set_enabled(bool b) { m_enabled = b; }
	bool get_enabled() { return m_enabled; }
	bool& get_enabled_ref() { return m_enabled; }// 返回引用，可从外部修改私有成员

private:
	GLFWwindow* m_window = nullptr;
	ImGuiIO* m_io = nullptr;
	Control* m_control = nullptr;
	Scene* m_scene = nullptr;
	bool m_enabled = true;

	void panel_ParticleSystem();
	void panel_FVMR();
	void update_imgui();

	void button_ColorPicker(ImVec4& color);
	void buttons_ColorTest();// 用于测试color相关控件
};