// stl
#include <stdio.h>
#include <iostream>
#include <vector>
// opengl
#define GL_SILENCE_DEPRECATION
#if defined(IMGUI_IMPL_OPENGL_ES2)
#include <GLES2/gl2.h>
#endif
#include <glad/glad.h> // used by glad.c
#include <GLFW/glfw3.h>// used by glfw3.lib
// opengl math
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#if defined(_MSC_VER) && (_MSC_VER >= 1900) && !defined(IMGUI_DISABLE_WIN32_FUNCTIONS)
#pragma comment(lib, "legacy_stdio_definitions")
#endif

#include "WinGui.h"
#include "sph/ParticleSystem.h"
#include "cuda/CudaGrid.h"


void WinGui::init() {
    if (!m_enabled)return;
    m_control = new Control();
    m_scene = new Scene();

    const char* window_title = "OpenGL Project";

    glfwSetErrorCallback(Control::glfw_error_callback);

    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

#ifdef __APPLE__
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif

    // Create window with graphics context
    m_window = glfwCreateWindow(m_control->SCR_WIDTH, m_control->SCR_HEIGHT, window_title, nullptr, nullptr);
    if (m_window == nullptr) {
        std::cout << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
        //return -1;
        exit(-1);
    }
    glfwSetWindowTitle(m_window, window_title);
    glfwMakeContextCurrent(m_window);
    glfwSetFramebufferSizeCallback(m_window, Control::framebuffer_resize_callback);
    glfwSetKeyCallback(m_window, Control::key_callback);
    glfwSetCursorPosCallback(m_window, Control::mouse_callback);
    glfwSetScrollCallback(m_window, Control::scroll_callback);
    glfwSwapInterval(1); // double buffer exchange interval time 双缓冲交换间隔时间

    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)) {
        std::cout << "Failed to initialize GLAD" << std::endl;
        //return -1;
        exit(-1);
    }

    if (m_control->enable_depthtest)glEnable(GL_DEPTH_TEST);




    m_scene->init(m_control);

    // setup imgui context
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO();
    m_io = &io;
    io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;     // Enable Keyboard Controls
    io.ConfigFlags |= ImGuiConfigFlags_NavEnableGamepad;      // Enable Gamepad Controls
    ImGui::StyleColorsDark();
    // setup imgui platform/renderer backends
    ImGui_ImplGlfw_InitForOpenGL(m_window, true);
    ImGui_ImplOpenGL3_Init(m_control->glsl_version);
    // setup imgui font "microsoft yahei"
    ImFont* font = io.Fonts->AddFontFromFileTTF("c:\\Windows\\Fonts\\msyh.ttc", 20.0f, nullptr, io.Fonts->GetGlyphRangesChineseFull());
    IM_ASSERT(font != nullptr);// 判断是否成功
}



void WinGui::render() {
    if (!m_enabled)return;
    auto& control = *m_control;
    float current_time = static_cast<float>(glfwGetTime());
    control.deltaTime = current_time - control.lastTime;
    control.lastTime = current_time;
    control.processInput(m_window);
    control.clearScreen();

    m_scene->update();
    m_scene->render();

    update_imgui();
    ImGui::Render();
    ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

    glfwSwapBuffers(m_window);
    glfwPollEvents();

}

void WinGui::cleanup() {
    if (!m_enabled)return;

    m_scene->cleanup();

    // cleanup imgui
    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();

    glfwDestroyWindow(m_window);
    glfwTerminate();

    delete m_control;
    delete m_scene;
}

void WinGui::panel_ParticleSystem() {
    // find scene
    ParticleSystem* ps = nullptr;
    if (m_scene != nullptr) {
        auto it = m_scene->objectLists.find("ParticleSystem");
        if (it != m_scene->objectLists.end()) {
            auto& list = it->second;
            if (!list.empty()) {
                ps = (ParticleSystem*)list[0];
            }
        }
    }
    if (ps == nullptr)return;

    if (ImGui::CollapsingHeader("Particle System")) {
        // gravity

        ImGui::Checkbox("Enable Rendering", &ps->enable_rendering);
        auto& para = ps->physicalParameters;
        ImGui::SliderFloat("gravity", &para.gravity, 0.0f, 20.0f);
        float g_dir[3]{ para.gravity_direction.x,para.gravity_direction.y,para.gravity_direction.z };
        ImGui::SliderFloat3("gravity.dir", g_dir, -1, 1);
        para.gravity_direction.x = g_dir[0];
        para.gravity_direction.y = g_dir[1];
        para.gravity_direction.z = g_dir[2];
        // 
        ImGui::Checkbox("Print on console", &para.printOnConsole);
        ImGui::Checkbox("Rigid Collision", &para.applyRigidCollision);
        ImGui::Checkbox("pressure force", &para.applyPressureForce);
        ImGui::Checkbox("viscous force", &para.applyViscousForce);
        ImGui::SliderFloat("air drag", &para.air_drag, 0.0f, 20.0f);
        int bounceType = para.bounceType;
        ImGui::DragInt("bounce type", &bounceType);
        if (bounceType > 1)bounceType = 1;
        if (bounceType < 0)bounceType = 0;
        para.bounceType = bounceType;
        // 
        //ImGui::InputFloat("pointRadius_render", &ps->pointRadius_render, 0.01f, 1.0f, "%.3f");
        ImGui::SliderFloat("pointRadius_render", &ps->pointRadius_render, 0, 100);
        ImGui::SliderFloat("rigidCollision_radius", &para.rigidCollision_radius, 0, 1);
        ImGui::SliderFloat("rigidCollision_spring", &para.rigidCollision_spring, 0, 1e4);

        ImGui::SliderFloat("bounce0.spring", &para.bounce0_spring, 0, 1e5);
        ImGui::SliderFloat("bounce0.norm_drag", &para.bounce0_norm_drag, 0, 1e2);
        ImGui::SliderFloat("bounce0.wall_friction", &para.bounce0_wall_friction, 0, 1e2);
        ImGui::SliderFloat("bounce1.bounce1_damping_coefficient", &para.bounce1_damping_coefficient, 0, 1);
        ImGui::SliderFloat("bounce1.bounce1_shearSpeedAddition", &para.bounce1_shearSpeedAddition, 0, 1);
        ImGui::Text("The following items take effect after reset");
        ImGui::SliderFloat("fluid.m_standard_distance", &para.m_standard_distance, 1e-3, 1);
        ImGui::SliderFloat("fluid.smooth_length", &para.smooth_length, 1e-3, 1);
        ImGui::Text("fluid.m_rest_density: %f", para.m_rest_density);
        ImGui::Checkbox("fluid.m_clip_smooth_function", &para.m_clip_smooth_function);

    }

}

// Helper to wire demo markers located in code to an interactive browser
typedef void (*ImGuiDemoMarkerCallback)(const char* file, int line, const char* section, void* user_data);
extern ImGuiDemoMarkerCallback      GImGuiDemoMarkerCallback;
extern void* GImGuiDemoMarkerCallbackUserData;
//ImGuiDemoMarkerCallback             GImGuiDemoMarkerCallback = NULL;
//void* GImGuiDemoMarkerCallbackUserData = NULL;
#define IMGUI_DEMO_MARKER(section)  do { if (GImGuiDemoMarkerCallback != NULL) GImGuiDemoMarkerCallback(__FILE__, __LINE__, section, GImGuiDemoMarkerCallbackUserData); } while (0)
#define IM_ARRAYSIZE(_ARR)          ((int)(sizeof(_ARR) / sizeof(*(_ARR))))     // Size of a static C-style array. Don't use on pointers!

static void HelpMarker(const char* desc) {
    ImGui::TextDisabled("(?)");
    if (ImGui::BeginItemTooltip()) {
        ImGui::PushTextWrapPos(ImGui::GetFontSize() * 35.0f);
        ImGui::TextUnformatted(desc);
        ImGui::PopTextWrapPos();
        ImGui::EndTooltip();
    }
}



void WinGui::panel_FVMR() {
    // find scene
    CudaGrid* cudaGrid = nullptr;
    if (m_scene != nullptr) {
        auto it = m_scene->objectLists.find("CudaGrid");
        if (it != m_scene->objectLists.end()) {
            auto& list = it->second;
            if (!list.empty()) {
                cudaGrid = (CudaGrid*)list[0];
            }
        }
    }
    if (cudaGrid == nullptr)return;

    // ui
    if (ImGui::CollapsingHeader("FVMR")) {
        IMGUI_DEMO_MARKER("Widgets/Color");
        const int num_control_point = cudaGrid->colormap_host.num_control_point;

        std::vector<ImVec4> colors(num_control_point);
        std::vector<float> controlPoints(num_control_point);

        for (int i = 0; i < num_control_point; i++) {
            const auto& data = cudaGrid->colormap_host.data;
            colors[i] = { data[i].x,data[i].y,data[i].z,1.0};
            controlPoints[i] = data->w;

        }

        for (int i = 0; i < num_control_point; i++) {
            ImGui::PushID(i);
            button_ColorPicker(colors[i]);
            ImGui::PopID();
        }

        for (int i = 0; i < num_control_point; i++) {
            auto& data = cudaGrid->colormap_host.data;
            data[i] = { colors[i].x,colors[i].y,colors[i].z,controlPoints[i] };
        }

    }


}

void WinGui::update_imgui() {
    Control& control = *m_control;
    bool& show_demo_window = control.show_demo_window;
    bool& show_another_window = control.show_another_window;
    glm::vec4& clear_color = control.clear_color;
    ImGuiIO& io = *m_io;


    // Start the Dear ImGui frame
    ImGui_ImplOpenGL3_NewFrame();
    ImGui_ImplGlfw_NewFrame();
    ImGui::NewFrame();

    // 1. Show the big demo window (Most of the sample code is in ImGui::ShowDemoWindow()! You can browse its code to learn more about Dear ImGui!).
    if (show_demo_window)
        ImGui::ShowDemoWindow(&show_demo_window);

    // 2. Show a simple window that we create ourselves. We use a Begin/End pair to create a named window.
    {

        ImGui::Begin("Control");                          // Create a window called "Hello, world!" and append into it.

        ImGui::Text("Application average %.3f ms/frame (%.1f FPS)", 1000.0f / io.Framerate, io.Framerate);
        ImGui::Text("control.dt %.3f ", control.deltaTime);

        if (ImGui::CollapsingHeader("Help")) {
            ImGui::Text("This is some useful text.");               // Display some text (you can use a format strings too)
            ImGui::Checkbox("Demo Window", &show_demo_window);      // Edit bools storing our window open/close state
            ImGui::Checkbox("Another Window", &show_another_window);
        }

        if (ImGui::CollapsingHeader("Camera")) {
            Camera& camera = control.camera;
            ImGui::ColorEdit3("clear color", (float*)&clear_color); // Edit 3 floats representing a color
            ImGui::Checkbox("update aspect ratio on window resize", &control.update_fovy_aspect_on_resize);
            ImGui::SliderFloat("glm_fovy", &camera.fov, 0.0f, 100.0f);
            ImGui::SliderFloat("camera.aspect_ratio", &camera.aspect_ratio, 0.0f, 5.0f);
            const char* move_mode_str = (control.moveMode == control.moveMode_normal) ? "normal" : "swimming";
            ImGui::Text("Camera move mode: %s", move_mode_str);
            ImGui::SameLine();
            if (ImGui::Button("Change move mode")) {
                if ((control.moveMode == control.moveMode_normal))control.moveMode = control.moveMode_swimming;
                if ((control.moveMode == control.moveMode_swimming))control.moveMode = control.moveMode_normal;
            }
            ImGui::Text("Camera position: %.3f, %.3f, %.3f", camera.position.x, camera.position.y, camera.position.z);
            ImGui::Text("Camera front: %.3f, %.3f, %.3f", camera.front.x, camera.front.y, camera.front.z);
            ImGui::Text("Camera up: %.3f, %.3f, %.3f", camera.up.x, camera.up.y, camera.up.z);
        }
        panel_ParticleSystem();
        panel_FVMR();
        ImGui::End();
    }
}

void WinGui::button_ColorPicker(ImVec4& color) {
    static int button_count = 0;
    button_count++;
    std::string button_id = "button_ColorPicker_" + std::to_string(button_count);

    // Generate a default palette. The palette will persist and can be edited.
    static bool saved_palette_init = true;
    static ImVec4 saved_palette[32] = {};
    if (saved_palette_init) {
        for (int n = 0; n < IM_ARRAYSIZE(saved_palette); n++) {
            ImGui::ColorConvertHSVtoRGB(n / 31.0f, 0.8f, 0.8f,
                saved_palette[n].x, saved_palette[n].y, saved_palette[n].z);
            saved_palette[n].w = 1.0f; // Alpha
        }
        saved_palette_init = false;
    }

    // color
    static ImVec4 backup_color;
    ImGuiColorEditFlags misc_flags = ImGuiColorEditFlags_AlphaPreviewHalf | ImGuiColorEditFlags_AlphaPreview | ImGuiColorEditFlags_NoOptions;
    bool open_popup = ImGui::ColorButton("MyColor##3b", color, misc_flags);
    ImGui::SameLine(0, ImGui::GetStyle().ItemInnerSpacing.x);
    if (open_popup) {
        ImGui::OpenPopup("color_picker");
        backup_color = color;
    }
    if (ImGui::BeginPopup("color_picker")) {
        ImGui::Text("MY CUSTOM COLOR PICKER WITH AN AMAZING PALETTE!");
        ImGui::Separator();
        ImGui::ColorPicker4("##picker", (float*)&color, misc_flags | ImGuiColorEditFlags_NoSidePreview | ImGuiColorEditFlags_NoSmallPreview);
        ImGui::SameLine();

        ImGui::BeginGroup(); // Lock X position
        ImGui::Text("Current");
        ImGui::ColorButton("##current", color, ImGuiColorEditFlags_NoPicker | ImGuiColorEditFlags_AlphaPreviewHalf, ImVec2(60, 40));
        ImGui::Text("Previous");
        if (ImGui::ColorButton("##previous", backup_color, ImGuiColorEditFlags_NoPicker | ImGuiColorEditFlags_AlphaPreviewHalf, ImVec2(60, 40)))
            color = backup_color;// 点击previous按钮会还原上次的颜色
        if (ImGui::IsItemHovered())
            ImGui::SetMouseCursor(7);// 若当前元素被悬浮，cursor改为hand
        ImGui::Separator();
        ImGui::Text("Palette");
        for (int n = 0; n < IM_ARRAYSIZE(saved_palette); n++) {
            ImGui::PushID(n);
            if ((n % 8) != 0)
                ImGui::SameLine(0.0f, ImGui::GetStyle().ItemSpacing.y);

            ImGuiColorEditFlags palette_button_flags = ImGuiColorEditFlags_NoAlpha | ImGuiColorEditFlags_NoPicker | ImGuiColorEditFlags_NoTooltip;
            if (ImGui::ColorButton("##palette", saved_palette[n], palette_button_flags, ImVec2(20, 20)))
                color = ImVec4(saved_palette[n].x, saved_palette[n].y, saved_palette[n].z, color.w); // Preserve alpha!

            // Allow user to drop colors into each palette entry. Note that ColorButton() is already a
            // drag source by default, unless specifying the ImGuiColorEditFlags_NoDragDrop flag.
            if (ImGui::BeginDragDropTarget()) {
                if (const ImGuiPayload* payload = ImGui::AcceptDragDropPayload(IMGUI_PAYLOAD_TYPE_COLOR_3F))
                    memcpy((float*)&saved_palette[n], payload->Data, sizeof(float) * 3);
                if (const ImGuiPayload* payload = ImGui::AcceptDragDropPayload(IMGUI_PAYLOAD_TYPE_COLOR_4F))
                    memcpy((float*)&saved_palette[n], payload->Data, sizeof(float) * 4);
                ImGui::EndDragDropTarget();
            }

            ImGui::PopID();
        }
        ImGui::EndGroup();
        ImGui::EndPopup();
    }
}

void WinGui::buttons_ColorTest() {

    static ImVec4 color = ImVec4(114.0f / 255.0f, 144.0f / 255.0f, 154.0f / 255.0f, 200.0f / 255.0f);

    static bool alpha_preview = true;
    static bool alpha_half_preview = true;
    static bool drag_and_drop = true;
    static bool options_menu = true;
    static bool hdr = false;
    ImGui::SeparatorText("Options");
    ImGui::Checkbox("With Drag and Drop", &drag_and_drop);
    ImGui::Checkbox("With Options Menu", &options_menu); ImGui::SameLine(); HelpMarker("Right-click on the individual color widget to show options.");
    ImGui::Checkbox("With HDR", &hdr); ImGui::SameLine(); HelpMarker("Currently all this does is to lift the 0..1 limits on dragging widgets.");
    ImGuiColorEditFlags misc_flags = (hdr ? ImGuiColorEditFlags_HDR : 0) | (drag_and_drop ? 0 : ImGuiColorEditFlags_NoDragDrop) | (alpha_half_preview ? ImGuiColorEditFlags_AlphaPreviewHalf : (alpha_preview ? ImGuiColorEditFlags_AlphaPreview : 0)) | (options_menu ? 0 : ImGuiColorEditFlags_NoOptions);



    IMGUI_DEMO_MARKER("Widgets/Color/ColorPicker");
    ImGui::SeparatorText("Color picker");
    static bool alpha = true;
    static bool alpha_bar = true;
    static bool side_preview = true;
    static bool ref_color = false;
    static ImVec4 ref_color_v(1.0f, 0.0f, 1.0f, 0.5f);
    static int display_mode = 0;
    static int picker_mode = 0;
    ImGui::Checkbox("With Alpha", &alpha);
    ImGui::Checkbox("With Alpha Bar", &alpha_bar);
    ImGui::Checkbox("With Side Preview", &side_preview);
    if (side_preview) {
        ImGui::SameLine();
        ImGui::Checkbox("With Ref Color", &ref_color);
        if (ref_color) {
            ImGui::SameLine();
            ImGui::ColorEdit4("##RefColor", &ref_color_v.x, ImGuiColorEditFlags_NoInputs | misc_flags);
        }
    }
    ImGui::Combo("Display Mode", &display_mode, "Auto/Current\0None\0RGB Only\0HSV Only\0Hex Only\0");
    ImGui::SameLine(); HelpMarker(
        "ColorEdit defaults to displaying RGB inputs if you don't specify a display mode, "
        "but the user can change it with a right-click on those inputs.\n\nColorPicker defaults to displaying RGB+HSV+Hex "
        "if you don't specify a display mode.\n\nYou can change the defaults using SetColorEditOptions().");
    ImGui::SameLine(); HelpMarker("When not specified explicitly (Auto/Current mode), user can right-click the picker to change mode.");
    ImGuiColorEditFlags flags = misc_flags;
    if (!alpha)            flags |= ImGuiColorEditFlags_NoAlpha;        // This is by default if you call ColorPicker3() instead of ColorPicker4()
    if (alpha_bar)         flags |= ImGuiColorEditFlags_AlphaBar;
    ImGui::ColorPicker4("MyColor##4", (float*)&color, flags, ref_color ? &ref_color_v.x : NULL);


    // Always display a small version of both types of pickers
    // (that's in order to make it more visible in the demo to people who are skimming quickly through it)
    ImGui::Text("Both types:");
    float w = (ImGui::GetContentRegionAvail().x - ImGui::GetStyle().ItemSpacing.y) * 0.40f;
    ImGui::SetNextItemWidth(w);
    ImGui::ColorPicker3("##MyColor##5", (float*)&color, ImGuiColorEditFlags_PickerHueBar | ImGuiColorEditFlags_NoSidePreview | ImGuiColorEditFlags_NoInputs | ImGuiColorEditFlags_NoAlpha);
    ImGui::SameLine();
    ImGui::SetNextItemWidth(w);
    ImGui::ColorPicker3("##MyColor##6", (float*)&color, ImGuiColorEditFlags_PickerHueWheel | ImGuiColorEditFlags_NoSidePreview | ImGuiColorEditFlags_NoInputs | ImGuiColorEditFlags_NoAlpha);
}



bool WinGui::windowShouldClose() {
    if (!m_enabled)return true;
    return glfwWindowShouldClose(m_window);
}
