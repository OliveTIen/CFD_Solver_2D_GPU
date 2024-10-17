#include <glad/glad.h> 
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include "Line.h"

float Line::vertices[14]{
	// x y z R G B A
	0.0f,0.0f,0.0f,0.0f,0.0f,1.0f,1.0f,
	1.0f,0.0f,0.0f,1.0f,0.0f,0.0f,1.0f
};

unsigned int Line::vertex_number = 14;
unsigned int Line::vao = 0;
unsigned int Line::vbo = 0;
ShaderLine* Line::shader = nullptr;
bool Line::static_buffer_allocated = false;
int Line::instance_count = 0;

void Line::init_static_buffer() {
    if (static_buffer_allocated) {
        printf("Warning: static buffer already allocated, cannot initialize.\n");
        return;
    }
    shader = new ShaderLine("./shader/line.vert", "./shader/line.frag");

    glGenVertexArrays(1, &vao);
    glGenBuffers(1, &vbo);

    glBindVertexArray(vao);
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferData(GL_ARRAY_BUFFER, sizeof(float) * Line::vertex_number, Line::vertices, GL_STATIC_DRAW);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 7 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(1, 4, GL_FLOAT, GL_FALSE, 7 * sizeof(float), (void*)(3 * sizeof(float)));
    glEnableVertexAttribArray(1);

    static_buffer_allocated = true;
}

void Line::delete_static_buffer() {
    if (!static_buffer_allocated) {
        printf("Warning: static buffer not initialized, cannot delete.\n");
        return;
    }
    glDeleteVertexArrays(1, &vao);
    glDeleteBuffers(1, &vbo);
    delete shader;
    static_buffer_allocated = false;
}

void Line::init() {
    if (instance_count == 0)init_static_buffer();
    instance_count++;
}

void Line::render(const glm::mat4& projection, const glm::mat4& view) const {
    glBindVertexArray(vao);
    shader->use();
    shader->setMat4("projection", projection);
    shader->setMat4("view", view);
    shader->setMat4("model", matrix_model);
    shader->setBool("use_user_color", use_user_color);
    shader->setVec4("userColor", userColor);
    glDrawArrays(GL_LINES, 0, 2);
    glBindVertexArray(0);
}

void Line::cleanup() {
    instance_count--;
    if (instance_count == 0)delete_static_buffer();
}
