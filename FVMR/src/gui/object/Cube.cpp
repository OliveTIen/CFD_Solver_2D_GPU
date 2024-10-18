
#include <glad/glad.h> // 和glad.c对应
#include <GLFW/glfw3.h>// 和glfw3.lib对应

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include "Cube.h"

float Cube::vertices[] = {
    // each vertex has 5 floats: 3 of aPos and 2 of aTexCoord
    // layout (location = 0) in vec3 aPos;
    // layout (location = 1) in vec2 aTexCoord;

    -0.5f, -0.5f, -0.5f,  0.0f, 0.0f,
     0.5f, -0.5f, -0.5f,  1.0f, 0.0f,
     0.5f,  0.5f, -0.5f,  1.0f, 1.0f,
     0.5f,  0.5f, -0.5f,  1.0f, 1.0f,
    -0.5f,  0.5f, -0.5f,  0.0f, 1.0f,
    -0.5f, -0.5f, -0.5f,  0.0f, 0.0f,

    -0.5f, -0.5f,  0.5f,  0.0f, 0.0f,
     0.5f, -0.5f,  0.5f,  1.0f, 0.0f,
     0.5f,  0.5f,  0.5f,  1.0f, 1.0f,
     0.5f,  0.5f,  0.5f,  1.0f, 1.0f,
    -0.5f,  0.5f,  0.5f,  0.0f, 1.0f,
    -0.5f, -0.5f,  0.5f,  0.0f, 0.0f,

    -0.5f,  0.5f,  0.5f,  1.0f, 0.0f,
    -0.5f,  0.5f, -0.5f,  1.0f, 1.0f,
    -0.5f, -0.5f, -0.5f,  0.0f, 1.0f,
    -0.5f, -0.5f, -0.5f,  0.0f, 1.0f,
    -0.5f, -0.5f,  0.5f,  0.0f, 0.0f,
    -0.5f,  0.5f,  0.5f,  1.0f, 0.0f,

     0.5f,  0.5f,  0.5f,  1.0f, 0.0f,
     0.5f,  0.5f, -0.5f,  1.0f, 1.0f,
     0.5f, -0.5f, -0.5f,  0.0f, 1.0f,
     0.5f, -0.5f, -0.5f,  0.0f, 1.0f,
     0.5f, -0.5f,  0.5f,  0.0f, 0.0f,
     0.5f,  0.5f,  0.5f,  1.0f, 0.0f,

    -0.5f, -0.5f, -0.5f,  0.0f, 1.0f,
     0.5f, -0.5f, -0.5f,  1.0f, 1.0f,
     0.5f, -0.5f,  0.5f,  1.0f, 0.0f,
     0.5f, -0.5f,  0.5f,  1.0f, 0.0f,
    -0.5f, -0.5f,  0.5f,  0.0f, 0.0f,
    -0.5f, -0.5f, -0.5f,  0.0f, 1.0f,

    -0.5f,  0.5f, -0.5f,  0.0f, 1.0f,
     0.5f,  0.5f, -0.5f,  1.0f, 1.0f,
     0.5f,  0.5f,  0.5f,  1.0f, 0.0f,
     0.5f,  0.5f,  0.5f,  1.0f, 0.0f,
    -0.5f,  0.5f,  0.5f,  0.0f, 0.0f,
    -0.5f,  0.5f, -0.5f,  0.0f, 1.0f
};

unsigned int Cube::vertex_number = 180;
unsigned int Cube::VAO = 0;
unsigned int Cube::VBO = 0;
ShaderCube* Cube::shader = nullptr;
unsigned int Cube::texture_stone = 0;
unsigned int Cube::texture_dirt = 0;
unsigned int Cube::texture_diamond_ore = 0;
unsigned int Cube::texture_obsidian = 0;
int Cube::instance_count = 0;
bool Cube::static_buffer_initialized = false;

void Cube::init_static_buffer() {
    if (static_buffer_initialized) {
        printf("Warning: static buffer already initialized, cannot initialize.\n");
        return;
    }

    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);
    glBindVertexArray(VAO);
    glBindBuffer(GL_ARRAY_BUFFER, VBO);// 将全局指针GL_ARRAY_BUFFER的内容设为VBO的地址
    glBufferData(GL_ARRAY_BUFFER, sizeof(float) * Cube::vertex_number, Cube::vertices, GL_STATIC_DRAW);// 数据拷贝到显存。通过全局指针修改VBO的内容
    // set attributes: aPos and aTexCoord
    // 0-attribute index, 3-array length, GL_FLOAT-type, GL_FALSE-normalize, 5*sizeof(float)-step, (void*)0-bias
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 5 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);// enable this attribute manually, because it's disabled by default
    glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 5 * sizeof(float), (void*)(3 * sizeof(float)));
    glEnableVertexAttribArray(1);

    // load textures
    shader = new ShaderCube("./shader/cube.vert", "./shader/cube.frag");
    //shader->loadTexture(texture1, "./img/container.jpg", GL_RGB, GL_RGB);
    shader->loadTexture(texture_stone, "./img/stone.png", GL_RGB, GL_RGBA, false);
    shader->loadTexture(texture_dirt, "./img/dirt.png", GL_RGB, GL_RGBA, false);
    shader->loadTexture(texture_diamond_ore, "./img/diamond_ore.png", GL_RGB, GL_RGBA, false);
    shader->loadTexture(texture_obsidian, "./img/obsidian.png", GL_RGB, GL_RGBA, false);
    // tell opengl for each sampler to which texture unit it belongs to (only has to be done once)
    shader->use();
    shader->setInt("texture1", 0);// glsl: uniform sampler2D texture1;
    shader->setInt("texture2", 1);

    static_buffer_initialized = true;
}

void Cube::delete_static_buffer() {
    if (static_buffer_initialized) {
        glDeleteVertexArrays(1, &VAO);
        glDeleteBuffers(1, &VBO);
        delete shader;
    }
    else {
        printf("Warning: static buffer not initialized, cannot delete.\n");
    }
}


void Cube::init() {
    if (instance_count == 0)init_static_buffer();
    instance_count++; 
}

void Cube::render(const glm::mat4& projection, const glm::mat4& view) const {
    // render cubes
    glActiveTexture(GL_TEXTURE0);
    switch (block_type) {
    case _type_stone:
        glBindTexture(GL_TEXTURE_2D, texture_stone);
        break;
    case _type_dirt:
        glBindTexture(GL_TEXTURE_2D, texture_dirt);
        break;
    case _type_diamond_ore:
        glBindTexture(GL_TEXTURE_2D, texture_diamond_ore);
        break;
    default:
        glBindTexture(GL_TEXTURE_2D, texture_obsidian);
    }
    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_2D, texture_diamond_ore);
    glBindVertexArray(VAO);
    ShaderCube& shader = *Cube::shader;
    shader.use();
    shader.setFloat("mix_weight", mix_weight);
    shader.setMat4("projection", projection);// assign "uniform" variables in shader file
    shader.setMat4("view", view);
    
    glm::mat4 model = glm::mat4(1.0f);
    model = glm::scale(model, scale);
    model = glm::translate(model, position);
    model = glm::rotate(model, glm::radians(rotate_angle), rotate_axis);
    shader.setMat4("model", model);
    glDrawArrays(GL_TRIANGLES, 0, 36);
}

void Cube::cleanup() {
    instance_count--;
    if (instance_count == 0)delete_static_buffer();
}
