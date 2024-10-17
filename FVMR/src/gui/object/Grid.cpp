#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <glm/gtc/type_ptr.hpp>
#include "Grid.h"
#include "../io/GridDataAdapter.h"

int Grid::instance_count = 0;

void Grid::initRenderingData_StructuredGrid() {
    int n = slices + 1;// num of nodes per dimension
    vertices.resize(n * n);
    reset_position();

    for (int j = 0; j < slices; ++j) {
        for (int i = 0; i < slices; ++i) {
            int row1 = j * (slices + 1);
            int row2 = (j + 1) * (slices + 1);
            // 我不知道为什么要用uvec4，但经过列举法验证，确实一次push_back了两条边
            //  1----2      row1+i,  row1+i+1
            //  |    |
            //  4----3      row2+i,  row2+i+1
            // line 12, line 23
            indices.push_back(glm::uvec2(row1 + i, row1 + i + 1));
            indices.push_back(glm::uvec2(row1 + i + 1, row2 + i + 1));
            // line 34, line 41
            indices.push_back(glm::uvec2(row2 + i + 1, row2 + i));
            indices.push_back(glm::uvec2(row2 + i, row1 + i));
        }
    }
    // sizeof(glm::uvec4)/sizeof(int)=4, 4 means each uvec4 has 4 ints 表示每个uvec4有4个整型数
    index_length = (GLuint)indices.size() * sizeof(glm::uvec2)/sizeof(int);// 绘制直线时，是2个点2个点地绘制的。例如index_length个点，会绘制index_length/2条边
}

void Grid::initPhysicsPointData_StructuredGrid() {

    ///////////////////
    // 顶点物理变量
    velocities.resize(vertices.size());
    X_hat.resize(vertices.size());
    G.resize(vertices.size());
    for (auto& v : velocities) {
        v *= 0;
    }

    
}

void Grid::initPhysicsEdgeData_StructuredGrid() {
    //////////////////
    // 边物理变量
    // 三角形所有顶点的集合
    std::vector<int> triangles;
    triangles.resize(slices * slices * 6);
    size_t t = 0;
    int n = slices + 1;
    for (int j = 0; j < slices; ++j) {
        for (int i = 0; i < slices; ++i) {
            triangles[t * 6 + 0] = j * n + i;
            triangles[t * 6 + 1] = j * n + i + 1;
            triangles[t * 6 + 2] = (j + 1) * n + i + 1;
            triangles[t * 6 + 3] = j * n + i;
            triangles[t * 6 + 4] = (j + 1) * n + i + 1;
            triangles[t * 6 + 5] = (j + 1) * n + i;
            t++;
        }
    }
    // 构造_E，是三角形所有边的集合
    std::vector<int> _E;// edge for all triangles, storing node indices
    _E.resize(triangles.size() * 2);
    for (size_t i = 0; i < triangles.size(); i += 3) {// 每次记录1个三角形的3条边。i每次+3地递增，因为1个三角形占用3个元素
        _E[i * 2 + 0] = triangles[i + 0];
        _E[i * 2 + 1] = triangles[i + 1];
        _E[i * 2 + 2] = triangles[i + 1];
        _E[i * 2 + 3] = triangles[i + 2];
        _E[i * 2 + 4] = triangles[i + 2];
        _E[i * 2 + 5] = triangles[i + 0];
    }
    // reorder the list
    for (int i = 0; i < _E.size(); i += 2)
        if (_E[i] > _E[i + 1])
            Swap(_E[i], _E[i + 1]);
    // sort the original edge list using quicksort
    Quick_Sort(_E.data(), 0, _E.size() / 2 - 1);

    // 构造边数组E，计算边长度L
    int e_number = 0;
    for (int i = 0; i < _E.size(); i += 2)
        if (i == 0 || _E[i + 0] != _E[i - 2] || _E[i + 1] != _E[i - 1])
            e_number++;

    E.resize(e_number * 2);
    for (int i = 0, e = 0; i < _E.size(); i += 2)
        if (i == 0 || _E[i + 0] != _E[i - 2] || _E[i + 1] != _E[i - 1]) {
            E[e * 2 + 0] = _E[i + 0];
            E[e * 2 + 1] = _E[i + 1];
            e++;
        }

    L.resize(E.size() / 2);
    const auto& X = vertices;
    for (int e = 0; e < E.size() / 2; e++) {
        int v0 = E[e * 2 + 0];
        int v1 = E[e * 2 + 1];
        L[e] = glm::distance(X[v0], X[v1]);
    }
}

void Grid::initRenderingData_UnstructuredGrid(std::string filepath) {
    GUI::GridDataAdapter adapter;
    adapter.setPointer(vertices, indices, boundaries);
    adapter.readData(filepath);
    // 绘制直线时，是2个点2个点地绘制的。例如index_length个点，会绘制index_length/2条边
    index_length = (GLuint)indices.size() * sizeof(glm::uvec2) / sizeof(int);

}

void Grid::initPhysicsPointData_UnstructuredGrid() {
    ///////////////////
    // 顶点物理变量
    velocities.resize(vertices.size());
    X_hat.resize(vertices.size());
    G.resize(vertices.size());
    for (auto& v : velocities) {
        v *= 0;
    }

}

void Grid::initPhysicsEdgeData_UnstructuredGrid() {
    //////////////////
    // 边物理变量
    // 三角形所有顶点的集合。长度为
    std::vector<int> triangles;
    triangles.resize(indices.size());
    for (size_t i = 0; i < indices.size(); i++) {
        triangles[i] = indices[i].x;
    }

    // 构造_E，是三角形所有边的集合
    std::vector<int> _E;// edge for all triangles, storing node indices
    _E.resize(triangles.size() * 2);
    for (size_t i = 0; i < triangles.size(); i += 3) {// 每次记录1个三角形的3条边。i每次+3地递增，因为1个三角形占用3个元素
        _E[i * 2 + 0] = triangles[i + 0];
        _E[i * 2 + 1] = triangles[i + 1];
        _E[i * 2 + 2] = triangles[i + 1];
        _E[i * 2 + 3] = triangles[i + 2];
        _E[i * 2 + 4] = triangles[i + 2];
        _E[i * 2 + 5] = triangles[i + 0];
    }
    // reorder the list
    for (int i = 0; i < _E.size(); i += 2)
        if (_E[i] > _E[i + 1])
            Swap(_E[i], _E[i + 1]);
    // sort the original edge list using quicksort
    Quick_Sort(_E.data(), 0, _E.size() / 2 - 1);

    // 构造边数组E，计算边长度L
    int e_number = 0;
    for (int i = 0; i < _E.size(); i += 2)
        if (i == 0 || _E[i + 0] != _E[i - 2] || _E[i + 1] != _E[i - 1])
            e_number++;

    E.resize(e_number * 2);
    for (int i = 0, e = 0; i < _E.size(); i += 2)
        if (i == 0 || _E[i + 0] != _E[i - 2] || _E[i + 1] != _E[i - 1]) {
            E[e * 2 + 0] = _E[i + 0];
            E[e * 2 + 1] = _E[i + 1];
            e++;
        }

    L.resize(E.size() / 2);
    const auto& X = vertices;
    for (int e = 0; e < E.size() / 2; e++) {
        int v0 = E[e * 2 + 0];
        int v1 = E[e * 2 + 1];
        L[e] = glm::distance(X[v0], X[v1]);
    }

}

void Grid::init() {
    /*
    https://stackoverflow.com/questions/58494179/how-to-create-a-grid-in-opengl-and-drawing-it-with-lines
    */
    instance_count++;
    if (buffer_allocated) {// not static_buffer
        printf("buffer already allocated.\n");
        return;
    }

    shader = new ShaderLine(vert_path, frag_path);

    if (b_structured) {
        initRenderingData_StructuredGrid();
        initPhysicsPointData_StructuredGrid();
        initPhysicsEdgeData_StructuredGrid();
    }
    else {
        initRenderingData_UnstructuredGrid(file_path);
        initPhysicsPointData_UnstructuredGrid();
        initPhysicsEdgeData_UnstructuredGrid();
        vertices_copy = vertices;
    }

    glGenVertexArrays(1, &vao);
    glBindVertexArray(vao);

    glGenBuffers(1, &vbo);
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(glm::vec3), glm::value_ptr(vertices[0]), GL_DYNAMIC_DRAW);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, nullptr);

    glGenBuffers(1, &ibo);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(glm::uvec2), glm::value_ptr(indices[0]), GL_STATIC_DRAW);

    glBindVertexArray(0);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);

    buffer_allocated = true;
}

void Grid::render(const glm::mat4& projection, const glm::mat4& view) const {
    glBindVertexArray(vao);
    shader->use();
    shader->setMat4("projection", projection);
    shader->setMat4("view", view);
    shader->setMat4("model", matrix_model);
    shader->setVec4("userColor", userColor);
    glLineWidth(2.0f);
    glDrawElements(GL_LINES, index_length, GL_UNSIGNED_INT, NULL);
    glBindVertexArray(0);
}

void Grid::cleanup() {
    if (!buffer_allocated) {
        printf("buffer not allocated, cannot delete buffer\n");
        return;
    }

    glDeleteVertexArrays(1, &vao);
    glDeleteBuffers(1, &vbo);
    delete shader;

    instance_count--;
}

void Grid::update(float dt) {
    float dt_inv = 1.0f / dt;
    float dt2_inv = 1.0f / (dt * dt);
    auto& X = vertices;
    auto& V = velocities;

    for (size_t i = 0; i < X.size(); i++) {
        V[i] *= velocity_damping;
        X_hat[i] = X[i] + dt * V[i];
        X[i] = X_hat[i];
    }

    for (int k_substep = 0; k_substep < substep_num; k_substep++) {
        Get_Gradient(dt);

        // Update X by gradient.
        // 遍历点，施加动量和重力 Momentum and Gravity
        for (int i = 0; i < X.size(); i++) {
            float hessian_inv = 1.0f / (dt2_inv * mass + 4 * spring_k);
            X[i] -= hessian_inv * G[i];
        }
    }

    // update V
    for (int i = 0; i < X.size(); i++) {
        V[i] += dt_inv * (X[i] - X_hat[i]);
    }


    for (auto& boundary : boundaries) {
        auto& edges = boundary.second;
        for (auto& edge : edges) {
            Set_Boundary_Condition_Fixed(edge.x);
            Set_Boundary_Condition_Fixed(edge.y);
        }
    }


    X_hat = vertices;

    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(glm::vec3), glm::value_ptr(vertices[0]), GL_DYNAMIC_DRAW);

}

void Grid::Get_Gradient(float dt) {
    float dt2_inv = 1.0f / (dt * dt);
    auto& X = vertices;
    auto& V = velocities;


    for (int i = 0; i < X.size(); i++) {
        G[i] = dt2_inv * mass * (X[i] - X_hat[i]);// momentum; 动量(惯性力) 
        G[i] -= gravity_vector;// gravity
    }

    for (int iedge = 0; iedge < E.size() / 2; iedge++) {
        int vi = E[iedge * 2 + 0];
        int vj = E[iedge * 2 + 1];
        float Le = L[iedge];// original length
        glm::vec3 xi_xj = X[vi] - X[vj];// real length
        float abs_xi_xj = glm::distance(X[vi], X[vj]);
        glm::vec3 spring_force = spring_k * (1.0f - Le / abs_xi_xj) * xi_xj;
        G[vi] += spring_force;
        G[vj] -= spring_force;
    }
}
