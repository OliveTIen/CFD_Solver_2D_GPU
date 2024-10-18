#include "CudaGrid.h"
#include "CudaGrid_kernels.cuh"
#include "../../solvers/SolverDataGetter.h"
#include "../../solvers/GPUSolver2.h"
#include "../io/GridDataAdapter.h"
#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <glm/gtc/type_ptr.hpp>
#include <cuda_gl_interop.h>



void CudaGrid::init() {
	GPU::GPUSolver2* pSolver = SolverDataGetter::getSolverInstance();
	if (!pSolver->isDeviceDataInitialized()) {
		std::cout << "Error: Device data not initialized. @CudaGrid::init()\n";
		return;
	}

    static_assert(sizeof(myfloat) == sizeof(float), "opengl requires single float type");

	// initialize shader
    shader_line = new ShaderBasic(vert_line_path, frag_line_path);
    shader_triangle = new ShaderBasic(vert_tri_path, frag_tri_path);


	// initialize host data
	initRenderingData_UnstructuredGrid();
	initPhysicsPointData_UnstructuredGrid();
	initPhysicsEdgeData_UnstructuredGrid();

    // initialize device colormap data
    std::vector< std::pair<float3, float>> colormap{
        { {1.0f, 0.0f, 0.0f}, 0.00f },
        { {1.0f, 1.0f, 0.0f}, 0.25f },
        { {0.0f, 1.0f, 0.0f}, 0.50f },
        { {0.0f, 0.0f, 1.0f}, 0.75f },
        { {1.0f, 0.0f, 1.0f}, 1.00f }
    };
    colormap_host.alloc(colormap.size());
    pSolver->setGPUDevice();
    colormap_device.cuda_alloc(colormap.size());
    nodeField_device.cuda_alloc(pSolver->node_device.num_node);
    getLastCudaError("CudaGrid device data alloc failed.");
    for (int i = 0; i < colormap.size(); i++) {
        colormap_host.data[i].x = colormap[i].first.x;
        colormap_host.data[i].y = colormap[i].first.y;
        colormap_host.data[i].z = colormap[i].first.z;
        colormap_host.data[i].w = colormap[i].second;
    }
    GPU::ColorMap::cuda_memcpy(&colormap_device, &colormap_host, cudaMemcpyHostToDevice);
    getLastCudaError("GPU::ColorMap::cuda_memcpy failed.");


    // initialize opengl data
    if (enable_render_line) {
        // data buffer of line
        glGenVertexArrays(1, &vao_line);
        glBindVertexArray(vao_line);

        glGenBuffers(1, &vbo_line);
        glBindBuffer(GL_ARRAY_BUFFER, vbo_line);
        glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(glm::vec3), glm::value_ptr(vertices[0]), GL_DYNAMIC_DRAW);
        glEnableVertexAttribArray(0);
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, nullptr);

        glGenBuffers(1, &ibo_line);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo_line);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices_line.size() * sizeof(glm::uvec2), glm::value_ptr(indices_line[0]), GL_STATIC_DRAW);

        glBindVertexArray(0);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
    }

    glGenVertexArrays(1, &vao_triangle);
    glBindVertexArray(vao_triangle);

    glGenBuffers(1, &vbo_triangle);
    glBindBuffer(GL_ARRAY_BUFFER, vbo_triangle);
    glBufferData(GL_ARRAY_BUFFER, vertices_triangle.size() * sizeof(float), vertices_triangle.data(), GL_DYNAMIC_DRAW);

    // position attribute. location=0, size=3, type=float, is_normalized=false, stride=6, offset=3. Ref: https://www.zaiwen.top/
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)0);
    // color attribute
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)(3 * sizeof(float)));

    glGenBuffers(1, &ibo_triangle);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo_triangle);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices_triangle.size() * sizeof(int), indices_triangle.data(), GL_STATIC_DRAW);

    glBindVertexArray(0);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);


    // set flag
    m_initialized = true;
}

void CudaGrid::reset() {
    if (!m_initialized)return;
    reset_position();
}

void CudaGrid::update(float dt) {
    // delay init
    if (!m_initialized) {
        init();
        if (!m_initialized)return;
    }

    ///////////////////////////////////////
    // physics
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


    if (enable_render_line) {
        // lines: update GPU buffer of vbo
        glBindBuffer(GL_ARRAY_BUFFER, vbo_line);
        glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(glm::vec3), glm::value_ptr(vertices[0]), GL_DYNAMIC_DRAW);
    }

    // triangles
    setTriangleVertexPosition(vertices_triangle, vertices);
    if (!render_using_device_data) setTriangleVertexColor();
    //cudaGraphicsUnregisterResource(cuda_vbo_resource);
    //getLastCudaError("cudaGraphicsUnregisterBuffer failed");
    glBindBuffer(GL_ARRAY_BUFFER, vbo_triangle);
    glBufferData(GL_ARRAY_BUFFER, vertices_triangle.size() * sizeof(float), vertices_triangle.data(), GL_DYNAMIC_DRAW);
    //cudaGraphicsGLRegisterBuffer(&cuda_vbo_resource, vbo_triangle, cudaGraphicsMapFlagsNone);
    //getLastCudaError("cudaGraphicsGLRegisterBuffer failed");
    if (render_using_device_data) setTriangleVertexColor_GPU();
}

void CudaGrid::render(const glm::mat4& projection, const glm::mat4& view) const {
    // GL_TRIANGLES ref: https://cloud.tencent.com/developer/article/1472451

    if (!m_initialized)return;

    if (enable_render_line) {
        glBindVertexArray(vao_line);
        shader_line->use();
        shader_line->setMat4("projection", projection);
        shader_line->setMat4("view", view);
        shader_line->setMat4("model", matrix_model);
        shader_line->setVec4("userColor", userColor);
        glLineWidth(2.0f);
        glDrawElements(GL_LINES, (GLuint)indices_line.size() * sizeof(glm::uvec2) / sizeof(int), GL_UNSIGNED_INT, NULL);

    }

    glBindVertexArray(vao_triangle);
    shader_triangle->use();
    shader_triangle->setMat4("projection", projection);
    shader_triangle->setMat4("view", view);
    shader_triangle->setMat4("model", matrix_model);
    shader_triangle->setFloat("alpha", render_blend_alpha);
    if (enable_render_blend) {
        glEnable(GL_BLEND);// transparent
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    }
    if (polygonmode_line) glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    glDrawElements(GL_TRIANGLES, (GLuint)indices_triangle.size(), GL_UNSIGNED_INT, NULL);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    glDisable(GL_BLEND);
    glBindVertexArray(0);

}

void CudaGrid::cleanup() {
    if (!m_initialized) {
        printf("Cannot cleanup. Not initialized.\n");
        return;
    }

    if (enable_render_line) {
        glDeleteVertexArrays(1, &vao_line);
        glDeleteBuffers(1, &vbo_line);
    }
    //cudaGraphicsUnregisterResource(cuda_vbo_resource);
    glDeleteVertexArrays(1, &vao_triangle);
    glDeleteBuffers(1, &vbo_triangle);

    delete shader_line;
    delete shader_triangle;

    cudaFree(colormap_device.data);
    cudaFree(nodeField_device.data);
    getLastCudaError("CudaGrid cudaFree failed.");
}


void CudaGrid::reset_position() {
    vertices = vertices_copy;
    setTriangleVertexPosition(vertices_triangle, vertices);
}

void CudaGrid::initRenderingData_UnstructuredGrid() {
    /*
    This function initializes vertices, vertices_copy, indices_line, boundaries
    */
    GPU::GPUSolver2* pSolver = SolverDataGetter::getSolverInstance();
    const auto& element_host = pSolver->element_host;
    const auto& node_host = pSolver->node_host;
    size_t num_element = element_host.num_element;
    size_t num_node = node_host.num_node;

    // vertices, vertices_copy of line
    vertices.resize(num_node);
    for (size_t i = 0; i < num_node; i++) {
        vertices[i].x = node_host.xy[0][i];
        vertices[i].z = node_host.xy[1][i];// exchange y z
        vertices[i].y = 0;
    }
    vertices_copy = vertices;

    // vertices_triangle
    setTriangleVertexPosition(vertices_triangle, vertices);

    // indices of line and triangle
    //// an edge contains two points, and a triangle contains three edges
    // 绘制时按照每两个点组成一条边来绘制
    // 1--2    12,23; 31,12
    // | /
    // |/
    // 3
    for (size_t i = 0; i < num_element; i++) {
        const auto& nodes = element_host.nodes;
        indices_line.push_back(glm::uvec2(nodes[0][i], nodes[1][i]));
        indices_line.push_back(glm::uvec2(nodes[1][i], nodes[2][i]));
        indices_line.push_back(glm::uvec2(nodes[2][i], nodes[0][i]));

        //indices_triangle.push_back(glm::uvec3(nodes[0][i], nodes[1][i], nodes[2][i]));
        indices_triangle.push_back(nodes[0][i]);
        indices_triangle.push_back(nodes[1][i]);
        indices_triangle.push_back(nodes[2][i]);
    }

    // boundary
    const auto& boundaryEdgeSets = pSolver->boundary_host_new.edgeSets;
    for (const auto& edgeSet : boundaryEdgeSets) {
        auto boundary = std::pair<std::string, std::vector<glm::uvec2>>(edgeSet.name, 0);
        boundaries.push_back(boundary);
        for (const auto& edgeID : edgeSet.edge_vector) {
            const auto& edge_host = pSolver->edge_host;
            int vertex0 = edge_host.nodes[0][edgeID];
            int vertex1 = edge_host.nodes[1][edgeID];
            auto& boundary = boundaries[boundaries.size() - 1];
            boundary.second.push_back(glm::uvec2(vertex0, vertex1));
        }
    }

    // triangle color
    for (size_t i = 0; i < vertices.size(); i++) {
        //float input = float(i) / float(12);
        //vertices_triangle[i * 2 + 1] = { 0.7f,0.7f,0.7f };
        vertices_triangle[i * 6 + 3] = 0.7f;
        vertices_triangle[i * 6 + 4] = 0.7f;
        vertices_triangle[i * 6 + 5] = 0.0f;
    }




}

void CudaGrid::initPhysicsPointData_UnstructuredGrid() {
    /*
    This function initializes velocities, X_hat and G;
    */
    velocities.resize(vertices.size());
    X_hat.resize(vertices.size());
    G.resize(vertices.size());
    for (auto& v : velocities) {
        v *= 0;
    }
}

void CudaGrid::initPhysicsEdgeData_UnstructuredGrid() {
    /*
    This function initializes E and L.
    */
    GPU::GPUSolver2* pSolver = SolverDataGetter::getSolverInstance();
    const auto& edge_host = pSolver->edge_host;
    size_t num_edge = edge_host.num_edge;
    E.resize(num_edge * 2);
    L.resize(num_edge);
    for (size_t edgeID = 0; edgeID < num_edge; edgeID++) {
        myint vertex0 = edge_host.nodes[0][edgeID];
        myint vertex1 = edge_host.nodes[1][edgeID];
        E[edgeID * 2 + 0] = vertex0;
        E[edgeID * 2 + 1] = vertex1;

        const auto& X = vertices;
        L[edgeID] = glm::distance(X[vertex0], X[vertex1]);
    }
}

void CudaGrid::Get_Gradient(float dt) {
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

void CudaGrid::Set_Boundary_Condition_Fixed(size_t index) {
    vertices[index] = X_hat[index];
    velocities[index] *= 0;
}

void CudaGrid::setTriangleVertexPosition(std::vector<float>& vertices_triangle, const std::vector<glm::vec3>& vertices) {
    if (vertices_triangle.empty()) {
        vertices_triangle.resize(vertices.size() * 6);
    }
    for (size_t i = 0; i < vertices.size(); i++) {
        vertices_triangle[i * 6 + 0] = vertices[i].x;
        vertices_triangle[i * 6 + 1] = vertices[i].y;
        vertices_triangle[i * 6 + 2] = vertices[i].z;
    }
}

void CudaGrid::setTriangleVertexColor() {
    auto pSolver = SolverDataGetter::getSolverInstance();
    int num_element = pSolver->element_host.num_element;
    enum ContourType {
        _pure_color,
        _color_bar
    };
    ContourType contour_type = ContourType::_color_bar;
    
    if (contour_type == _pure_color) {
        for (size_t i = 0; i < vertices.size(); i++) {
            glm::vec3 color = { 0.7f,0.7f,0.0f };
            vertices_triangle[i * 6 + 3] = color.x;
            vertices_triangle[i * 6 + 4] = color.y;
            vertices_triangle[i * 6 + 5] = color.z;
        }
    }
    else if (contour_type == _color_bar) {
        int num_triangle_test = num_element;
        for (size_t i = 0; i < num_triangle_test; i++) {
            for (int j = 0; j < 3; j++) {
                int nodeID = indices_triangle[i * 3 + j];

                float input = float(i) / float(num_triangle_test);
                input = vertices_triangle[nodeID * 6 + 0];
                glm::vec3 color = getColor(colormap_host, input);

                vertices_triangle[nodeID * 6 + 3] = color.x;
                vertices_triangle[nodeID * 6 + 4] = color.y;
                vertices_triangle[nodeID * 6 + 5] = color.z;
            }
        }
    }
}

void CudaGrid::setTriangleVertexColor_GPU() {
    GPU::setVerticesColor(vbo_triangle, cuda_vbo_resource, this);
}

glm::vec3 CudaGrid::getColor(const GPU::ColorMap& colormap, float input) {
    // Ref: zaiwen.top
    if (input < 0.0f) input = 0.0f;
    if (input > 1.0f) input = 1.0f;

    // 查找区间
    for (size_t i = 0; i < colormap.num_control_point - 1; ++i) {
        if (input >= colormap.data[i].w && input <= colormap.data[i + 1].w) {
            // get current interval data
            const glm::vec3& color0 = float4_to_vec3(colormap.data[i]);
            const glm::vec3& color1 = float4_to_vec3(colormap.data[i + 1]);
            float t0 = colormap.data[i].w;
            float t1 = colormap.data[i + 1].w;

            float normalizedInput = (input - t0) / (t1 - t0);

            // Catmull-Rom interpolation 三次插值
            glm::vec3 color =
                (1 - normalizedInput) * (1 - normalizedInput) * (1 - normalizedInput) * color0 +
                3 * (1 - normalizedInput) * (1 - normalizedInput) * normalizedInput * (color0 + color1) / 2.0f +
                3 * (1 - normalizedInput) * normalizedInput * normalizedInput * (color0 + color1) / 2.0f +
                normalizedInput * normalizedInput * normalizedInput * color1;
            return color;
        }
    }

    // default color #fc03fd = rgb(252,3,253). Purple color in Unity
    return glm::vec3(252.0f, 3.0f, 253.0f) / 255.0f;
    
}

