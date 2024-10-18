#include "CudaGrid_kernels.cuh"
#include "../../gpu/GPUGlobalFunction.h"
#include "../../solvers/SolverDataGetter.h"
#include "../../solvers/GPUSolver2.h"
#include <cuda_gl_interop.h>
#include <cuda_runtime.h>
#include "StructColorMap.h"
#include "../../math/ReduceGPU.h"
#include "../../global/GlobalPara.h"
#include "CudaGrid.h"



__device__ float3 operator*(float a, const float3&b) {
    return { a * b.x,a * b.y,a * b.z };
}

__device__ float3 operator+(const float3& a, const float3& b) {
    return { a.x + b.x, a.y + b.y, a.z + b.z };
}

__device__ float3 getColor_device(GPU::ColorMap colormap, float input) {
    /*
    input: must between 0 and 1
    Ref: zaiwen.top
    */
    if (input < 0.0f) input = 0.0f;
    if (input > 1.0f) input = 1.0f;

    for (size_t i = 0; i < colormap.num_control_point - 1; ++i) {
        if (input >= colormap.data[i].w && input <= colormap.data[i + 1].w) {
            // get current interval data
            float3 color0{}, color1{}, color_mean{};
            float t0{}, t1{};
            color0.x = colormap.data[i].x;
            color0.y = colormap.data[i].y;
            color0.z = colormap.data[i].z;
            t0 = colormap.data[i].w;
            color1.x = colormap.data[i + 1].x;
            color1.y = colormap.data[i + 1].y;
            color1.z = colormap.data[i + 1].z;
            t1 = colormap.data[i + 1].w;
            color_mean = 0.5 * (color0 + color1);

            float normalizedInput = (input - t0) / (t1 - t0);
            const float& n = normalizedInput;
            float omn = 1 - n;// one_minus_n
            float omn2 = omn * omn;
            float omn3 = omn2 * omn;
            float n2 = n * n;
            float n3 = n2 * n;

            // Catmull-Rom interpolation
            float3 color =
                omn3 * color0 +
                3 * omn2 * n * color_mean +
                3 * omn * n2 * color_mean +
                n3 * color1;
            return color;
            
        }
    }

    // default color #fc03fd = rgb(252,3,253). Purple color in Unity
    return { 252.0f / 255.0f, 3.0f / 255.0f, 253.0f / 255.0f };

}

__global__ void resetVerticesColor_k(float* vertices, GPU::NodeSoA node_device) {
    // reset color to zero, preparing to be accumulated
    const int iNode = blockIdx.x * blockDim.x + threadIdx.x;
    if (iNode >= node_device.num_node || iNode < 0) return;

    vertices[iNode * 6 + 3] = 0.0f;
    vertices[iNode * 6 + 4] = 0.0f;
    vertices[iNode * 6 + 5] = 0.0f;
}

__global__ void accumulateElementDataToNodeData_k(float* element_data, float* node_data, GPU::ElementSoA element_device, GPU::NodeSoA node_device) {
    /*
    thread id: element id
    element_data: array, length = num_element
    node_data: array,  length = num_node
    todo: 没有保证加法的原子性，求和会出问题
    */
    const int iElement = blockIdx.x * blockDim.x + threadIdx.x;
    if (iElement >= element_device.num_element || iElement < 0) return;

    for (int j = 0; j < 3; j++) {
        const int jNode = element_device.nodes[j][iElement];
        float input_var = element_data[iElement];
        //node_data[jNode] += input_var;
        input_var *= node_device.neighbor_element_count[jNode];
        float value = node_data[jNode];
        node_data[jNode] = (value > input_var) ? value : input_var;
        //atomicAdd(&node_data[jNode], input_var);
    }

    //__syncthreads();
}

__global__ void divideNodeData_k(float* node_data, GPU::NodeSoA node_device) {
    const int iNode = blockIdx.x * blockDim.x + threadIdx.x;
    if (iNode >= node_device.num_node || iNode < 0) return;

    float div = 1.0 / node_device.neighbor_element_count[iNode];
    if (div > 1e5)div = 1.0f;// avoid being divided by 0
    node_data[iNode] *= div;
}

__global__ void nodeDataToVertexColor_k(float* node_data, float* vertices, float2 input_range, GPU::ColorMap colormap, GPU::NodeSoA node_device) {
    const int iNode = blockIdx.x * blockDim.x + threadIdx.x;
    if (iNode >= node_device.num_node || iNode < 0) return;

    float input_var = node_data[iNode];

    // normalize
    float interval_length = abs(input_range.y - input_range.x);
    if (interval_length == 0)input_var = 0;// avoid being divided by zero
    else {
        input_var = (input_var - input_range.x) / interval_length;
    }

    float3 color = getColor_device(colormap, input_var);
    vertices[iNode * 6 + 3] = color.x;
    vertices[iNode * 6 + 4] = color.y;
    vertices[iNode * 6 + 5] = color.z;
}

__global__ void accumulateVerticesColor_k(float* vertices, float* dev_input, float2 input_range, GPU::ColorMap colormap, GPU::ElementSoA element_device, GPU::ElementFieldSoA elementField_device) {
    /*
    Function: add element value to vertices color
    dev_input: array, length = num_element
    */

    const int iElement = blockIdx.x * blockDim.x + threadIdx.x;
    if (iElement >= elementField_device.num || iElement < 0) return;

    float rho_ele = elementField_device.ruvp[0][iElement];
    for (int j = 0; j < 3; j++) {
        const int jNode = element_device.nodes[j][iElement];
        float input_var = dev_input[iElement];
        float interval_length = abs(input_range.y - input_range.x);
        if (interval_length == 0)input_var = 0;
        else {
            input_var = (input_var - input_range.x) / interval_length;
        }
        float3 color = getColor_device(colormap, input_var);

        vertices[jNode * 6 + 3] = color.x;
        vertices[jNode * 6 + 4] = color.y;
        vertices[jNode * 6 + 5] = color.z;
    }

}

__global__ void divideVerticesColor_k(float* vertices, GPU::NodeSoA node_device) {
    const int iNode = blockIdx.x * blockDim.x + threadIdx.x;
    if (iNode >= node_device.num_node || iNode < 0) return;

    float div = 1.0 / node_device.neighbor_element_count[iNode];
    if (div > 1e5)div = 1.0f;// avoid being divided by 0
    vertices[iNode * 6 + 3] *= div;
    vertices[iNode * 6 + 4] *= div;
    vertices[iNode * 6 + 5] *= div;

}

float2 getDeviceValueRange(myfloat* dev_input) {
    /*
    获取数组的最大、最小值。输入数组长度必须为num_element，因为临时数组大小有限
    规约运算。此处借用elementFieldVariable_dt_device的临时数组dev_output存储规约结果
    */
    myfloat out_min{}, out_max{};
    
    auto pSolver = SolverDataGetter::getSolverInstance();
    auto& elementField_device = pSolver->elementField_device;
    auto& elementFieldVariable_dt_device = pSolver->elementFieldVariable_dt_device;
    myint num_reduce = elementFieldVariable_dt_device.num_reduce;
    myfloat* dev_output = elementFieldVariable_dt_device.dev_output;

    GPU::Math::reduce_device(num_reduce, dev_input, dev_output, false, GPU::Math::ReduceType::reduceType_min);
    cudaMemcpy(&out_min, dev_output, sizeof(myfloat), cudaMemcpyDeviceToHost);
    getLastCudaError("get min value failed.");

    GPU::Math::reduce_device(num_reduce, dev_input, dev_output, false, GPU::Math::ReduceType::reduceType_max);
    cudaMemcpy(&out_max, dev_output, sizeof(myfloat), cudaMemcpyDeviceToHost);
    getLastCudaError("get max value failed.");

    return{ out_min, out_max };
}

void GPU::setVerticesColor(GLuint vbo, cudaGraphicsResource* cuda_vbo_resource, CudaGrid* cudaGrid) {
    auto pSolver = SolverDataGetter::getSolverInstance();
    int num_element = pSolver->element_device.num_element;
    int num_node = pSolver->node_device.num_node;
    
    // block and grid
    int block_size = GPU::MY_BLOCK_SIZE;
    int grid_size_ele = (num_element + block_size - 1) / block_size;
    int grid_size_node = (num_node + block_size - 1) / block_size;
    dim3 block(block_size, 1, 1);
    dim3 grid_ele(grid_size_ele, 1, 1);
    dim3 grid_node(grid_size_node, 1, 1);

    // get opengl resources
    cudaGraphicsGLRegisterBuffer(&cuda_vbo_resource, vbo, cudaGraphicsMapFlagsNone);
    getLastCudaError("cudaGraphicsGLRegisterBuffer failed");

    cudaGraphicsMapResources(1, &cuda_vbo_resource, 0);
    getLastCudaError("cudaGraphicsMapResources failed");

    float* vertices = nullptr;
    size_t num_bytes;
    cudaGraphicsResourceGetMappedPointer((void**)&vertices, &num_bytes, cuda_vbo_resource);
    getLastCudaError("cudaGraphicsResourceGetMappedPointer failed");

    // modify opengl vertices color
    {
        float* element_data = pSolver->elementField_device.ruvp[0];
        float* node_data = cudaGrid->nodeField_device.data;
        float2 input_range{ GlobalPara::render::range_1,GlobalPara::render::range_2 };
        
        cudaMemset(node_data, 0, num_node * sizeof(myfloat));
        getLastCudaError("reset node data failed.");

        accumulateElementDataToNodeData_k <<< grid_ele, block >>> (element_data, node_data, pSolver->element_device, pSolver->node_device);
        getLastCudaError("accumulate node data failed.");

        divideNodeData_k <<< grid_node, block >>> (node_data, pSolver->node_device);
        getLastCudaError("divide node data failed.");

        resetVerticesColor_k <<< grid_node, block >>> (vertices, pSolver->node_device);
        getLastCudaError("reset vertices color failed.");

        nodeDataToVertexColor_k <<< grid_node, block >>> (node_data, vertices, input_range, cudaGrid->colormap_device, pSolver->node_device);
        getLastCudaError("node data to vertex failed.");
    }



    //{
    //    // 
    //    resetVerticesColor_k <<< grid_node, block >>> (vertices, pSolver->node_device);
    //    getLastCudaError("reset vertices color failed.");

    //    float* dev_input = pSolver->elementField_device.ruvp[0];
    //    //float2 input_range = getDeviceValueRange(dev_input);
    //    //getLastCudaError("get dev_input range failed.");
    //    float2 input_range{ 0,2 };

    //    accumulateVerticesColor_k <<< grid_ele, block >>> (vertices, dev_input, input_range, *colormap, pSolver->element_device, pSolver->elementField_device);
    //    getLastCudaError("accumulate vertices color failed.");

    //    divideVerticesColor_k <<< grid_node, block >>> (vertices, pSolver->node_device);
    //    getLastCudaError("divide vertices color failed.");
    //}

    // release opengl resources
    cudaGraphicsUnmapResources(1, &cuda_vbo_resource, 0);
    getLastCudaError("cudaGraphicsUnmapResources failed");

    cudaGraphicsUnregisterResource(cuda_vbo_resource);
    getLastCudaError("cudaGraphicsUnregisterBuffer failed");
}