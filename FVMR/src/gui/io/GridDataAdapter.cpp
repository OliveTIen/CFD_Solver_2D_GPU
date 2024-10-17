#include "GridDataAdapter.h"
#include "UGridReader.h"
using namespace GUI;

void GridDataAdapter::setPointer(std::vector<glm::vec3>& vertices, std::vector<glm::uvec2>& indices,
	std::vector<std::pair<std::string, std::vector<glm::uvec2>>>& boundaries) {
	p_vertices = &vertices;
	p_indices = &indices;
	p_boundaries = &boundaries;
}

void GridDataAdapter::readData(const std::string& filepath) {
	UGridReader::readData(filepath, true);
	if (p_vertices == nullptr || p_indices == nullptr) {
		printf("Error: has not set pointer. @GridDataAdapter::readData\n");
		exit(getchar());
	}
	const auto& nodes = UGridReader::getData()->nodes;
	const auto& elements = UGridReader::getData()->elements;
	std::vector<glm::vec3>& vertices = *p_vertices;
	std::vector<glm::uvec2>& indices = *p_indices;
	auto& boundaries = *p_boundaries;
	vertices.resize(nodes.size());
	for (size_t i = 0; i < nodes.size(); i++) {
		vertices[i].x = nodes[i].x;
		vertices[i].z = nodes[i].y;// exchange y z
		vertices[i].y = 0;
	}
	//// 绘制时实际上是按照每两个点组成一条边来绘制的，不一定非要四个四个地去填充
	// 只是对于结构网格，绘制正方形网格比较好看，所以按照四个四个地填充
	// 1--2    12,23; 31,12
	// | /
	// |/
	// 3
	for (size_t i = 0; i < elements.size(); i++) {
		// triangles default
		indices.push_back(glm::uvec2(elements[i].nodes[0], elements[i].nodes[1]));
		indices.push_back(glm::uvec2(elements[i].nodes[1], elements[i].nodes[2]));
		indices.push_back(glm::uvec2(elements[i].nodes[2], elements[i].nodes[0]));
	}

	boundaries = UGridReader::getData()->boundaries;
	
}
