#pragma once
#include <glm/glm.hpp>
#include <vector>
#include <string>

namespace GUI {
	class GridDataAdapter {
	public:
		std::vector<glm::vec3>* p_vertices = nullptr;
		std::vector<glm::uvec2>* p_indices = nullptr;
		std::vector<std::pair<std::string, std::vector<glm::uvec2>>>* p_boundaries = nullptr;

		void setPointer(std::vector<glm::vec3>& vertices, std::vector<glm::uvec2>& indices,
			std::vector<std::pair<std::string, std::vector<glm::uvec2>>>& boundaries);
		void readData(const std::string& filepath);
	};
}

