#pragma once

#include <vector>
#include <string>
#include <map>
#include <unordered_map>
#include <memory>
#include "Control.h"

class Object;
class Scene {
public:
	std::map<std::string, std::vector<Object*>> objectLists;
	
	void init(Control* control);
	void update();
	void render() const;
	void cleanup();

private:
	Control* m_control = nullptr;

};

