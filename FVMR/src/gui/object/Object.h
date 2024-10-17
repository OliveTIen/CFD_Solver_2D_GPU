#pragma once

class Object {
public:
	virtual ~Object() {};
	virtual void init() = 0;
	virtual void reset() = 0;
	virtual void update(float dt) = 0;
	virtual void render(const glm::mat4& projection, const glm::mat4& view) const = 0;
	virtual void cleanup() = 0;
};