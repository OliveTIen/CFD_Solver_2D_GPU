#version 330 core
layout (location = 0) in vec3 aPos;
layout (location = 1) in vec3 aColor;

out vec4 ourColor;

uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;
uniform float alpha;

void main()
{
	gl_Position = projection * view * model * vec4(aPos, 1.0);
	ourColor = vec4(aColor, alpha);
}