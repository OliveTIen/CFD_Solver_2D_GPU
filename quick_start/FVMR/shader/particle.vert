#version 330 core
layout (location = 0) in vec3 aPos;

out vec4 ourColor;
out float fs_pointSize;// 顶点球半径。源自顶点着色器
out vec3 fs_PosEye;// 相机空间的坐标
out mat4 u_Persp;// 投影矩阵

uniform mat4 projection;
uniform mat4 view;// 从世界坐标系变换到摄像机坐标系
uniform mat4 model;
uniform vec4 userColor;
uniform float pointRadius;// point size in world space
uniform float refDistance;// the distance where gl_PointSize = pointRadius

/*
opengl内置变量有
gl_PointSize：点渲染模式，渲染的方形像素大小
varying：旧版glsl使用的类型。现在用in out代替
*/

void main()
{
    vec4 view_position = view * model * vec4(aPos, 1.0);
    gl_Position = projection * view_position;

    float distance = length(view_position.xyz);// distance from particle to camera
    gl_PointSize = pointRadius * refDistance / distance;

	ourColor = userColor;
    fs_pointSize = gl_PointSize;
    fs_PosEye = view_position.xyz;
    u_Persp = projection;
}