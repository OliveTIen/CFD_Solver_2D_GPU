#version 330 core

in vec4 ourColor;
in float fs_pointSize;// 顶点球半径。源自顶点着色器
in vec3 fs_PosEye;// 相机空间的坐标
in mat4 u_Persp;// 投影矩阵
out vec4 FragColor;

void main()
{
	FragColor = ourColor;

    // 根据圆还原法向量坐标(半球坐标)
    // 根据纹理坐标计算球面法向量，然后计算片段的三维坐标pixelPos
    // gl_PointCoord告诉你当前片段(fragment)相对于所属顶点的坐标(类似纹理坐标)
    vec3 N;
    N.xy = gl_PointCoord.xy*vec2(2.0, -2.0) + vec2(-1.0, 1.0);
    // 根据球半径公式求z坐标，z = sqrt(1-x*x-y*y) 
    // 如果把1.0改为2.0，会导致圆球贴图被正方形裁剪
    // 使用discard清除圆外的黑边，否则在正方形内、圆外会有黑边，给人的观感就是一堆正方形贴图在移动
    float mag = dot(N.xy, N.xy);
    if (mag > 1.0) discard;// 清除圆外的黑边
    N.z = sqrt(1.0-mag);

    //calculate depth
    // vec4 pixelPos = vec4(fs_PosEye + normalize(N)*fs_pointSize,1.0f);
    // vec4 clipSpacePos = u_Persp * pixelPos;// 
    // gl_FragDepth = clipSpacePos.z / clipSpacePos.w;

	// output color 三维球体阴影效果
    FragColor = vec4(exp(-mag*mag)*ourColor.rgb,1.0f);// 添加边缘阴影，增强立体感
}