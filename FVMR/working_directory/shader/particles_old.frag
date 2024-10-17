// Copyright (C) 2019 Xiao Zhai
// 
// This file is part of CPP-Fluid-Particles.
// 
// CPP-Fluid-Particles is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// CPP-Fluid-Particles is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with CPP-Fluid-Particles.  If not, see <http://www.gnu.org/licenses/>.

#version 120


varying float fs_pointSize;// 顶点球半径。源自顶点着色器
varying vec3 fs_PosEye;// 相机空间的坐标
varying mat4 u_Persp;// 投影矩阵

void main(void)
{
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
    vec4 pixelPos = vec4(fs_PosEye + normalize(N)*fs_pointSize,1.0f);
    vec4 clipSpacePos = u_Persp * pixelPos;
    //gl_FragDepth = clipSpacePos.z / clipSpacePos.w;
    
	// gl_FragColor = vec4(vec3(0.7f),1.0f);// 原，单一颜色
    // gl_FragColor = vec4(gl_Color.rgb,1.0f);// 颜色由压力决定
    gl_FragColor = vec4(exp(-mag*mag)*gl_Color.rgb,1.0f);// 添加边缘阴影，增强立体感
}
