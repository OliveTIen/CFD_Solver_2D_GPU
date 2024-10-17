#version 330 core
out vec4 FragColor;

in vec2 TexCoord;

// texture samplers
uniform sampler2D texture1;
uniform sampler2D texture2;
uniform float mix_weight;

void main()
{
	// linearly interpolate between both textures (80% container, mix_weight=20% awesomeface)
	FragColor = mix(texture(texture1, TexCoord), texture(texture2, TexCoord), mix_weight);
}