#version 330

uniform vec3 color;
out vec4 fragColor;

void main() {
	fragColor = vec4(color.xyz, 1.0);
}
