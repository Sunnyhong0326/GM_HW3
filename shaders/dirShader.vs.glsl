#version 330

layout (location = 0) in vec3 aPos;

uniform mat4 projection_matrix;	// projection matrix
uniform mat4 view_matrix;	// camera viewing transformation matrix
uniform vec3 trans;
uniform vec3 scale;
uniform vec3 dir;

void main() {
   vec3 lx = normalize(dir);
   vec3 ly = vec3(0.0f, 1.0f, 0.0f);
   if (1.0f - abs(dot(lx, ly)) < 1e-6) { ly = vec3(0.0f, 0.0f, 1.0f); }
   vec3 lz = cross(lx, ly);
   lz = normalize(lz);
   ly = cross(lx, lz);
   ly = normalize(ly);
   mat4 model_matrix = mat4(scale.x * lx.x, scale.x * lx.y, scale.x * lx.z, 0,
						    scale.y * ly.x, scale.y * ly.y, scale.y * ly.z, 0,
						    scale.z * lz.x, scale.z * lz.y, scale.z * lz.z, 0,
						    trans.x,trans.y,trans.z,1);
   gl_Position = projection_matrix * view_matrix * model_matrix * vec4(aPos, 1.0);
}
