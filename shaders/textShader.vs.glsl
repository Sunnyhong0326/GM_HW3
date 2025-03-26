#version 330 core
layout (location = 0) in vec2 vertex;
layout (location = 1) in vec2 texture_coord;

out vec2 TexCoords;

uniform mat4 projection_matrix;
uniform mat4 view_matrix;

// Selected Vertex Normal & Translation (a.k.a. Position)
uniform vec3 up;
uniform vec3 normal;
uniform vec3 translation;


void main()
{
    // tangent (T)
    vec3 tangent = normalize(cross(normal, -up));
    // binormal (B)
    vec3 binormal = normalize(cross(normal, tangent));
    // rotate
    vec3 r_vertex = mat3(tangent, binormal, normal) * vec3(vertex, 0.0);

    // final pos
    vec3 pos = r_vertex + translation + normal * 0.01;

    // apply final pos
    gl_Position = projection_matrix * view_matrix * vec4(pos, 1.0);

    TexCoords = texture_coord;
}
