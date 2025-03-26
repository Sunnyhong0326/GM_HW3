#pragma once
#include "../miscellaneous/shader.h"
#include <Eigen/Dense>
#include <glad/glad.h>
#include <vector>

class EdgeHelper {
public:
  EdgeHelper();
  ~EdgeHelper();

  void initialize(Shader* edge_shader);
  void use(const Eigen::Matrix4f& projectionMatrix, const Eigen::Matrix4f& viewMatrix);
  void draw(const Eigen::Vector3f &start, const Eigen::Vector3f &end, const Eigen::Vector3f& color);

private:
  GLuint VAO, VBO, EBO;
  Shader *edgeShader;
  size_t indicesCount;

  void createCylinder(float radius = 0.003f, int segments = 8);
};
