#pragma once
#include "learnply.h"
#include <glad/glad.h>
#include "miscellaneous/tinycolormap.h"


class MeshRenderer {
public:
  MeshRenderer();
  ~MeshRenderer();

  void setupBuffers(Polyhedron *poly);

  void drawPLY();
  void drawWireframe();

  void setNormalMode(Polyhedron *poly, bool using_face_normal);
  void setColors(Polyhedron* poly, int mode, double max_value = 100.0);
  void updateColors(Polyhedron *poly);

private:
  // For regular rendering
  GLuint VAO, VBO, NBO, CBO;
  size_t vertexCount;

  // For wireframe rendering
  GLuint wireVAO, wireVBO;
  size_t edgeCount;
};
