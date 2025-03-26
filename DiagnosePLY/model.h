#pragma once
#include <memory>
#include <string>
#include "learnply.h"
#include "meshrenderer.h"

class Model {
public:
  Model();
  Model(const std::string &path);

  // Setters for internal objects
  void setPolyhedron(Polyhedron *polyhedron);

  // Getters for internal objects
  Polyhedron *getPolyhedron() { return polyhedron.get(); }
  MeshRenderer *getMeshRenderer() { return meshRenderer.get(); }

private:
  std::unique_ptr<Polyhedron> polyhedron;
  std::unique_ptr<MeshRenderer> meshRenderer;
};
