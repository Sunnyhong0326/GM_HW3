#pragma once

#include <glad/glad.h>
#include <iostream>
#include <memory>

#include "helpers/axeshelper.h"
#include "helpers/edgehelper.h"
#include "helpers/vertexhelper.h"
#include "scene.h"

class Renderer {
public:
  Renderer(int width, int height, Scene *scenePtr);
  ~Renderer();

  void resize(int width, int height);
  void render();

private:
	void setupShaders();
	void releaseShaders();

private:
  int screenWidth;
  int screenHeight;

  // State
  Scene *scene;

  //Shader
  Shader* plyShader;
  Shader* colorShader;
  Shader* dirShader;
  Shader* textShader;
  Shader* axesShader;

  // Helpers
  std::unique_ptr<VertexHelper> vertexHelper;
  std::unique_ptr<EdgeHelper> edgeHelper;
  std::unique_ptr<AxesHelper> axesHelper;
};
