#pragma once
#include "learnply.h"
#include <string>
#include <vector>

namespace MeshProcessor {
  //Normalize Mesh
  void normalizeMesh(Polyhedron* poly, double scale = 1.0);

  // Ray Testing
  bool rayIntersectsTriangle(Eigen::Vector3f &rayOrigin, Eigen::Vector3f &rayDirection,
                             Eigen::Vector3f &v0, Eigen::Vector3f &v1, Eigen::Vector3f &v2,
                             Eigen::Vector3f &out);

  // Gemotry 
  bool isNonManifoldVert(Vertex* vert);
  bool isNonManifoldEdge(Edge* edge);

  bool findHoles(Polyhedron* poly, std::vector<std::vector<int>>& holes);

  void calcVertNormals(Polyhedron* poly);
  void calcVertArea(Polyhedron* poly);
  void calcFaceNormalsAndArea(Polyhedron* poly);
  void calcEdgeLength(Polyhedron* poly);
  void calcInteriorAngle(Polyhedron* poly);
  void calcDihedralAngle(Polyhedron* poly);
  double calcVolume(Polyhedron* poly);
  double calcTotalFaceArea(Polyhedron* poly);
  double calcTotalVertexArea(Polyhedron* poly);

  int calcEulerCharacteristic(Polyhedron* poly);
  double calcAngleDeficit(Vertex* vert);
  double calcTotalAngleDeficit(Polyhedron* poly);
  int calcValenceDeficit(Vertex* vert);
  int calcTotalValenceDeficit(Polyhedron* poly);

  void calcGaussCurvature(Polyhedron* poly);
  Eigen::Vector3d calcMeanCurvatureNormal(Vertex* vert);
  void calcMeanCurvature(Polyhedron* poly);
  void calcPrincipalCurvature(Polyhedron* poly);

  void calcCurvatureTensor(Polyhedron* poly);
  void calcVertLocalframe(Vertex* vert, Eigen::Vector3d& local_u, Eigen::Vector3d& local_v);

};
