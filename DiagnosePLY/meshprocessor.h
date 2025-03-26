#pragma once
#include "learnply.h"
#include <string>
#include <vector>

class MeshProcessor {
public:
  MeshProcessor() = delete;

  // Getters
  static std::vector<Triangle *> &getTriangles(Polyhedron *poly) { return poly->tlist; }
  static std::vector<Vertex *> &getVertices(Polyhedron *poly) { return poly->vlist; }
  static std::vector<Edge *> &getEdges(Polyhedron *poly) { return poly->elist; }
  static std::vector<Corner *> &getCorners(Polyhedron *poly) { return poly->clist; }

  // Ray Testing
  static bool rayIntersectsTriangle(Eigen::Vector3f &rayOrigin, Eigen::Vector3f &rayDirection,
                                    Eigen::Vector3f &v0, Eigen::Vector3f &v1, Eigen::Vector3f &v2,
                                    Eigen::Vector3f &out);

  // Gemotry 
  static bool isNonManifoldVert(Vertex* vert);
  static bool isNonManifoldEdge(Edge* edge);

  static bool detectHole(Polyhedron* poly, std::vector<std::vector<int>>& holes);

  static void calcVertNormals(Polyhedron* poly);
  static void calcFaceNormalsAndArea(Polyhedron* poly);
  static void calcEdgeLength(Polyhedron* poly);
  static void calcInteriorAngle(Polyhedron* poly);
  static void calcDihedralAngle(Polyhedron* poly);
  static void calcVertexArea(Polyhedron* poly);

  static double calcVolume(Polyhedron* poly);
  static int calcEulerCharacteristic(Polyhedron* poly);
  static double calcAngularDeficit(Vertex* vert);
  static double calcTotalAngularDeficit(Polyhedron* poly);

  static void calcGaussCurvature(Polyhedron* poly);
  static Eigen::Vector3d calcMeanCurvatureNormal(Vertex* vert);
  static void calcMeanCurvature(Polyhedron* poly);
  static void calcPrincipalCurvature(Polyhedron* poly);

  static void calcCurvatureTensor(Polyhedron* poly);
  static void calcVertLocalframe(Vertex* vert, Eigen::Vector3d& local_u, Eigen::Vector3d& local_v);

};
