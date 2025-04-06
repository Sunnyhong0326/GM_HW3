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

  //Normalize Mesh
  static void normalizeMesh(Polyhedron* poly, double scale = 1.0);

  // Ray Testing
  static bool rayIntersectsTriangle(Eigen::Vector3f &rayOrigin, Eigen::Vector3f &rayDirection,
                                    Eigen::Vector3f &v0, Eigen::Vector3f &v1, Eigen::Vector3f &v2,
                                    Eigen::Vector3f &out);

  // Gemotry 
  static bool isNonManifoldVert(Vertex* vert);
  static bool isNonManifoldEdge(Edge* edge);

  static bool findHoles(Polyhedron* poly, std::vector<std::vector<int>>& holes);

  static void calcVertNormals(Polyhedron* poly);
  static void calcVertArea(Polyhedron* poly);
  static void calcFaceNormalsAndArea(Polyhedron* poly);
  static void calcEdgeLength(Polyhedron* poly);
  static void calcInteriorAngle(Polyhedron* poly);
  static void calcDihedralAngle(Polyhedron* poly);
  static double calcVolume(Polyhedron* poly);
  static double calcTotalFaceArea(Polyhedron* poly);
  static double calcTotalVertexArea(Polyhedron* poly);

  static int calcEulerCharacteristic(Polyhedron* poly);
  static double calcAngleDeficit(Vertex* vert);
  static double calcTotalAngleDeficit(Polyhedron* poly);
  static int calcValenceDeficit(Vertex* vert);
  static int calcTotalValenceDeficit(Polyhedron* poly);

  static void calcGaussCurvature(Polyhedron* poly);
  static Eigen::Vector3d calcMeanCurvatureNormal(Vertex* vert);
  static void calcMeanCurvature(Polyhedron* poly);
  static void calcPrincipalCurvature(Polyhedron* poly);

  static void calcCurvatureTensor(Polyhedron* poly);
  static void calcVertLocalframe(Vertex* vert, Eigen::Vector3d& local_u, Eigen::Vector3d& local_v);

};
