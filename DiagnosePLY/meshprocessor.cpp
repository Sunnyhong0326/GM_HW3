#include "meshprocessor.h"
#include "learnply.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

void MeshProcessor::normalizeMesh(Polyhedron* poly, double scale)
{
    Eigen::Vector3d centroid(0.0, 0.0, 0.0);
    for (int i = 0; i < poly->nverts(); i++) {
        centroid += poly->vlist[i]->pos;
    }
    centroid /= (double)poly->nverts();
    double max_dist = 0.0;
    for (int i = 0; i < poly->nverts(); i++) {
        double dist = (poly->vlist[i]->pos - centroid).squaredNorm();
        if (dist > max_dist) {
            max_dist = dist;
        }
    }
    max_dist = std::max(sqrt(max_dist), 1e-6);
    for (int i = 0; i < poly->nverts(); i++) {
        poly->vlist[i]->pos = centroid + scale * (poly->vlist[i]->pos - centroid) / max_dist;
    }
}

/******************************************************************************
Check if the given ray intersects the triangle

Entry:
  rayOrigin - origin of the ray
  rayDirection - direction of the ray
  v0,v1,v2 - three vertices of the triangle
  out - output of (u,v,t)

Exit:
  return true if the ray intersects the triangle; otherwise, return false
******************************************************************************/
#define EPSILON_RAY  1e-8f
bool MeshProcessor::rayIntersectsTriangle(Eigen::Vector3f &rayOrigin, Eigen::Vector3f &rayDirection,
                                          Eigen::Vector3f &v0, Eigen::Vector3f &v1, Eigen::Vector3f& v2,
                                          Eigen::Vector3f &out) {
  float u = 0.0f, v = 0.0f, t = 0.0f;

  Eigen::Vector3f edge1 = v1 - v0;
  Eigen::Vector3f edge2 = v2 - v0;
  Eigen::Vector3f h = rayDirection.cross(edge2);
  float a = edge1.dot(h);

  if (fabs(a) < EPSILON_RAY)
      return false;

  float f = 1.0f / a;
  Eigen::Vector3f s = rayOrigin - v0;
  u = f * s.dot(h);

  if (u < 0.0f || u > 1.0f)
      return false;

  Eigen::Vector3f q = s.cross(edge1);
  v = f * rayDirection.dot(q);

  if (v < 0.0f || u + v > 1.0f)
      return false;

  t = f * edge2.dot(q);
  out << u, v, t;
  return t > EPSILON_RAY;
}

/******************************************************************************
Calculate vertex normals for the given polyhedron

Entry:
  poly - pointer to the polyhedron

Exit:
  Vertex normals are calculated and stored in the polyhedron
******************************************************************************/
void MeshProcessor::calcVertNormals(Polyhedron* poly) {
	/// TODO: Weighted the face normals by the angle at the vertex or the face area
    for (int i = 0; i < poly->nverts(); i++) {
        poly->vlist[i]->normal = Eigen::Vector3d(0.0, 0.0, 0.0);
        for (int j = 0; j < poly->vlist[i]->ntris(); j++) {
            poly->vlist[i]->normal += poly->vlist[i]->tris[j]->normal;
        }
        poly->vlist[i]->normal.normalize();
    }
    ///
}

/******************************************************************************
Calculate face normals and area for the given polyhedron

Entry:
  poly - pointer to the polyhedron

Exit:
  Face normals and areas are calculated and stored in the triangles
******************************************************************************/
void MeshProcessor::calcFaceNormalsAndArea(Polyhedron* poly) {
    for (int i = 0; i < poly->ntris(); i++) {
        Vertex* v0 = poly->tlist[i]->verts[0];
        Vertex* v1 = poly->tlist[i]->verts[1];
        Vertex* v2 = poly->tlist[i]->verts[2];
        poly->tlist[i]->normal = (v2->pos - v0->pos).cross(v1->pos - v0->pos);
        poly->tlist[i]->area = poly->tlist[i]->normal.norm() * 0.5;
        poly->tlist[i]->normal.normalize();
    }
}

/******************************************************************************
Calculate edge lengths for the given polyhedron

Entry:
  poly - pointer to the polyhedron

Exit:
  Edge lengths are calculated and stored in the edges
******************************************************************************/
void MeshProcessor::calcEdgeLength(Polyhedron* poly) {
    for (int i = 0; i < poly->nedges(); i++) {
        Vertex* v0 = poly->elist[i]->verts[0];
        Vertex* v1 = poly->elist[i]->verts[1];
        poly->elist[i]->length = (v1->pos - v0->pos).norm();
    }
}

/******************************************************************************
Check if the given vertex is non-manifold

Entry:
  vert - pointer to the vertex

Exit:
  return true if the vertex is non-manifold; otherwise, return false
******************************************************************************/
bool MeshProcessor::isNonManifoldVert(Vertex* vert)
{
	/// Implement:
	/// 1. Check if the vertex is a cutting vertex
	/// 2. Check if the vertex is a dangling vertex
    return false;
}

/******************************************************************************
Check if the given edge is non-manifold

Entry:
  edge - pointer to the edge

Exit:
  return true if the edge is non-manifold; otherwise, return false
******************************************************************************/
bool MeshProcessor::isNonManifoldEdge(Edge* edge)
{
    /// Implement:
    /// 1. Check if the vertex is a cutting edge
    /// 2. Check if the vertex is a dangling edge (no incident face)
    return false;
}

/******************************************************************************
Detect holes in the given polyhedron

Entry:
  poly - pointer to the polyhedron
  holes - reference to a vector of vectors to store hole indices

Exit:
  return true if holes are detected; otherwise, return false
******************************************************************************/
bool MeshProcessor::detectHole(Polyhedron* poly, std::vector<std::vector<int>>& holes)
{
    holes.clear();
    /// Implement:
	/// 1. Mark an edge if the edge is on the boundary
    /// 2. Connect the marked edges into holes
	/// 3. Store the edges indices in the vector (i.e., std::vector<int>)
    /// 4. Store the hole in the vector (i.e., std::vector<std::vector<int>>)
    return !holes.empty();
}

/******************************************************************************
Calculate interior angles for the given polyhedron

Entry:
  poly - pointer to the polyhedron

Exit:
  Interior angles are calculated and stored in the corners
******************************************************************************/
void MeshProcessor::calcInteriorAngle(Polyhedron* poly) {
    for (int i = 0; i < poly->ncorners(); i++) 
    {
        Corner* c = poly->clist[i];
        /// Implement:
		/// Calculate the interior angle of the corner
        /// Note use tan2
        double interior_angle = 0.0;
        ///
        c->interior_angle = interior_angle;
    }
}

/******************************************************************************
Calculate dihedral angles for the given polyhedron

Entry:
  poly - pointer to the polyhedron

Exit:
  Dihedral angles are calculated and stored in the edges
******************************************************************************/
void MeshProcessor::calcDihedralAngle(Polyhedron* poly) {
    for (int i = 0; i < poly->nedges(); i++) {
        Edge* e = poly->elist[i];
        /// Implement:
        /// Calculate the dihedral angle of the corner
        /// Note use tan2
        double dihedral_angle = 0.0;
        ///
        e->dihedral_angle = dihedral_angle;
    }
}

/******************************************************************************
Calculate vertex areas for the given polyhedron

Entry:
  poly - pointer to the polyhedron

Exit:
  Vertex areas are calculated and stored in the vertices
******************************************************************************/
void MeshProcessor::calcVertArea(Polyhedron* poly) {
    for (int i = 0; i < poly->nverts(); i++) {
        Vertex* vert_i = poly->vlist[i];
        /// Implement:
        /// Calculate the vertex area
        double area = 0.0;
        ///
        vert_i->area = area;
    }
}

/******************************************************************************
Calculate the average vertex degree for the given polyhedron

Entry:
  poly - pointer to the polyhedron

Exit:
  return the average vertex degree of the polyhedron
******************************************************************************/
double MeshProcessor::calcAvgVertDegree(Polyhedron* poly)
{
	double avgDegree = 0.0;
    for (int i = 0; i < poly->nverts(); i++) {
        Vertex* vert_i = poly->vlist[i];
        /// Implement:
        /// Count the number of incient edges of the vertex
        avgDegree += 0.0;
        ///
    }
    return avgDegree;
}

/******************************************************************************
Calculate the volume of the given polyhedron

Entry:
  poly - pointer to the polyhedron

Exit:
  return the calculated volume of the polyhedron
******************************************************************************/
double MeshProcessor::calcVolume(Polyhedron* poly) {
    /// Implement:
    /// Calculate the volume
    double volume = 0.0;
    ///
    return volume;
}

/******************************************************************************
Calculate the total face area of the given polyhedron

Entry:
  poly - pointer to the polyhedron

Exit:
  return the total face area of the polyhedron
******************************************************************************/
double MeshProcessor::calcTotalFaceArea(Polyhedron* poly)
{
    double area = 0.0;
    for (int i = 0; i < poly->ntris(); i++) {
        area += poly->tlist[i]->area;
    }
    return area;
}

/******************************************************************************
Calculate the total vertex area of the given polyhedron

Entry:
  poly - pointer to the polyhedron

Exit:
  return the total vertex area of the polyhedron
******************************************************************************/
double MeshProcessor::calcTotalVertexArea(Polyhedron* poly)
{
    double area = 0.0;
    for (int i = 0; i < poly->nverts(); i++) {
        area += poly->vlist[i]->area;
    }
    return area;
}

/******************************************************************************
Calculate the Euler characteristic of the given polyhedron

Entry:
  poly - pointer to the polyhedron

Exit:
  return the calculated Euler characteristic of the polyhedron
******************************************************************************/
int MeshProcessor::calcEulerCharacteristic(Polyhedron* poly) {
    /// Implement:
	/// Calculate the Euler characteristic
    int euler = 0;
    ///
    return euler;
}

/******************************************************************************
Calculate the angular deficit of the given vertex

Entry:
  vert - pointer to the vertex

Exit:
  return the calculated angular deficit of the vertex
******************************************************************************/
double MeshProcessor::calcAngularDeficit(Vertex* vert) {
    /// Implement:
	/// Calculate the angular deficit of the vertex
    double deficit = 0.0;
    /// 
    return deficit;
}

/******************************************************************************
Calculate the total angular deficit of the given polyhedron

Entry:
  poly - pointer to the polyhedron

Exit:
  return the total angular deficit of the polyhedron
******************************************************************************/
double MeshProcessor::calcTotalAngularDeficit(Polyhedron* poly) {
    /// Implement:
	/// Calculate the total angular deficit of the polyhedron
    double totalAngularDeficit = 0.0;
    ///
    return totalAngularDeficit;
}

/******************************************************************************
Calculate the Gaussian curvature of the given polyhedron

Entry:
  poly - pointer to the polyhedron

Exit:
  Gaussian curvature is calculated and stored in the vertices
******************************************************************************/
void MeshProcessor::calcGaussCurvature(Polyhedron* poly) {
    for (int i = 0; i < poly->nverts(); i++) {
        Vertex* vert_i = poly->vlist[i];
        /// Implement:
		/// Calculate the Gaussian curvature of the vertex
        double gauss = 0.0;
        ///
        vert_i->gaussCurvature = gauss;
    }
}

/******************************************************************************
Calculate the mean curvature normal of the given vertex

Entry:
  vert - pointer to the vertex

Exit:
  return the calculated mean curvature normal of the vertex
******************************************************************************/
Eigen::Vector3d MeshProcessor::calcMeanCurvatureNormal(Vertex* vert)
{
    /// Implement 
	/// Calculate the mean curvature normal of the vertex
    Eigen::Vector3d normal(0.0, 0.0, 0.0);
    ///
    return normal;
}

/******************************************************************************
Calculate the mean curvature of the given polyhedron

Entry:
  poly - pointer to the polyhedron

Exit:
  Mean curvature is calculated and stored in the vertices
******************************************************************************/
void MeshProcessor::calcMeanCurvature(Polyhedron* poly) {
    for (int i = 0; i < poly->nverts(); i++) {
        Vertex* vert_i = poly->vlist[i];
        Eigen::Vector3d normal = calcMeanCurvatureNormal(vert_i);
        /// Implement:
		/// Calculate the mean curvature of the vertex
        double meanCurvature = 0.0;
        ///
        vert_i->meanCurvature = meanCurvature;
    }
}

/******************************************************************************
Calculate the principal curvature of the given polyhedron

Entry:
  poly - pointer to the polyhedron

Exit:
  Principal curvature is calculated and stored in the vertices
******************************************************************************/
void MeshProcessor::calcPrincipalCurvature(Polyhedron* poly) {
    for (int i = 0; i < poly->nverts(); i++) {
        Vertex* vert_i = poly->vlist[i];
        /// Implement:
		/// Calculate the principal curvature of the vertex
        double k1 = 0.0;
        double k2 = 0.0;
        ///
        vert_i->minPrincCurvature = k1;
        vert_i->maxPrincCurvature = k2;
    }
}

/******************************************************************************
Calculate the curvature tensor of the given polyhedron

Entry:
  poly - pointer to the polyhedron

Exit:
  Curvature tensor and principal direction are calculated and stored in the vertex
******************************************************************************/
void MeshProcessor::calcCurvatureTensor(Polyhedron* poly)
{
    for (int i = 0; i < poly->nverts(); i++) {
        Vertex* vi = poly->vlist[i];
        //Least Square Fitting
        Eigen::MatrixXd matA(vi->corners.size(), 3);
        Eigen::VectorXd k(vi->corners.size()); //kij
        Eigen::Vector3d e1, e2;
        calcVertLocalframe(vi, e1, e2);
        /// Implement:
        /// 1. Assign A and k
        /// 2. Slove Ax = k
		/// 3. Compute Principal Direction
        //vi->tensor(0, 0) = l;
        //vi->tensor(1, 1) = n;
        //vi->tensor(0, 1) = m;
        //vi->tensor(1, 0) = m;
        //vi->princDir2D[0] = v1;
        //vi->princDir2D[1] = v2;
        //vi->princDir3D[0] = d1;
        //vi->princDir3D[1] = d2;
    }
}

/******************************************************************************
Calculate the local frame of the given vertex

Entry:
  vi - pointer to the vertex
  local_u - reference to store the local u vector
  local_v - reference to store the local v vector

Exit:
  Local frame is calculated and stored in the references
******************************************************************************/
void MeshProcessor::calcVertLocalframe(Vertex* vi, Eigen::Vector3d& e1, Eigen::Vector3d& e2) {
    if (vi->corners.size() == 0) { return; }
    // Find edge with minimum projection distance to the vertex normal
    double  minProj = DBL_MAX;
    Vertex* minVj = NULL;
    for (int j = 0; j < vi->corners.size(); j++)
    {
        Vertex* vj = vi->corners[j]->next->vertex;
		/// Implement:
		/// 1. Calculate the projection distance
		/// 2. Find the edge with the minimum projection distance
        /// 
    }
    // Calculate the local frame perpendicular to the vertex normal 
    if (minVj != NULL)
    {
        /// Implement:
		/// Calculate the local frame with vi and minVj
        /// 
    }
    else
    {
		e1 << 1.0, 0.0, 0.0;
        e2 = vi->normal.cross(e1);
        e1 = vi->normal.cross(e2);
    }
}