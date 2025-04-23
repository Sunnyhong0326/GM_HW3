#include "meshprocessor.h"
#include "learnply.h"
#include <unordered_set>
#include <queue>
#include <iostream>

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
        Vertex* vert = poly->vlist[i];

        // Initialize the normal as zero
        vert->normal = Eigen::Vector3d(0.0, 0.0, 0.0);

        // Loop through all adjacent faces (triangles) to the vertex
        for (Triangle* tri : vert->tris) {
            if (!tri) continue;

            // Calculate the face normal (triangle normal)
            Eigen::Vector3d v0 = tri->corners[0]->vertex->pos;
            Eigen::Vector3d v1 = tri->corners[1]->vertex->pos;
            Eigen::Vector3d v2 = tri->corners[2]->vertex->pos;

            Eigen::Vector3d face_normal = (v1 - v0).cross(v2 - v0).normalized();

            // Optionally weight by face area or angle
            double area = 0.5 * (v1 - v0).cross(v2 - v0).norm();
            vert->normal += area * face_normal;  // Weighting by area

            // Or weight by angle (if desired)
            // double angle_weight = calcAngleAtVertex(vert, tri);
            // vert->normal += angle_weight * face_normal;
        }

        // Normalize the vertex normal
        vert->normal.normalize();
    }
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
    if (!vert) return false;

    // 1. Dangling vertex: fewer than 2 triangles
    //if (vert->tris.size() < 2)
    //    return true;

    // 2. Cutting vertex: check if connected triangles form one fan
    // Use corner connectivity to traverse the triangle fan

    std::unordered_set<Triangle*> visited;
    std::queue<Corner*> q;

    // Start with the first corner
    // if (vert->corners.empty()) return true;  // No corners: invalid
    q.push(vert->corners[0]);
    visited.insert(vert->corners[0]->tri);

    while (!q.empty()) {
        Corner* c = q.front(); q.pop();

        // Look around the fan: prev and next corners lead to adjacent corners
        for (Corner* dir : { c->prev, c->next }) {
            if (dir && dir->vertex == vert) {
                Corner* opp = dir->oppsite;
                if (opp && opp->vertex == vert) {
                    Triangle* neighborTri = opp->tri;
                    if (visited.count(neighborTri) == 0) {
                        visited.insert(neighborTri);
                        q.push(opp);
                    }
                }
            }
        }
    }

    // If we didn't visit all triangles the vertex belongs to, it's a cutting vertex
    if (visited.size() != vert->tris.size())
        return true;

    return false; // It's manifold
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
bool MeshProcessor::findHoles(Polyhedron* poly, std::vector<std::vector<int>>& holes)
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
    for (int i = 0; i < poly->ncorners(); i++) {
        Corner* c = poly->clist[i];

        // Get the two edges at the corner (each corner belongs to one edge)
        Vertex* v0 = c->vertex;               // Vertex of the corner
        Edge* edge1 = c->next->edge;                 // First edge at the corner
        Edge* edge2 = c->next->next->edge;           // Next edge at the corner (adjacent corner's edge)

        // Calculate the direction vectors of the two edges
        Eigen::Vector3d vector1 = edge1->verts[1]->pos - v0->pos;  // Direction of edge 1
        Eigen::Vector3d vector2 = edge2->verts[1]->pos - v0->pos;  // Direction of edge 2

        // Calculate the interior angle between the two vectors
        double cross_product = vector1.cross(vector2).norm();  // |u x v| (magnitude of cross product)
        double dot_product = vector1.dot(vector2);            // u . v (dot product)

        // Calculate the angle using atan2 (cross product magnitude, dot product)
        double interior_angle = std::atan2(cross_product, dot_product);

        // Store the calculated interior angle in the corner
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

        // Get the two triangles (faces) that share the edge
        Triangle* t1 = e->tris[0];  // First triangle sharing the edge
        Triangle* t2 = e->tris[1];  // Second triangle sharing the edge

        // Get the vertices of the two faces that share the edge
        Vertex* v1 = e->verts[0];  // First vertex of the edge
        Vertex* v2 = e->verts[1];  // Second vertex of the edge

        // Compute the normal vector for each triangle (face)
        Eigen::Vector3d normal1 = (t1->verts[1]->pos - t1->verts[0]->pos).cross(t1->verts[2]->pos - t1->verts[0]->pos);
        Eigen::Vector3d normal2 = (t2->verts[1]->pos - t2->verts[0]->pos).cross(t2->verts[2]->pos - t2->verts[0]->pos);

        // Normalize the normal vectors
        normal1.normalize();
        normal2.normalize();

        // Compute the edge vector (from v1 to v2)
        Eigen::Vector3d edge_vector = v2->pos - v1->pos;
        Eigen::Vector3d e_hat = edge_vector.normalized();  // Unit vector along the edge

        // Calculate the cross product between the normals
        Eigen::Vector3d cross_product = normal1.cross(normal2);

        // Calculate the dot product between the normals
        double dot_product = normal1.dot(normal2);

        // Calculate the dihedral angle using atan2
        double dihedral_angle = std::atan2(cross_product.dot(e_hat), dot_product); // Angle in radians

        // Store the dihedral angle in the edge
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
        Vertex* vi = poly->vlist[i];
        double area = 0.0;
        /*
        
        for (Corner* corner_i : vi->corners) {
            Corner* cornerbeta = corner_i->next;
            Corner* corneralpha = corner_i->next->oppsite;
            Edge* eij = corner_i->next->edge;
            Triangle* tri = corner_i->tri;

            double angle_i = corner_i->interior_angle;
            double angle_alpha = corneralpha->interior_angle;
            double angle_beta = cornerbeta->interior_angle;

            double cot_beta = 1 / tan(angle_beta); // cot(£])
            double cot_alpha = 1 / tan(angle_alpha); // cot(alpha)

            bool is_obtuse = (angle_i > M_PI / 2 || angle_alpha > M_PI / 2 || angle_beta > M_PI / 2);

            if (is_obtuse) {
                if (angle_i > M_PI / 2) {
                    area += tri->area / 2.0;
                }
                else {
                    area += tri->area / 4.0;
                }
            }
            else {
                area += (cot_alpha +  cot_beta) * eij->length / 8.0;
            }
        }
        */
        for (int j = 0; j < vi->ntris(); j++) {
            Triangle* tri = vi->tris[j];
            Vertex* v0 = tri->verts[0];
            Vertex* v1 = tri->verts[1];
            Vertex* v2 = tri->verts[2];

            // Identify which one is vi
            Vertex* vj, * vk;
            if (vi == v0) { vj = v1; vk = v2; }
            else if (vi == v1) { vj = v2; vk = v0; }
            else { vj = v0; vk = v1; }

            Eigen::Vector3d pij = vj->pos - vi->pos;
            Eigen::Vector3d pik = vk->pos - vi->pos;
            Eigen::Vector3d pjk = vk->pos - vj->pos;
            Eigen::Vector3d pkj = vj->pos - vk->pos;

            double angle_i = acos(pij.normalized().dot(pik.normalized()));
            double angle_j = acos((-pij).normalized().dot(pjk.normalized()));
            double angle_k = acos((-pik).normalized().dot(pkj.normalized()));

            bool is_obtuse = (angle_i > M_PI / 2 || angle_j > M_PI / 2 || angle_k > M_PI / 2);

            if (is_obtuse) {
                if (angle_i > M_PI / 2) {
                    area += tri->area / 2.0;
                }
                else {
                    area += tri->area / 4.0;
                }
            }
            else {
                double cot_beta = pik.dot(pjk) / pik.cross(pjk).norm(); // cot(£])
                double cot_gamma = pij.dot(-pjk) / pij.cross(-pjk).norm(); // cot(£^)
                area += (pij.squaredNorm() * cot_gamma + pik.squaredNorm() * cot_beta) / 8.0;
            }
        }
        
        vi->area = area;
    }
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
    if (!poly) return 0.0;

    double volume = 0.0;

    for (Triangle* tri : poly->tlist) {
        if (!tri ) continue;

        Eigen::Vector3d v0 = tri->corners[0]->vertex->pos;
        Eigen::Vector3d v1 = tri->corners[1]->vertex->pos;
        Eigen::Vector3d v2 = tri->corners[2]->vertex->pos;

        // Signed volume of tetrahedron formed with the origin
        double vol = v0.dot(v1.cross(v2)) / 6.0;
        volume += vol;
    }

    return std::abs(volume);  // Use absolute value to return positive volume
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
    if (!poly) return 0;

    int V = static_cast<int>(poly->vlist.size());
    int E = static_cast<int>(poly->elist.size());
    int F = static_cast<int>(poly->tlist.size());

    int euler = V - E + F;
    return euler;
    return euler;
}

/******************************************************************************
Calculate the angular deficit of the given vertex

Entry:
  vert - pointer to the vertex

Exit:
  return the calculated angular deficit of the vertex
******************************************************************************/
double MeshProcessor::calcAngleDeficit(Vertex* vert) {
    /// Implement:
	/// Calculate the angular deficit of the vertex
    if (!vert || vert->corners.empty())
        return 0.0;

    double angleSum = 0.0;

    for (Corner* c : vert->corners) {
        if (!c || !c->prev || !c->next) continue;

        double angle = c->interior_angle;
        angleSum += angle;
    }

    const double TWO_PI = 2.0 * M_PI;
    double deficit = TWO_PI - angleSum;
    return deficit;
}

/******************************************************************************
Calculate the total angular deficit of the given polyhedron

Entry:
  poly - pointer to the polyhedron

Exit:
  return the total angular deficit of the polyhedron
******************************************************************************/
double MeshProcessor::calcTotalAngleDeficit(Polyhedron* poly) {
    /// Implement:
	/// Calculate the total angular deficit of the polyhedron
    if (!poly) return 0.0;

    double total_deficit = 0.0;
            
    for (Vertex* vert : poly->vlist) {
        total_deficit += calcAngleDeficit(vert);  // use the function you implemented earlier
    }

    return total_deficit;
}

/******************************************************************************
Calculate the vertex valence deficit for the given vertex

Entry:
  vert - pointer to the vertex

Exit:
  return the average vertex valence of the polyhedron
******************************************************************************/
int MeshProcessor::calcValenceDeficit(Vertex* vert)
{
    /// Implement:
	/// Count the number of incient edges of the vertex and minus 6
    if (!vert) return 0;

    std::unordered_set<Vertex*> neighborVerts;

    for (Corner* c : vert->corners) {
        if (!c || !c->prev || !c->next) continue;

        // Add the two neighboring vertices from the triangle
        if (c->prev->vertex != vert)
            neighborVerts.insert(c->prev->vertex);
        if (c->next->vertex != vert)
            neighborVerts.insert(c->next->vertex);
    }

    int valence = static_cast<int>(neighborVerts.size());
    return valence - 6;
}

/******************************************************************************
Calculate the vertex valence for the given vertex

Entry:
  vert - pointer to the vertex

Exit:
  return the average vertex valence of the polyhedron
******************************************************************************/
int MeshProcessor::calcTotalValenceDeficit(Polyhedron* poly)
{
    /// Implement:
    /// Calculate the total vertex valence deficit of the polyhedron
    if (!poly) return 0;

    int total_deficit = 0;

    for (Vertex* vert : poly->vlist) {
        total_deficit += calcValenceDeficit(vert);
    }

    return total_deficit;
}

/******************************************************************************
Calculate the Gaussian curvature of the given polyhedron

Entry:
  poly - pointer to the polyhedron

Exit:
  Gaussian curvature is calculated and stored in the vertices
******************************************************************************/
void MeshProcessor::calcGaussCurvature(Polyhedron* poly) {
    if (!poly) return;
    for (int i = 0; i < poly->nverts(); i++) {
        Vertex* vert_i = poly->vlist[i];
        /// Implement:
		/// Calculate the Gaussian curvature of the vertex
        if (!vert_i) continue;

        double gauss = calcAngleDeficit(vert_i);  // Use angle deficit as curvature
        ///
        vert_i->gaussCurvature = gauss / vert_i->area;
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
    Eigen::Vector3d H(0.0, 0.0, 0.0);

    if (!vert || vert->corners.empty())
        return H;

    double mixedArea = 0.0;

    for (Corner* c : vert->corners) {
        Triangle* tri = c->tri;
        if (!tri) continue;

        // Vertices of the triangle
        Eigen::Vector3d vi = vert->pos;
        Eigen::Vector3d vj = c->next->vertex->pos;
        Eigen::Vector3d vk = c->prev->vertex->pos;

        // Vectors
        Eigen::Vector3d e0 = vj - vi;
        Eigen::Vector3d e1 = vk - vi;
        Eigen::Vector3d e2 = vk - vj;

        // Cotangent weights
        double cot_alpha = e0.dot(e2) / (e0.cross(e2)).norm(); // angle at vk
        double cot_beta = e1.dot(-e2) / (e1.cross(-e2)).norm(); // angle at vj

        // Contribution from edge (i, j) and (i, k)
        H += (cot_beta * e0 + cot_alpha * e1) * 0.5;

    }

    return  H /= (2.0 * vert->area);
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
        double meanCurvature = (normal.dot(vert_i->normal.normalized())) / (vert_i->area * 4);
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
        Eigen::VectorXd vecK(vi->corners.size());
        Eigen::Vector3d e1, e2;
        calcVertLocalframe(vi, e1, e2);
        /// Implement:
        /// 1. Assign A and k
        /// 2. Slove Ax = k
		/// 3. Compute Principal Direction
        /// 
        // Curvature Tensor
        //vi->tensor(0, 0) = l;
        //vi->tensor(1, 1) = n;
        //vi->tensor(0, 1) = m;
        //vi->tensor(1, 0) = m;
        // Eigenvectors
        //vi->princDir2D[0] = v1;
        //vi->princDir2D[1] = v2;
        // Principal Direction
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