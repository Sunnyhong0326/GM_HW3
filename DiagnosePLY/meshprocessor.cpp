#include "meshprocessor.h"
#include "learnply.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

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

void MeshProcessor::calcVertNormals(Polyhedron* poly) {
    /// TODO: Weighted 
    for (int i = 0; i < poly->nverts(); i++) {
        poly->vlist[i]->normal = Eigen::Vector3d(0.0, 0.0, 0.0);
        for (int j = 0; j < poly->vlist[i]->ntris(); j++) {
            poly->vlist[i]->normal += poly->vlist[i]->tris[j]->normal;
        }
        poly->vlist[i]->normal.normalize();
    }
    ///
}

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

void MeshProcessor::calcEdgeLength(Polyhedron* poly) {
    for (int i = 0; i < poly->nedges(); i++) {
        Vertex* v0 = poly->elist[i]->verts[0];
        Vertex* v1 = poly->elist[i]->verts[1];
        poly->elist[i]->length = (v1->pos - v0->pos).norm();
    }
}

bool MeshProcessor::isNonManifoldVert(Vertex* vert)
{
    /// TODO: Implement 
    ///
    return false;
}

bool MeshProcessor::isNonManifoldEdge(Edge* edge)
{
    /// TODO: Implement 
    ///
    return false;
}

bool MeshProcessor::detectHole(Polyhedron* poly, std::vector<std::vector<int>>& holes)
{
    holes.clear();
    /// TODO: Implement 
    ///
    return !holes.empty();
}

void MeshProcessor::calcInteriorAngle(Polyhedron* poly) {
    for (int i = 0; i < poly->ncorners(); i++) 
    {
        Corner* c = poly->clist[i];
        /// TODO: Implement 
        double interior_angle = 0.0;
        ///
        c->interior_angle = 0.0;
    }
}

void MeshProcessor::calcDihedralAngle(Polyhedron* poly) {
    for (int i = 0; i < poly->nedges(); i++) {
        Edge* e = poly->elist[i];
        /// TODO: Implement 
        double dihedral_angle = 0.0;
        ///
        e->dihedral_angle = dihedral_angle;
    }
}

void MeshProcessor::calcVertexArea(Polyhedron* poly) {
    for (int i = 0; i < poly->nverts(); i++) {
        Vertex* vert_i = poly->vlist[i];
        /// TODO: Implement 
        double area = 0.0;
        ///
        vert_i->area = area;
    }
}

double MeshProcessor::calcVolume(Polyhedron* poly) {
    /// TODO: Implement 
    double volume = 0.0;
    ///
    return volume;
}

double MeshProcessor::calcArea(Polyhedron* poly)
{
    double area = 0.0;
    for (int i = 0; i < poly->ntris(); i++) {
        area += poly->tlist[i]->area;
    }
    return area;
}

int MeshProcessor::calcEulerCharacteristic(Polyhedron* poly) {
    /// TODO: Implement 
    int euler = 0;
    ///
    return euler;
}

double MeshProcessor::calcAngularDeficit(Vertex* vert) {
    ///  TODO: Implement 
    double deficit = 0.0;
    /// 
    return deficit;
}

double MeshProcessor::calcTotalAngularDeficit(Polyhedron* poly) {
    /// TODO: Implement 
    double totalAngularDeficit = 0.0;
    ///
    return totalAngularDeficit;
}

void MeshProcessor::calcGaussCurvature(Polyhedron* poly) {
    for (int i = 0; i < poly->nverts(); i++) {
        Vertex* vert_i = poly->vlist[i];
        /// TODO: Implement 
        double gauss = 0.0;
        ///
        vert_i->gaussCurvature = gauss;
    }
}


Eigen::Vector3d MeshProcessor::calcMeanCurvatureNormal(Vertex* vert)
{
    /// TODO: Implement 
    ///
    return Eigen::Vector3d(0.0, 0.0, 0.0);
}

void MeshProcessor::calcMeanCurvature(Polyhedron* poly) {
    for (int i = 0; i < poly->nverts(); i++) {
        Vertex* vert_i = poly->vlist[i];
        Eigen::Vector3d meanCurvatureNormal = calcMeanCurvatureNormal(vert_i);
        // TODO: Implement 
        double meanCurvature = 0.0;
        //
        vert_i->meanCurvature = meanCurvature;
    }
}

void MeshProcessor::calcPrincipalCurvature(Polyhedron* poly) {
    for (int i = 0; i < poly->nverts(); i++) {
        Vertex* vert_i = poly->vlist[i];
        /// TODO: Implement 
        double k1 = 0.0;
        double k2 = 0.0;
        ///
        vert_i->minPrincCurvature = k1;
        vert_i->maxPrincCurvature = k2;
    }
}

void MeshProcessor::calcCurvatureTensor(Polyhedron* poly)
{
    for (int i = 0; i < poly->nverts(); i++) {
        Vertex* vi = poly->vlist[i];
        //Least Square Fitting
        Eigen::MatrixXd matA(vi->corners.size(), 3);
        Eigen::VectorXd k(vi->corners.size()); //kij
        Eigen::Vector3d local_u, local_v;
        calcVertLocalframe(vi, local_u, local_v);
        //Assign A and k
        for (int j = 0; j < vi->corners.size(); j++)
        {
            Vertex* vj = vi->corners[j]->next->vertex;
            /// TODO: Implement 
            //set matrix A
            matA(j, 0) = 0.0;
            matA(j, 1) = 0.0;
            matA(j, 2) = 0.0;
            //set matrix K
            k(j) = 0.0;
            ///
        }
        //Slove Ax = k
        {
            //Tensor: [l m; m n]
            Eigen::VectorXd sol = matA.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(k);
            double l = sol[0], m = sol[1], n = sol[2];
            // Find two eigenvectors of matrix([[l,m], [m,n]]). They are the two principal directions.
            vi->tensor(0, 0) = l;
            vi->tensor(1, 1) = n;
            vi->tensor(0, 1) = m;
            vi->tensor(1, 0) = m;
            double delta = sqrt((l - n) * (l - n) + 4.0 * (m * m));
            vi->princDir2D[0] = Eigen::Vector2d(m, (-l + n - delta) * 0.5);
            vi->princDir2D[1] = Eigen::Vector2d(m, (-l + n + delta) * 0.5);
            vi->princDir3D[0] = vi->princDir2D[0].x() * local_u + vi->princDir2D[0].y() * local_v; // convert back to 3D
            vi->princDir3D[1] = vi->princDir2D[1].x() * local_u + vi->princDir2D[1].y() * local_v; // convert back to 3D
            vi->princDir2D[0].normalize();
            vi->princDir2D[1].normalize();
            vi->princDir3D[0].normalize();
            vi->princDir3D[1].normalize();
        }
    }
}

void MeshProcessor::calcVertLocalframe(Vertex* vi, Eigen::Vector3d& local_u, Eigen::Vector3d& local_v) {
    if (vi->corners.size() == 0) { return; }
    // Find edge with minimum projection distance to the vertex normal
    double  minProj = DBL_MAX;
    Vertex* minV = NULL;
    for (int j = 0; j < vi->corners.size(); j++)
    {
        Vertex* vj = vi->corners[j]->next->vertex;
        Eigen::Vector3d vec = vj->pos - vi->pos;
        double proj = fabs(vec.dot(vi->normal));
        if (proj < minProj)
        {
            minProj = proj;
            minV = vj;
        }
    }
    // Calculate the local frame perpendicular to the vertex normal 
    if (minV != NULL)
    {
        Eigen::Vector3d vec = minV->pos - vi->pos;
        local_u = vec - vec.dot(vi->normal) * vi->normal;
        local_u.normalize();
        local_v = vi->normal.cross(local_u);
        local_v.normalize();
    }
}