#include <algorithm>
#include <vector>

#include <OpenMesh/Core/IO/MeshIO.hh>

#include "MyViewer.h"

MyViewer::MyViewer(QWidget *parent) : QGLViewer(parent),
                                      mean_min(0.0), mean_max(0.0), cutoff_ratio(0.05)
{
}

MyViewer::~MyViewer()
{
}

void MyViewer::updateMeanMinMax()
{
  size_t n = mesh.n_vertices();
  if(n == 0)
    return;

  std::vector<double> mean;
  mean.reserve(n);
  for(MyMesh::VertexIter i = mesh.vertices_begin(), ie = mesh.vertices_end(); i != ie; ++i)  
    mean.push_back(mesh.data(i).mean);

  std::sort(mean.begin(), mean.end());
  size_t k = (double)n * cutoff_ratio;
  mean_min = std::min(mean[k-1], 0.0);
  mean_max = std::max(mean[n-k], 0.0);
}

void MyViewer::updateMeanCurvature()
{
  for(MyMesh::ConstFaceIter i = mesh.faces_begin(), ie = mesh.faces_end(); i != ie; ++i)
    mesh.data(i).area = -1;

  // Compute triangle strip areas
  for(MyMesh::VertexIter i = mesh.vertices_begin(), ie = mesh.vertices_end(); i != ie; ++i) {
    mesh.data(i).area = 0;
    mesh.data(i).mean = 0;
    for(MyMesh::ConstVertexFaceIter j(mesh, i); (bool)j; ++j) {
      if(mesh.data(j).area == -1) {
        MyMesh::HalfedgeHandle h1 = mesh.halfedge_handle(j);
        MyMesh::HalfedgeHandle h2 = mesh.next_halfedge_handle(h1);
        mesh.data(j).area = (halfedgeVector(h1) % halfedgeVector(h2)).norm() / 2.0;
      }
      mesh.data(i).area += mesh.data(j).area;
    }
  }

  // Compute mean values using normal difference angles
  for(MyMesh::VertexIter i = mesh.vertices_begin(), ie = mesh.vertices_end(); i != ie; ++i) {
    for(MyMesh::ConstVertexEdgeIter j(mesh, i); (bool)j; ++j) {
      double angle;
      MyMesh::HalfedgeHandle h1 = mesh.halfedge_handle(j, 0);
      MyMesh::HalfedgeHandle h2 = mesh.halfedge_handle(j, 1);
      Vector v = halfedgeVector(h1);
      if(mesh.is_boundary(h1) || mesh.is_boundary(h2)) {
        angle = 0.0;
      } else {
        Vector n1 = mesh.normal(mesh.face_handle(h1));
        Vector n2 = mesh.normal(mesh.face_handle(h2));
        angle = n1 | n2;
        angle *= ((n1 % n2) | v) > 0.0 ? 1.0 : -1.0;
      }
      mesh.data(i).mean += angle * v.norm();
    }
    mesh.data(i).mean *= 3.0 / 4.0 / mesh.data(i).area;
  }

  updateMeanMinMax();
}

void MyViewer::meanMapColor(double d, double *color) const
{
  if(d <= mean_min) {
    color[0] = 0.0;
    color[1] = 0.0;
    color[2] = 1.0;
  } else if(d >= mean_max) {
    color[0] = 1.0;
    color[1] = 0.0;
    color[2] = 0.0;
  } else if(d < 0) {
    double alpha = d / mean_min;
    color[0] = 0.0;
    color[1] = 1.0 - alpha;
    color[2] = alpha;
  } else {
    double alpha = d / mean_max;
    color[0] = alpha;
    color[1] = 1.0 - alpha;
    color[2] = 0;
  }
}

bool MyViewer::openMesh(std::string const &filename)
{
  if(!OpenMesh::IO::read_mesh(mesh, filename) || mesh.n_vertices() == 0)
    return false;
  mesh.request_face_normals();
  mesh.update_face_normals();

  updateMeanCurvature();

  // Set camera on the model
  MyMesh::Point box_min, box_max;
  box_min = box_max = mesh.point(mesh.vertices_begin());
  for(MyMesh::ConstVertexIter i = mesh.vertices_begin(), ie = mesh.vertices_end(); i != ie; ++i) {
    box_min.minimize(mesh.point(i));
    box_max.maximize(mesh.point(i));
  }
  camera()->setSceneBoundingBox(qglviewer::Vec(box_min[0], box_min[1], box_min[2]),
                                qglviewer::Vec(box_max[0], box_max[1], box_max[2]));
  camera()->showEntireScene();

  updateGL();
  return true;
}

void MyViewer::init()
{
  glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, 1);
}

void MyViewer::draw()
{
  std::vector<double> color(3);
  for(MyMesh::ConstFaceIter i = mesh.faces_begin(), ie = mesh.faces_end(); i != ie; ++i) {
    glBegin(GL_POLYGON);
    glNormal3fv(mesh.normal(i).data());
    for(MyMesh::ConstFaceVertexIter j(mesh, i); (bool)j; ++j) {
      meanMapColor(mesh.data(j).mean, &color[0]);
      glColor3dv(&color[0]);
      glVertex3fv(mesh.point(j).data());
    }
    glEnd();
  }
}
