#include <OpenMesh/Core/IO/MeshIO.hh>

#include "MyViewer.h"

MyViewer::MyViewer(QWidget *parent) : QGLViewer(parent)
{
}

MyViewer::~MyViewer()
{
}

bool MyViewer::openMesh(std::string const &filename)
{
  if(!OpenMesh::IO::read_mesh(mesh, filename) || mesh.n_vertices() == 0)
    return false;
  mesh.request_face_normals();
  mesh.update_face_normals();

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
  for(MyMesh::ConstFaceIter i = mesh.faces_begin(), ie = mesh.faces_end(); i != ie; ++i) {
    glBegin(GL_POLYGON);
    glNormal3fv(mesh.normal(i).data());
    for(MyMesh::ConstFaceVertexIter j(mesh, i); (bool)j; ++j)
      glVertex3fv(mesh.point(j).data());
    glEnd();
  }
}
