#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

#include <QKeyEvent>

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Tools/Smoother/JacobiLaplaceSmootherT.hh>

#include "MyViewer.h"

using qglviewer::Vec;

MyViewer::MyViewer(QWidget *parent) :
  QGLViewer(parent),
  mean_min(0.0), mean_max(0.0), cutoff_ratio(0.05),
  show_solid(true), show_wireframe(false), coloring(COLOR_PLAIN)
{
  setSelectRegionWidth(5);
  setSelectRegionHeight(5);
  axes.shown = false;
}

MyViewer::~MyViewer()
{
  glDeleteTextures(1, &isophote_texture);
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

void MyViewer::updateMeanCurvature(bool update_min_max)
{
  for(MyMesh::ConstFaceIter i = mesh.faces_begin(), ie = mesh.faces_end(); i != ie; ++i)
    mesh.data(i).area = -1;

  // Compute triangle strip areas
  for(MyMesh::VertexIter i = mesh.vertices_begin(), ie = mesh.vertices_end(); i != ie; ++i) {
    mesh.data(i).area = 0;
    mesh.data(i).mean = 0;
    for(MyMesh::ConstVertexFaceIter j(mesh, i); (bool)j; ++j) {
      if(mesh.data(j).area == -1) {
        MyMesh::HalfedgeHandle h1 = mesh.halfedge_handle((OpenMesh::FaceHandle const &)j);
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
      if(mesh.is_boundary(h1) || mesh.is_boundary(h2))
        angle = 0.0;
      else {
        Vector n1 = mesh.normal(mesh.face_handle(h1));
        Vector n2 = mesh.normal(mesh.face_handle(h2));
        angle = acos(std::min(std::max(n1 | n2, -1.0f), 1.0f));
        angle *= ((n1 % n2) | v) >= 0.0 ? 1.0 : -1.0;
      }
      mesh.data(i).mean += angle * v.norm();
    }
    mesh.data(i).mean *= 3.0 / 4.0 / mesh.data(i).area;
  }

  if(update_min_max)
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

void MyViewer::fairMesh()
{
  emit startComputation(tr("Fairing mesh..."));
  OpenMesh::Smoother::JacobiLaplaceSmootherT<MyMesh> smoother(mesh);
  smoother.initialize(OpenMesh::Smoother::SmootherT<MyMesh>::Normal, // or: Tangential_and_Normal
                      OpenMesh::Smoother::SmootherT<MyMesh>::C1);
  for(size_t i = 1; i <= 10; ++i) {
    smoother.smooth(10);
    emit midComputation(i * 10);
  }
  mesh.update_face_normals();
  mesh.update_vertex_normals();
  updateMeanCurvature(false);
  emit endComputation();
}

bool MyViewer::openMesh(std::string const &filename)
{
  if(!OpenMesh::IO::read_mesh(mesh, filename) || mesh.n_vertices() == 0)
    return false;
  mesh.request_face_normals(); mesh.request_vertex_normals();
  mesh.update_face_normals();  mesh.update_vertex_normals();

  updateMeanCurvature();

  // Set camera on the model
  MyMesh::Point box_min, box_max;
  box_min = box_max = mesh.point(mesh.vertices_begin());
  for(MyMesh::ConstVertexIter i = mesh.vertices_begin(), ie = mesh.vertices_end(); i != ie; ++i) {
    box_min.minimize(mesh.point(i));
    box_max.maximize(mesh.point(i));
  }
  camera()->setSceneBoundingBox(Vec(box_min[0], box_min[1], box_min[2]),
                                Vec(box_max[0], box_max[1], box_max[2]));
  camera()->showEntireScene();

  setSelectedName(-1);
  updateGL();
  return true;
}

void MyViewer::init()
{
  glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, 1);
  QImage img(":/isophotes.png");
  glGenTextures(1, &isophote_texture);
  glBindTexture(GL_TEXTURE_2D, isophote_texture);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, img.width(), img.height(), 0, GL_BGRA,
               GL_UNSIGNED_BYTE, img.convertToFormat(QImage::Format_ARGB32).constBits());
}

void MyViewer::draw()
{
  if(!show_solid && show_wireframe)
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
  else
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

  glEnable(GL_POLYGON_OFFSET_FILL);
  glPolygonOffset(1, 1);

  std::vector<double> color(3, 1.0);
  if(show_solid || show_wireframe) {
    if(coloring == COLOR_PLAIN)
      glColor3dv(&color[0]);
    else if(coloring == COLOR_ISOPHOTES) {
      glBindTexture(GL_TEXTURE_2D, isophote_texture);
      glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
      glEnable(GL_TEXTURE_2D);
      glTexGeni(GL_S, GL_TEXTURE_GEN_MODE, GL_SPHERE_MAP);
      glTexGeni(GL_T, GL_TEXTURE_GEN_MODE, GL_SPHERE_MAP);
      glEnable(GL_TEXTURE_GEN_S);
      glEnable(GL_TEXTURE_GEN_T);
    }
    for(MyMesh::ConstFaceIter i = mesh.faces_begin(), ie = mesh.faces_end(); i != ie; ++i) {
      glBegin(GL_POLYGON);
      for(MyMesh::ConstFaceVertexIter j(mesh, i); (bool)j; ++j) {
        if(coloring == COLOR_MEAN) {
          meanMapColor(mesh.data(j).mean, &color[0]);
          glColor3dv(&color[0]);
        }
        glNormal3fv(mesh.normal(j).data());
        glVertex3fv(mesh.point(j).data());
      }
      glEnd();
    }
    if(coloring == COLOR_ISOPHOTES) {
      glDisable(GL_TEXTURE_GEN_S);
      glDisable(GL_TEXTURE_GEN_T);
      glDisable(GL_TEXTURE_2D);
      glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
    }
  }

  if(show_solid && show_wireframe) {
    glPolygonMode(GL_FRONT, GL_LINE);
    glColor3d(0.0, 0.0, 0.0);
    glDisable(GL_LIGHTING);
    for(MyMesh::ConstFaceIter i = mesh.faces_begin(), ie = mesh.faces_end(); i != ie; ++i) {
      glBegin(GL_POLYGON);
      for(MyMesh::ConstFaceVertexIter j(mesh, i); (bool)j; ++j)
        glVertex3fv(mesh.point(j).data());
      glEnd();
    }
    glEnable(GL_LIGHTING);
  }

  if(axes.shown)
    drawAxes();
}

void MyViewer::drawAxes() const
{
  glDisable(GL_LIGHTING);
  glLineWidth(5.0);
  glBegin(GL_LINES);
  glColor3f(1.0, 0.0, 0.0);
  glVertex3f(axes.position[0], axes.position[1], axes.position[2]);
  glVertex3f(axes.position[0] + axes.size, axes.position[1], axes.position[2]);
  glColor3f(0.0, 1.0, 0.0);
  glVertex3f(axes.position[0], axes.position[1], axes.position[2]);
  glVertex3f(axes.position[0], axes.position[1] + axes.size, axes.position[2]);
  glColor3f(0.0, 0.0, 1.0);
  glVertex3f(axes.position[0], axes.position[1], axes.position[2]);
  glVertex3f(axes.position[0], axes.position[1], axes.position[2] + axes.size);
  glEnd();
  glLineWidth(1.0);
  glEnable(GL_LIGHTING);
}

void MyViewer::drawWithNames()
{
  if(axes.shown)
    drawAxesWithNames();
  else {
    if(!show_wireframe)
      return;

    int j = 0;
    for(MyMesh::ConstVertexIter i = mesh.vertices_begin(), ie = mesh.vertices_end(); i != ie; ++i) {
      glPushName(j++);
      glRasterPos3fv(mesh.point(i).data());
      glPopName();
    }
  }
}

void MyViewer::drawAxesWithNames() const
{
  glLineWidth(5.0);
  glPushName(0);
  glBegin(GL_LINES);
  glVertex3f(axes.position[0], axes.position[1], axes.position[2]);
  glVertex3f(axes.position[0] + axes.size, axes.position[1], axes.position[2]);
  glEnd();
  glPopName();
  glPushName(1);
  glBegin(GL_LINES);
  glVertex3f(axes.position[0], axes.position[1], axes.position[2]);
  glVertex3f(axes.position[0], axes.position[1] + axes.size, axes.position[2]);
  glEnd();
  glPopName();
  glPushName(2);
  glBegin(GL_LINES);
  glVertex3f(axes.position[0], axes.position[1], axes.position[2]);
  glVertex3f(axes.position[0], axes.position[1], axes.position[2] + axes.size);
  glEnd();
  glPopName();
  glLineWidth(1.0);
}

void MyViewer::postSelection(const QPoint &p)
{
  int sel = selectedName();
  if(sel == -1) {
    axes.shown = false;
    return;
  }

  if(axes.shown) {
    axes.selected_axis = sel;
    bool found;
    axes.grabbed_pos = camera()->pointUnderPixel(p, found);
    axes.original_pos[0] = axes.position[0];
    axes.original_pos[1] = axes.position[1];
    axes.original_pos[2] = axes.position[2];
    if(!found)
      axes.shown = false;
  } else {
    MyMesh::ConstVertexIter i = mesh.vertices_begin();
    for(int j = 0; j != sel; ++i, ++j);
    selected = i;
    axes.position[0] = mesh.point(i).data()[0];
    axes.position[1] = mesh.point(i).data()[1];
    axes.position[2] = mesh.point(i).data()[2];
    Vec const q1 = camera()->worldCoordinatesOf(Vec(0, 0, 0));
    Vec const q2 = camera()->worldCoordinatesOf(Vec(1.0, 1.0, 0));
    axes.size = std::sqrt(q1*q1+q2*q2) / 10.0;
    axes.shown = true;
  }
}

void MyViewer::keyPressEvent(QKeyEvent *e)
{
  if(e->modifiers() == Qt::NoModifier)
    switch(e->key()) {
    case Qt::Key_P:
      coloring = COLOR_PLAIN;
      updateGL();
      break;
    case Qt::Key_M:
      coloring = COLOR_MEAN;
      updateGL();
      break;
    case Qt::Key_I:
      coloring = COLOR_ISOPHOTES;
      updateGL();
      break;
    case Qt::Key_S:
      show_solid = !show_solid;
      updateGL();
      break;
    case Qt::Key_W:
      show_wireframe = !show_wireframe;
      updateGL();
      break;
    case Qt::Key_F:
      fairMesh();
      updateGL();
      break;
    default:
      QGLViewer::keyPressEvent(e);
    }
  else
    QGLViewer::keyPressEvent(e);
}

Vec MyViewer::intersectLines(Vec const &ap, Vec const &ad, Vec const &bp, Vec const &bd) const
{
  // This can be solved analitically as well (may be faster)
  float min = 2 * sceneRadius();
  Vec result;
  for(float t = 0.0, te = 2 * sceneRadius(), td = te / 100.0; t < te; t += td) {
    Vec p = bp + t * bd;
    Vec q = ap + ad * ((p - ap) * ad);
    float d = (p - q).norm();
    if(d < min) {
      min = d;
      result = q;
    }
  }
  return result;
}

void MyViewer::mouseMoveEvent(QMouseEvent *e)
{
  if(axes.shown && selectedName() != -1 &&
     e->modifiers() & Qt::ShiftModifier && e->buttons() & Qt::LeftButton) {
    Vec axis = Vec(axes.selected_axis == 0, axes.selected_axis == 1, axes.selected_axis == 2);
    Vec from, dir;
    camera()->convertClickToLine(e->pos(), from, dir);
    Vec p = intersectLines(axes.grabbed_pos, axis, from, dir);
    float d = (p - axes.grabbed_pos) * axis;
    axes.position[axes.selected_axis] = axes.original_pos[axes.selected_axis] + d;
    mesh.set_point(selected, MyMesh::Point(axes.position[0],
                                           axes.position[1],
                                           axes.position[2]));
    updateGL();
  } else
    QGLViewer::mouseMoveEvent(e);
}

QString MyViewer::helpString() const
{
  QString text("<h2>Sample Framework</h2>"
               "<p>This is a minimal framework for 3D mesh manipulation, which can be "
               "extended and used as a base for various projects, for example "
               "prototypes for fairing algorithms, or even displaying/modifying "
               "parametric surfaces, etc.</p>"
               "<p>The following hotkeys are available:</p>"
               "<ul>"
               "<li>&nbsp;P: Set plain map (no coloring)</li>"
               "<li>&nbsp;M: Set mean curvature map</li>"
               "<li>&nbsp;I: Set isophote line map</li>"
               "<li>&nbsp;S: Toggle solid (filled polygon) visualization</li>"
               "<li>&nbsp;W: Toggle wireframe visualization</li>"
               "<li>&nbsp;F: Fair mesh</li>"
               "</ul>"
               "<p>There is also a simple selection and movement interface, enabled "
               "only when the wireframe is displayed: a mesh vertex can be selected "
               "by shift-clicking, and it can be moved by shift-dragging one of the "
               "displayed axes.</p>"
               "<p>This is evidently of little practical use; it serves "
               "only to demonstrate the selection and movement process.</p>"
               "<p>Note that libQGLViewer is furnished with a lot of useful features, "
               "such as storing/loading view positions, or saving screenshots. "
               "OpenMesh also has a nice collection of tools for mesh manipulation: "
               "decimation, subdivision, smoothing, etc. These can provide "
               "good comparisons to the methods you implement.</p>"
               "<p>This software can be used as a sample GUI base for handling "
               "parametric or procedural surfaces, as well. The power of "
               "Qt and libQGLViewer makes it easy to set up a prototype application. "
               "Feel free to modify and explore!</p>"
               "<p align=\"right\">Peter Salvi</p>");
  return text;
}
