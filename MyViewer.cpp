#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <vector>

#include <QtGui/QKeyEvent>

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Tools/Smoother/JacobiLaplaceSmootherT.hh>

#define BETTER_MEAN_CURVATURE
#ifdef BETTER_MEAN_CURVATURE
#include "Eigen/Eigenvalues"
#include "Eigen/Geometry"
#include "Eigen/LU"
#include "Eigen/SVD"
#endif

#include "MyViewer.h"

#ifdef _WIN32
#define GL_CLAMP_TO_EDGE 0x812F
#define GL_BGRA 0x80E1
#endif

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
    mean.push_back(mesh.data(*i).mean);

  std::sort(mean.begin(), mean.end());
  size_t k = (double)n * cutoff_ratio;
  mean_min = std::min(mean[k-1], 0.0);
  mean_max = std::max(mean[n-k], 0.0);
}

void MyViewer::localSystem(MyViewer::Vector const &normal,
                           MyViewer::Vector &u, MyViewer::Vector &v) {
  // Generates an orthogonal (u,v) coordinate system in the plane defined by `normal`.
  int maxi = 0, nexti = 1;
  double max = fabs(normal[0]), next = fabs(normal[1]);
  if (max < next) {
    std::swap(max, next);
    std::swap(maxi, nexti);
  }
  if (fabs(normal[2]) > max) {
    nexti = maxi;
    maxi = 2;
  } else if (fabs(normal[2]) > next)
    nexti = 2;

  u.vectorize(0.0);
  u[nexti] = -normal[maxi];
  u[maxi] = normal[nexti];
  u /= u.norm();
  v = cross(normal, u);
}

double MyViewer::voronoiWeight(MyViewer::MyMesh::HalfedgeHandle in_he) {
  // Returns the area of the triangle bounded by in_he that is closest
  // to the vertex pointed to by in_he.

  // If we are at point `A` with angle `alpha` and opposite edge length `a`,
  // the circumradius `r` is
  //   a / (2 * sin(alpha)),
  // and the requested area is
  //   0.5 * (area(b, r, r) + area(c, r, r)),
  // where area(a, b, c) is the area of a triangle with edge length a, b, c.

  double alpha = mesh.calc_sector_angle(in_he);
  double a = mesh.calc_edge_vector(mesh.prev_halfedge_handle(in_he)).norm();
  double b = mesh.calc_edge_vector(in_he).norm();
  double c = mesh.calc_edge_vector(mesh.next_halfedge_handle(in_he)).norm();
  double r = a / (2 * sin(alpha));
  auto area = [](double a, double b) { // Isosceles triangle
    return 0.5 * std::pow(a, 2) * sqrt(std::pow(b / a, 2) - 0.25);
  };
  return 0.5 * (area(b, r) + area(c, r));
}

#ifndef BETTER_MEAN_CURVATURE
void MyViewer::updateMeanCurvature(bool update_min_max)
{
  for(MyMesh::ConstFaceIter i = mesh.faces_begin(), ie = mesh.faces_end(); i != ie; ++i) {
    MyMesh::HalfedgeHandle h1 = mesh.halfedge_handle(*i);
    MyMesh::HalfedgeHandle h2 = mesh.next_halfedge_handle(h1);
    mesh.data(*i).area = (halfedgeVector(h1) % halfedgeVector(h2)).norm() / 2.0;
  }

  // Compute triangle strip areas
  for(MyMesh::VertexIter i = mesh.vertices_begin(), ie = mesh.vertices_end(); i != ie; ++i) {
    mesh.data(*i).area = 0;
    mesh.data(*i).mean = 0;
    for(MyMesh::ConstVertexFaceIter j(mesh, *i); j.is_valid(); ++j)
      mesh.data(*i).area += mesh.data(*j).area;
  }

  // Compute mean values using normal difference angles
  for(MyMesh::VertexIter i = mesh.vertices_begin(), ie = mesh.vertices_end(); i != ie; ++i) {
    for(MyMesh::ConstVertexEdgeIter j(mesh, *i); j.is_valid(); ++j) {
      double angle;
      MyMesh::HalfedgeHandle h1 = mesh.halfedge_handle(*j, 0);
      MyMesh::HalfedgeHandle h2 = mesh.halfedge_handle(*j, 1);
      Vector v = halfedgeVector(h1);
      if(mesh.is_boundary(h1) || mesh.is_boundary(h2))
        angle = 0.0;
      else {
        Vector n1 = mesh.normal(mesh.face_handle(h1));
        Vector n2 = mesh.normal(mesh.face_handle(h2));
        angle = acos(std::min(std::max(n1 | n2, -1.0f), 1.0f));
        angle *= ((n1 % n2) | v) >= 0.0 ? 1.0 : -1.0;
      }
      mesh.data(*i).mean += angle * v.norm();
    }
    mesh.data(*i).mean *= 3.0 / 4.0 / mesh.data(*i).area;
  }

  if(update_min_max)
    updateMeanMinMax();
}
#else // BETTER_MEAN_CURVATURE
void MyViewer::updateMeanCurvature(bool update_min_max)
{
  // As in the paper
  //   S. Rusinkiewicz, Estimating curvatures and their derivatives on triangle meshes.
  //     3D Data Processing, Visualization and Transmission, IEEE, 2004.

  // (e,f,g): 2nd principal form
  // w: accumulated weight
  std::map<MyMesh::VertexHandle, Vector> efgp;
  std::map<MyMesh::VertexHandle, double> wp;

  // Initial setup
  for (MyMesh::ConstVertexIter i = mesh.vertices_begin(), ie = mesh.vertices_end(); i != ie; ++i) {
    efgp[*i].vectorize(0.0);
    wp[*i] = 0.0;
  }

  for(MyMesh::ConstFaceIter i = mesh.faces_begin(), ie = mesh.faces_end(); i != ie; ++i) {
    // Setup local edges, vertices and normals
    auto h0 = mesh.halfedge_handle(*i);
    auto h1 = mesh.next_halfedge_handle(h0);
    auto h2 = mesh.next_halfedge_handle(h1);
    Vector e0, e1, e2;
    mesh.calc_edge_vector(h0, e0);
    mesh.calc_edge_vector(h1, e1);
    mesh.calc_edge_vector(h2, e2);
    auto v0 = mesh.to_vertex_handle(h1);
    auto v1 = mesh.to_vertex_handle(h2);
    auto v2 = mesh.to_vertex_handle(h0);
    auto n0 = mesh.normal(v0), n1 = mesh.normal(v1), n2 = mesh.normal(v2);

    Vector n = mesh.normal(*i), u, v;
    localSystem(n, u, v);

    // Solve a LSQ equation for (e,f,g) of the face
    Eigen::MatrixXd A(6, 3);
    Eigen::VectorXd b(6);
    A(0,0) = e0 | u; A(0,1) = e0 | v; A(0,2) = 0.0;    b(0) = (n2 - n1) | u;
    A(1,0) = 0.0;    A(1,1) = e0 | u; A(1,2) = e0 | v; b(1) = (n2 - n1) | v;
    A(2,0) = e1 | u; A(2,1) = e1 | v; A(2,2) = 0.0;    b(2) = (n0 - n2) | u;
    A(3,0) = 0.0;    A(3,1) = e1 | u; A(3,2) = e1 | v; b(3) = (n0 - n2) | v;
    A(4,0) = e2 | u; A(4,1) = e2 | v; A(4,2) = 0.0;    b(4) = (n1 - n0) | u;
    A(5,0) = 0.0;    A(5,1) = e2 | u; A(5,2) = e2 | v; b(5) = (n1 - n0) | v;
    Eigen::Vector3d x = A.fullPivLu().solve(b);

    if (!(A*x).isApprox(b)) {
      // Singular case
      Eigen::JacobiSVD<Eigen::MatrixXd> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
      x = svd.solve(b);
    }

    Eigen::Matrix2d F;          // Fundamental matrix for the face
    F << x(0), x(1),
         x(1), x(2);

    auto h = h0;
    do {
      auto p = mesh.to_vertex_handle(h);

      // Rotate the (up,vp) local coordinate system to be coplanar with that of the face
      Vector up, vp;
      auto np = mesh.normal(p);
      localSystem(np, up, vp);
      auto axis = cross(np, n); axis /= axis.norm();
      double angle = acos(std::min(std::max(n | np, -1.0f), 1.0f));
      auto rotation = Eigen::AngleAxisd(angle, Eigen::Vector3d(axis[0], axis[1], axis[2]));
      Eigen::Vector3d up1(up[0], up[1], up[2]), vp1(vp[0], vp[1], vp[2]);
      up1 = rotation * up1; vp1 = rotation * vp1;
      up = Vector(up1(0), up1(1), up1(2));
      vp = Vector(vp1(0), vp1(1), vp1(2));

      // Compute the vertex-local (e,f,g)
      double e, f, g;
      Eigen::Vector2d upf, vpf;
      upf << (up | u), (up | v);
      vpf << (vp | u), (vp | v);
      e = upf.transpose() * F * upf;
      f = upf.transpose() * F * vpf;
      g = vpf.transpose() * F * vpf;

      // Accumulate the results with Voronoi weights
      double w = voronoiWeight(h);
      efgp[p] += Vector(e, f, g) * w;
      wp[p] += w;
      
      // Next halfedge
      h = mesh.next_halfedge_handle(h);
    } while (h != h0);
  }

  // Compute the principal curvatures
  for(MyMesh::VertexIter i = mesh.vertices_begin(), ie = mesh.vertices_end(); i != ie; ++i) {
    auto &efg = efgp[*i];
    efg /= wp[*i];
    Eigen::Matrix2d F;
    F << efg[0], efg[1],
         efg[1], efg[2];
    auto k = F.eigenvalues();   // always real, because F is a symmetric real matrix
    mesh.data(*i).mean = (k(0).real() + k(1).real()) / 2.0;
  }

  if(update_min_max)
    updateMeanMinMax();
}
#endif

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

MyViewer::MyMesh::Normal MyViewer::computeNormal(MyViewer::MyMesh::VertexHandle v) const {
  // Weights according to:
  //   N. Max, Weights for computing vertex normals from facet normals.
  //     Journal of Graphics Tools, Vol. 4(2), 1999.

  // Based on calc_vertex_normal{,_correct} in OpenMesh::PolyMeshT
  MyMesh::Normal n;

  n.vectorize(0.0);
  MyMesh::ConstVertexIHalfedgeIter cvih_it = mesh.cvih_iter(v);
  if (!cvih_it.is_valid())
    return n;

  MyMesh::Normal in_he_vec;
  mesh.calc_edge_vector(*cvih_it, in_he_vec);
  for (; cvih_it.is_valid(); ++cvih_it) {
    if (mesh.is_boundary(*cvih_it))
      continue;
    MyMesh::HalfedgeHandle out_heh(mesh.next_halfedge_handle(*cvih_it));
    MyMesh::Normal out_he_vec;
    mesh.calc_edge_vector(out_heh, out_he_vec);
    double w = in_he_vec.sqrnorm() * out_he_vec.sqrnorm();
    n += cross(in_he_vec, out_he_vec) / (w == 0.0 ? 1.0 : w);
    in_he_vec = out_he_vec;
    in_he_vec *= -1;
  }

  double len = n.length();
  if (len != 0.0)
    n *= 1 / len;
  return n;
}

bool MyViewer::openMesh(std::string const &filename)
{
  if(!OpenMesh::IO::read_mesh(mesh, filename) || mesh.n_vertices() == 0)
    return false;
  mesh.request_face_normals(); mesh.request_vertex_normals();
  mesh.update_face_normals();  //mesh.update_vertex_normals();

  // update vertex normals by hand
  for (MyMesh::VertexIter i = mesh.vertices_begin(), ie = mesh.vertices_end(); i != ie; ++i)
    mesh.set_normal(*i, computeNormal(*i));

  updateMeanCurvature();

  // Set camera on the model
  MyMesh::Point box_min, box_max;
  box_min = box_max = mesh.point(*mesh.vertices_begin());
  for(MyMesh::ConstVertexIter i = mesh.vertices_begin(), ie = mesh.vertices_end(); i != ie; ++i) {
    box_min.minimize(mesh.point(*i));
    box_max.maximize(mesh.point(*i));
  }
  camera()->setSceneBoundingBox(Vec(box_min[0], box_min[1], box_min[2]),
                                Vec(box_max[0], box_max[1], box_max[2]));
  camera()->showEntireScene();

  setSelectedName(-1);
  axes.shown = false;

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
               GL_UNSIGNED_BYTE, img.convertToFormat(QImage::Format_ARGB32).bits());
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
      for(MyMesh::ConstFaceVertexIter j(mesh, *i); j.is_valid(); ++j) {
        if(coloring == COLOR_MEAN) {
          meanMapColor(mesh.data(*j).mean, &color[0]);
          glColor3dv(&color[0]);
        }
        glNormal3fv(mesh.normal(*j).data());
        glVertex3fv(mesh.point(*j).data());
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
      for(MyMesh::ConstFaceVertexIter j(mesh, *i); j.is_valid(); ++j)
        glVertex3fv(mesh.point(*j).data());
      glEnd();
    }
    glEnable(GL_LIGHTING);
  }

  if(axes.shown)
    drawAxes();
}

void MyViewer::drawAxes() const
{
  Vec const p(axes.position[0], axes.position[1], axes.position[2]);
  glColor3f(1.0, 0.0, 0.0);
  drawArrow(p, p + Vec(axes.size, 0.0, 0.0), axes.size / 50.0);
  glColor3f(0.0, 1.0, 0.0);
  drawArrow(p, p + Vec(0.0, axes.size, 0.0), axes.size / 50.0);
  glColor3f(0.0, 0.0, 1.0);
  drawArrow(p, p + Vec(0.0, 0.0, axes.size), axes.size / 50.0);
  glEnd();
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
      glRasterPos3fv(mesh.point(*i).data());
      glPopName();
    }
  }
}

void MyViewer::drawAxesWithNames() const
{
  Vec const p(axes.position[0], axes.position[1], axes.position[2]);
  glPushName(0);
  drawArrow(p, p + Vec(axes.size, 0.0, 0.0), axes.size / 50.0);
  glPopName();
  glPushName(1);
  drawArrow(p, p + Vec(0.0, axes.size, 0.0), axes.size / 50.0);
  glPopName();
  glPushName(2);
  drawArrow(p, p + Vec(0.0, 0.0, axes.size), axes.size / 50.0);
  glPopName();
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
    axes.position[0] = mesh.point(*i).data()[0];
    axes.position[1] = mesh.point(*i).data()[1];
    axes.position[2] = mesh.point(*i).data()[2];
    Vec const pos(axes.position[0], axes.position[1], axes.position[2]);
    double const depth = camera()->projectedCoordinatesOf(pos)[2];
    Vec const q1 = camera()->unprojectedCoordinatesOf(Vec(0.0, 0.0, depth));
    Vec const q2 = camera()->unprojectedCoordinatesOf(Vec(width(), height(), depth));
    axes.size = (q1-q2).norm() / 10.0;
    axes.shown = true;
    axes.selected_axis = -1;
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
  // always returns a point on the (ap, ad) line
  double a = ad * ad, b = ad * bd, c = bd * bd;
  double d = ad * (ap - bp), e = bd * (ap - bp);
  if(a * c - b * b < 1.0e-7)
    return ap;
  double s = (b * e - c * d) / (a * c - b * b);
  return ap + s * ad;
}

void MyViewer::mouseMoveEvent(QMouseEvent *e)
{
  if(axes.shown && axes.selected_axis >= 0 &&
     e->modifiers() & Qt::ShiftModifier && e->buttons() & Qt::LeftButton) {
    Vec axis = Vec(axes.selected_axis == 0, axes.selected_axis == 1, axes.selected_axis == 2);
    Vec from, dir;
    camera()->convertClickToLine(e->pos(), from, dir);
    Vec p = intersectLines(axes.grabbed_pos, axis, from, dir);
    float d = (p - axes.grabbed_pos) * axis;
    axes.position[axes.selected_axis] = axes.original_pos[axes.selected_axis] + d;
    mesh.set_point(*selected, MyMesh::Point(axes.position[0],
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
