#pragma once

#include <string>

#include <QGLViewer/qglviewer.h>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>

using qglviewer::Vec;

class MyViewer : public QGLViewer
{
  Q_OBJECT

public:
  MyViewer(QWidget *parent);
  virtual ~MyViewer();

  inline double getCutoffRatio() const;
  inline void setCutoffRatio(double ratio);
  inline double getMeanMin() const;
  inline void setMeanMin(double min);
  inline double getMeanMax() const;
  inline void setMeanMax(double max);
  bool openMesh(std::string const &filename);

signals:
  void startComputation(QString message);
  void midComputation(int percent);
  void endComputation();

protected:
  virtual void init();
  virtual void draw();
  virtual void drawWithNames();
  virtual void postSelection(const QPoint &p);
  virtual void keyPressEvent(QKeyEvent *e);
  virtual void mouseMoveEvent(QMouseEvent *e);
  virtual QString helpString() const;

private:
  struct MyTraits : public OpenMesh::DefaultTraits {
    using Point  = OpenMesh::Vec3d; // the default would be Vec3f
    using Normal = OpenMesh::Vec3d;
    VertexTraits {
      double area;              // total area of the surrounding triangles
      double mean;              // approximated mean curvature
    };
    FaceTraits {
      double area;
    };
  };
  using MyMesh = OpenMesh::TriMesh_ArrayKernelT<MyTraits>;
  using Vector = OpenMesh::VectorT<double,3>;

  void updateMeanMinMax();
  void updateMeanCurvature(bool update_min_max = true);
  void meanMapColor(double d, double *color) const;
  void fairMesh();
  void drawAxes() const;
  void drawAxesWithNames() const;
  Vec intersectLines(Vec const &ap, Vec const &ad, Vec const &bp, Vec const &bd) const;
  Vector computeNormal(MyMesh::VertexHandle v) const;
  void localSystem(Vector const &normal, Vector &u, Vector &v);
  double voronoiWeight(MyMesh::HalfedgeHandle in_he);

  MyMesh mesh;
  double mean_min, mean_max;
  double cutoff_ratio;
  bool show_solid, show_wireframe;
  enum { COLOR_PLAIN, COLOR_MEAN, COLOR_ISOPHOTES } coloring;
  GLuint isophote_texture;
  MyMesh::ConstVertexIter selected;
  struct ModificationAxes {
    bool shown;
    float size;
    int selected_axis;
    Vec position, grabbed_pos, original_pos;
  } axes;
};

#include "MyViewer.hpp"
