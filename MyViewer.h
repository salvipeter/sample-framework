#pragma once

#include <string>

#include <QGLViewer/qglviewer.h>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>

class MyViewer : public QGLViewer
{
  Q_OBJECT

public:
  MyViewer(QWidget *parent);
  virtual ~MyViewer();
  bool openMesh(std::string const &filename);

protected:
  virtual void init();
  virtual void draw();

private:
  typedef OpenMesh::PolyMesh_ArrayKernelT<> MyMesh;
  MyMesh mesh;
};
