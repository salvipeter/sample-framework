double MyViewer::getCutoffRatio() const
{
  return cutoff_ratio;
}

void MyViewer::setCutoffRatio(double ratio)
{
  cutoff_ratio = ratio;
  updateMeanMinMax();
}

double MyViewer::getMeanMin() const
{
  return mean_min;
}

void MyViewer::setMeanMin(double min)
{
  mean_min = min;
}

double MyViewer::getMeanMax() const
{
  return mean_max;
}

void MyViewer::setMeanMax(double max)
{
  mean_max = max;
}

MyViewer::Vector MyViewer::halfedgeVector(MyMesh::HalfedgeHandle const &h) const
{
  return mesh.point(mesh.to_vertex_handle(h)) - mesh.point(mesh.from_vertex_handle(h));
}
