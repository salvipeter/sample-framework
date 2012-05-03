#pragma once

#include <QMainWindow>

#include "MyViewer.h"

class MyWindow : public QMainWindow
{
  Q_OBJECT

public:
  MyWindow();
  ~MyWindow();

private slots:
  void open();
  void setCutoff();
  void setRange();

private:
  MyViewer *viewer;
};
