// -*- mode: c++ -*-
#pragma once

#include <QtWidgets/QMainWindow>

#include "MyViewer.h"

class QApplication;
class QProgressBar;

class MyWindow : public QMainWindow {
  Q_OBJECT

public:
  explicit MyWindow(QApplication *parent);
  ~MyWindow();

private slots:
  void open();
  void save();
  void setCutoff();
  void setRange();
  void setSlicing();
  void startComputation(QString message);
  void midComputation(int percent);
  void endComputation();

private:
  QApplication *parent;
  MyViewer *viewer;
  QProgressBar *progress;
  QString last_directory;
};
