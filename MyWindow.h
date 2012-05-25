#pragma once

#include <QMainWindow>

#include "MyViewer.h"

class QApplication;
class QProgressBar;

class MyWindow : public QMainWindow
{
  Q_OBJECT

public:
  MyWindow(QApplication *parent);
  ~MyWindow();

private slots:
  void open();
  void setCutoff();
  void setRange();
  void startComputation(QString message);
  void midComputation(int percent);
  void endComputation();

private:
  QApplication *parent;
  MyViewer *viewer;
  QProgressBar *progress;
};
