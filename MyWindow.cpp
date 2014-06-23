#include <limits>

#include <QtWidgets>

#include "MyWindow.h"

MyWindow::MyWindow(QApplication *parent) : QMainWindow(), parent(parent)
{
  setWindowTitle(tr("Sample 3D Framework"));
  setStatusBar(new QStatusBar);
  progress = new QProgressBar;
  progress->setMinimum(0); progress->setMaximum(100);
  progress->hide();
  statusBar()->addPermanentWidget(progress);

  viewer = new MyViewer(this);
  connect(viewer, SIGNAL(startComputation(QString)), this, SLOT(startComputation(QString)));
  connect(viewer, SIGNAL(midComputation(int)), this, SLOT(midComputation(int)));
  connect(viewer, SIGNAL(endComputation()), this, SLOT(endComputation()));
  setCentralWidget(viewer);

  /////////////////////////
  // Setup actions/menus //
  /////////////////////////

  QAction *openAction = new QAction(tr("&Open"), this);
  openAction->setShortcut(tr("Ctrl+O"));
  openAction->setStatusTip(tr("Load a model from a file"));
  connect(openAction, SIGNAL(triggered()), this, SLOT(open()));

  QAction *quitAction = new QAction(tr("&Quit"), this);
  quitAction->setShortcut(tr("Ctrl+Q"));
  quitAction->setStatusTip(tr("Quit the program"));
  connect(quitAction, SIGNAL(triggered()), this, SLOT(close()));

  QAction *cutoffAction = new QAction(tr("Set &cutoff ratio"), this);
  cutoffAction->setStatusTip(tr("Set mean map cutoff ratio"));
  connect(cutoffAction, SIGNAL(triggered()), this, SLOT(setCutoff()));

  QAction *rangeAction = new QAction(tr("Set &range"), this);
  rangeAction->setStatusTip(tr("Set mean map range"));
  connect(rangeAction, SIGNAL(triggered()), this, SLOT(setRange()));

  QMenu *fileMenu = menuBar()->addMenu(tr("&File"));
  fileMenu->addAction(openAction);
  fileMenu->addAction(quitAction);

  QMenu *visMenu = menuBar()->addMenu(tr("&Visualization"));
  visMenu->addAction(cutoffAction);
  visMenu->addAction(rangeAction);
}

MyWindow::~MyWindow()
{
}

void MyWindow::open()
{
  QString fileName =
    QFileDialog::getOpenFileName(this, tr("Open File"), ".",
                                 tr("PLY Mesh (*.ply);;STL Mesh (*.stl);;All files (*.*)"));
  if(!fileName.isEmpty()) 
    if(!viewer->openMesh(fileName.toUtf8().data()))
      QMessageBox::warning(this, tr("Cannot open file"),
                           tr("Could not open file: ") + fileName + ".");
}

void MyWindow::setCutoff()
{
  QDialog *dlg = new QDialog(this);
  QHBoxLayout *hb1 = new QHBoxLayout, *hb2 = new QHBoxLayout;
  QVBoxLayout *vb = new QVBoxLayout;
  QLabel *text = new QLabel(tr("Cutoff ratio:"));
  QDoubleSpinBox *sb = new QDoubleSpinBox;
  QPushButton *cancel = new QPushButton(tr("Cancel"));
  QPushButton *ok = new QPushButton(tr("Ok"));

  sb->setDecimals(3);
  sb->setRange(0.001, 0.5);
  sb->setSingleStep(0.01);
  sb->setValue(viewer->getCutoffRatio());
  connect(cancel, SIGNAL(pressed()), dlg, SLOT(reject()));
  connect(ok, SIGNAL(pressed()), dlg, SLOT(accept()));
  ok->setDefault(true);

  hb1->addWidget(text);
  hb1->addWidget(sb);
  hb2->addWidget(cancel);
  hb2->addWidget(ok);
  vb->addLayout(hb1);
  vb->addLayout(hb2);

  dlg->setWindowTitle(tr("Set ratio"));
  dlg->setLayout(vb);

  if(dlg->exec() == QDialog::Accepted) {
    viewer->setCutoffRatio(sb->value());
    viewer->updateGL();
  }
}

void MyWindow::setRange()
{
  QDialog *dlg = new QDialog(this);
  QGridLayout *grid = new QGridLayout;
  QLabel *text1 = new QLabel(tr("Min:")), *text2 = new QLabel(tr("Max:"));
  QDoubleSpinBox *sb1 = new QDoubleSpinBox, *sb2 = new QDoubleSpinBox;
  QPushButton *cancel = new QPushButton(tr("Cancel"));
  QPushButton *ok = new QPushButton(tr("Ok"));

  double max = std::numeric_limits<double>::max();
  sb1->setDecimals(5);                 sb2->setDecimals(5);
  sb1->setRange(-max, 0.0);            sb2->setRange(0.0, max);
  sb1->setSingleStep(0.01);            sb2->setSingleStep(0.01);
  sb1->setValue(viewer->getMeanMin()); sb2->setValue(viewer->getMeanMax());
  connect(cancel, SIGNAL(pressed()), dlg, SLOT(reject()));
  connect(ok, SIGNAL(pressed()), dlg, SLOT(accept()));
  ok->setDefault(true);

  grid->addWidget( text1, 1, 1, Qt::AlignRight);
  grid->addWidget(   sb1, 1, 2);
  grid->addWidget( text2, 2, 1, Qt::AlignRight);
  grid->addWidget(   sb2, 2, 2);
  grid->addWidget(cancel, 3, 1);
  grid->addWidget(    ok, 3, 2);
  dlg->setWindowTitle(tr("Set range"));
  dlg->setLayout(grid);

  if(dlg->exec() == QDialog::Accepted) {
    viewer->setMeanMin(sb1->value());
    viewer->setMeanMax(sb2->value());
    viewer->updateGL();
  }
}

void MyWindow::startComputation(QString message)
{
  statusBar()->showMessage(message);
  progress->setValue(0);
  progress->show();
  parent->processEvents(QEventLoop::ExcludeUserInputEvents);
}

void MyWindow::midComputation(int percent)
{
  progress->setValue(percent);
  parent->processEvents(QEventLoop::ExcludeUserInputEvents);
}

void MyWindow::endComputation()
{
  progress->hide();
  statusBar()->clearMessage();
}
