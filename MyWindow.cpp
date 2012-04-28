#include <QFileDialog>
#include <QMenuBar>
#include <QMessageBox>
#include <QStatusBar>

#include "MyWindow.h"

MyWindow::MyWindow() : QMainWindow()
{
  setWindowTitle(tr("Sample 3D Framework"));
  setStatusBar(new QStatusBar);

  viewer = new MyViewer(this);
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

  QMenu *fileMenu = menuBar()->addMenu(tr("&File"));
  fileMenu->addAction(openAction);
  fileMenu->addAction(quitAction);
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
