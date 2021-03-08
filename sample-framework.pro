# -*- mode: Makefile -*-

TARGET = sample-framework
CONFIG += c++14 qt opengl debug
QT += gui widgets opengl xml

HEADERS = MyWindow.h MyViewer.h MyViewer.hpp
SOURCES = MyWindow.cpp MyViewer.cpp main.cpp

INCLUDEPATH += /usr/include/eigen3
LIBS *= -lQGLViewer-qt5 -L/usr/lib/OpenMesh -lOpenMeshCore -lGL -lGLU

# Optional
# DEFINES += BETTER_MEAN_CURVATURE USE_JET_NORMALS
# LIBS += -lCGAL # this library will be header-only from version 5

RESOURCES = sample-framework.qrc
