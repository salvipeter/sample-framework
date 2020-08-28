# -*- mode: Makefile -*-

TARGET = sample-framework
CONFIG += c++14 qt opengl debug
QT += gui widgets opengl xml
QMAKE_CXX = clang++

HEADERS = MyWindow.h MyViewer.h MyViewer.hpp
SOURCES = MyWindow.cpp MyViewer.cpp main.cpp

INCLUDEPATH += /usr/include/eigen3
LIBS *= -lQGLViewer-qt5 -L/usr/lib/OpenMesh -lOpenMeshCore -lGL -lGLU

# Optional
# DEFINES += BETTER_MEAN_CURVATURE USE_JET_FITTING
# LIBS += -lCGAL

RESOURCES = sample-framework.qrc
