# -*- mode: Makefile -*-

TARGET = onallo-alap
CONFIG *= qt opengl release
QT *= opengl xml

HEADERS = MyWindow.h MyViewer.h
SOURCES = MyWindow.cpp MyViewer.cpp main.cpp

LIBS *= -lQGLViewer -L/usr/lib/OpenMesh -lOpenMeshCore
