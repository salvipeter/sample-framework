# -*- mode: Makefile -*-

TARGET = sample-framework
CONFIG *= qt opengl release
QT *= opengl xml

HEADERS = MyWindow.h MyViewer.h MyViewer.hpp
SOURCES = MyWindow.cpp MyViewer.cpp main.cpp

LIBS *= -lQGLViewer -L/usr/lib/OpenMesh -lOpenMeshCore

RESOURCES = sample-framework.qrc
