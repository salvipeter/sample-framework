# -*- mode: Makefile -*-

TARGET = sample-framework
CONFIG += c++14 qt opengl debug
QT += gui widgets opengl xml

HEADERS = MyWindow.h MyViewer.h MyViewer.hpp
SOURCES = MyWindow.cpp MyViewer.cpp main.cpp

unix:INCLUDEPATH += /usr/include/eigen3
unix:LIBS *= -lQGLViewer-qt5 -L/usr/lib/OpenMesh -lOpenMeshCore -lGL -lGLU

# WIN32 instructions only

# Replace this variable to the install path of the OpenMesh lib install
OPENMESH_INSTALL_PATH = 'C:\\Program Files\OpenMesh'
# If your OpenMesh source is separate from the lib install replace this variable
OPENMESH_SRC_INSTALL_PATH = 'C:\\Program Files\OpenMesh'
# Replace this variable to the install path of libQGLViewer
LIBQGLVIEWER_INSTALL_PATH = 'C:\\Progam Files\libQGLViewer'

win32:
{
    Release:LIBS += -lOpenMeshCore -lQGLViewer2
    else:Debug:LIBS += -lOpenMeshCored -lQGLViewerd2
    LIBS += -lOpenGL32 -lGLU32 
    LIBS += -L'$$OPENMESH_INSTALL_PATH\lib' -L'$$LIBQGLVIEWER_INSTALL_PATH\QGLViewer'
    INCLUDEPATH += '$$OPENMESH_SRC_INSTALL_PATH\src' '$$LIBQGLVIEWER_INSTALL_PATH'
    DEFINES += NOMINMAX 
    DEFINES += _USE_MATH_DEFINES
}


# Optional
# DEFINES += BETTER_MEAN_CURVATURE USE_JET_NORMALS
# LIBS += -lCGAL # this library will be header-only from version 5

RESOURCES = sample-framework.qrc
