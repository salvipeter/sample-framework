# -*- mode: Makefile -*-

TARGET = sample-framework
CONFIG += c++14 qt opengl debug
QT += gui widgets opengl xml

HEADERS = MyWindow.h MyViewer.h MyViewer.hpp
SOURCES = MyWindow.cpp MyViewer.cpp main.cpp jet-wrapper.cpp

QMAKE_CXXFLAGS += -O3

unix:INCLUDEPATH += /usr/include/eigen3
unix:LIBS *= -lQGLViewer-qt5 -lOpenMeshCore -lGL -lGLU

RESOURCES = sample-framework.qrc

# Optional
# DEFINES += BETTER_MEAN_CURVATURE
# DEFINES += USE_JET_FITTING
# LIBS += -lCGAL # this library will be header-only from version 5

###########################
# WIN32 instructions only #
###########################

win32 {
    # Replace this variable to the install path of the OpenMesh lib install
    OPENMESH_INSTALL_PATH = 'C:\Program Files\OpenMesh'

    # Replace this variable to the install path of libQGLViewer
    LIBQGLVIEWER_INSTALL_PATH = 'C:\Program Files\libQGLViewer'

    # If your OpenMesh source is separate from the lib install replace this variable
    OPENMESH_SRC_INSTALL_PATH = $$OPENMESH_INSTALL_PATH\include\

    DEFINES += NOMINMAX _USE_MATH_DEFINES

    INCLUDEPATH += '$$OPENMESH_SRC_INSTALL_PATH' '$$LIBQGLVIEWER_INSTALL_PATH'

    LIBS += -lOpenGL32 -lGLU32
    LIBS += -L'$$OPENMESH_INSTALL_PATH\lib' -L'$$LIBQGLVIEWER_INSTALL_PATH\QGLViewer'
    Release:LIBS += -lOpenMeshCore -lQGLViewer2
    Debug:LIBS += -lOpenMeshCored -lQGLViewerd2
}
