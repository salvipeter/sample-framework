# Usage

Pressing `h` displays a help window that contains all the necessary information.
If you want to use this as a base for your own project, examine the source code
and read the documentation of the various libraries (see "Dependencies" below).

# INSTALLATION

Dependencies:

- Qt5 (http://qt-project.org/)
- libQGLViewer (http://www.libqglviewer.com/)
- OpenMesh (http://www.openmesh.org/)

## Linux

    qmake && make

Assumes that:

- qmake (Qt 5.x) is in the path
- QGLViewer include files are in `<standard include directory>/QGLViewer/`
- OpenMesh  include files are in `<standard include directory>/OpenMesh/`
- QGLViewer library file  is  in `<standard library directory>/`
- OpenMesh  library files are in `/usr/lib/OpenMesh/`

If any of the above is not satisfied, edit sample-framework.pro accordingly.

## Windows / Visual Studio

1. Install the Qt SDK, which should integrate itself into Visual Studio.

1. Download the source package for libQGLViewer, and put it somewhere,
   e.g. at c:\Program Files\libQGLViewer. Open Visual Studio,
   in the Qt menu select "Open Qt project file (*.pro)",
   and open the project file in the QGLViewer subdirectory
   (in my case c:\Program Files\libQGLViewer\libQGLViewer.pro).
   Compile a release version of the library.

1. Now open sample-framework.pro, and modify the following line:

        LIBQGLVIEWER_INSTALL_PATH = 'C:\Program Files\libQGLViewer'
   (using the correct path on your system, of course).

1. Download the source package for OpenMesh, and put it somewhere,
   e.g. at c:\Program Files\OpenMesh. Open the solution file
   in Visual Studio and build a release version of the core library.
   Then open sample-framework.pro, and modify the following line:

        OPENMESH_INSTALL_PATH = 'C:\Program Files\OpenMesh'
   (using the correct path on your system, of course).

1. Open Visual Studio, in the Qt menu select "Open Qt project file (*.pro)",
   and open sample-framework.pro.


1. You should be able to build the project, but it won't start. Solution:
   copy QGLViewer2.dll (found in QGLViewer\QGLViewer\release\)
   into c:\Windows\System32 or into the project's directory.

1. If you also want to build a debug version of sample-framework,
   you are still not ready! You have to build debug versions of libQGLViewer
   and OpenMesh first, then change the library names in the project properties
   dialog window (and don't forget to copy QGLViewerd2.dll to a location
   in the path).
