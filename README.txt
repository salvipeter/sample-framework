=========
= USAGE =
=========

Pressing 'h' displays a help window that contains all the necessary information.
If you want to use this as a base for your own project, examine the source code
and read the documentation of the various libraries (see "Dependencies" below).

================
= INSTALLATION =
================

Dependencies:
- Qt4 (http://qt.nokia.com/)
- libQGLViewer (http://www.libqglviewer.com/)
- OpenMesh (http://www.openmesh.org/)

Linux:
======
qmake && make

Assumes that:
- qmake (Qt 4.x) is in the path
- QGLViewer include files are in <standard include directory>/QGLViewer/
- OpenMesh  include files are in <standard include directory>/OpenMesh/
- QGLViewer library file  is  in <standard library directory>/
- OpenMesh  library files are in /usr/lib/OpenMesh/

If any of the above is not satisfied, edit sample-framework.pro accordingly.

Windows / Visual Studio:
========================
1. Install the Qt SDK, which should integrate itself into Visual Studio.

2. Download the source package for libQGLViewer, and put it somewhere,
   e.g. at c:\Program Files\libQGLViewer-2.3.16. Open Visual Studio,
   in the Qt menu select "Open Qt project file (*.pro)",
   and open the project file in the QGLViewer subdirectory
   (in my case c:\Program Files\libQGLViewer-2.3.16\libQGLViewer.pro).
   Compile a release version of the library.

3. Now open sample-framework.pro, and insert the following line:

INCLUDEPATH += "c:\Program Files\libQGLViewer-2.3.16"

  (using the correct path on your system, of course).

3. Download the source package for OpenMesh, and put it somewhere,
   e.g. at c:\Program Files\OpenMesh-2.1.1. Open the solution file
   in Visual Studio and build a release version of the core library.
   Then open sample-framework.pro, and insert the following line:

INCLUDEPATH += "c:\Program Files\OpenMesh-2.1.1\src"

  (using the correct path on your system, of course).

4. The LIBS variable of sample-framework.pro should be updated as well,
   replace the LIBS line with the following two:

LIBS += -L"c:\Program Files\libQGLViewer-2.3.16\QGLViewer\release" -lQGLViewer2
LIBS += -L"c:\Program Files\OpenMesh-2.1.1\lib" -llibOpenMeshCore

5. Open Visual Studio, in the Qt menu select "Open Qt project file (*.pro)",
   and open sample-framework.pro.

6. In the project's properties, under C/C++ / Preprocessor, add the following
   preprocessor definitions: _USE_MATH_DEFINES, NOMINMAX

7. You should be able to build the project, but it won't start. Solution:
   copy QGLViewer2.dll (found in QGLViewer-2.3.16\QGLViewer\release\)
   into c:\Windows\System32 or into the project's directory.

8. If you also want to build a debug version of sample-framework,
   you are still not ready! You have to build debug versions of libQGLViewer
   and OpenMesh first, then change the library names in the project properties
   dialog window (and don't forget to copy QGLViewerd2.dll to a location
   in the path).
