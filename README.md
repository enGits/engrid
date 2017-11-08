# enGrid
*enGrid* is a mesh generation software with CFD applications in mind. It supports automatic prismatic boundary layer grids for Navier-Stokes simulations and has a Qt based GUI.

A [Doxygen](http://www.stack.nl/~dimitri/doxygen/index.html) created source code documentation can be found here:

http://todtnau.engits.de/engrid-doc/master/html

The documentation is updated automatically and should contain the correct documentation latest one day after a commit to the **master branch**.

## Building enGrid

Due significant changes to the *enGrid* codebase, only the master branch version of *enGrid* is actively supported.  The code now makes use of the CMake build system, which should simplify compilation.

The main dependencies for *enGrid* are:
* Qt 4
* VTK 6.\*
* CMake
* CGAL 

VTK needs to be compiled with Qt support, as *enGrid* depends on QVtkWidget.  The plan is to move to qt5 in the near future, which would also allow to upgrade to the latest version of VTK.

*enGrid* was successfully compiled on **Ubuntu 16.04 (Xenial Xerus)** with the following dependency versions:
* Qt 4.8.7
* CMake 3.5.1
* VTK 6.2
* CGAL 4.7-4

As **Ubuntu 16.04** only has Qt 5 support for VTK, VTK had to be compiled locally. The VTK build was configured using the following command:

`cmake -DCMAKE_BUILD_TYPE=Release -DVTK_Group_Qt=ON -DCMAKE_INSTALL_PREFIX=/the/path/of/your/choice`

*enGrid* can then be configured and compiled in a separate build directory using:

`ccmake ..\src`

pressing `[c]` to configure, pressing `[c]` a second time to accept the changes, and pressing `[g]` to generate the Makefiles and exit. The code can then be compiled and installed using make:

`make -j8 install`
