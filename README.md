# enGrid
*enGrid* is a mesh generation software with CFD applications in mind. It supports automatic prismatic boundary layer grids for Navier-Stokes simulations and has a Qt based GUI.

A [Doxygen](http://www.stack.nl/~dimitri/doxygen/index.html) created source code documentation can be found here:

http://todtnau.engits.de/engrid-doc/master/html

The documentation is updated automatically and should contain the correct documentation latest one day after a commit to the **master branch**.

## Building enGrid

Due significant changes to the *enGrid* codebase, only the master branch version of *enGrid* is actively supported.  The code now makes use of the CMake build system, which should simplify compilation.

The main dependencies for *enGrid* are:

+ Qt 4
+ VTK 6.\*
+ CMake
+ CGAL 

VTK needs to be compiled with Qt support, as *enGrid* depends on QVtkWidget.  The plan is to move to qt5 in the near future, which would also allow to upgrade to the latest version of VTK.

*enGrid* was successfully compiled on **Ubuntu 16.04 (Xenial Xerus)** with the following dependency versions:

+ Qt 4.8.7
+ CMake 3.5.1
+ VTK 6.2
+ CGAL 4.7-4

As **Ubuntu 16.04** only has Qt 5 support for VTK, VTK had to be compiled locally. The VTK build was configured using the following command:

`cmake -DCMAKE_BUILD_TYPE=Release -DVTK_Group_Qt=ON -DCMAKE_INSTALL_PREFIX=/the/path/of/your/choice`

*enGrid* can then be configured and compiled in a separate build directory using:

`ccmake ..\src`

pressing `[c]` to configure, pressing `[c]` a second time to accept the changes, and pressing `[g]` to generate the Makefiles and exit. The code can then be compiled and installed using make:

`make -j8 install`


### Detailed building procedure on Ubuntu 18.04.3 LTS (Bionic Beaver)

The steps to complete the compilation are:

1. Installation of previous compiler version and cmake with interface
2. Setup to work with previous compiler version
3. Compile VTK 6.3 with flags and dependencies
4. Compile libboost 1.58.0
5. Compile CGAL 4.7 with flags and dependencies
6. Compile *enGrid* with the necessary flags


#### 1. Installation of previous compiler version and cmake with interface

`sudo apt-get install cmake-curses-gui gcc-5 g++-5`


#### 2. Setup to work with previous compiler version
In Ubuntu the command 'update-alternatives' allows to setup a different version of the gcc compiler:

    sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-7 10
    sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-5 20
    sudo update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-5 20

    sudo update-alternatives --install /usr/bin/cc cc /usr/bin/gcc 30
    sudo update-alternatives --set cc /usr/bin/gcc

    sudo update-alternatives --install /usr/bin/c++ c++ /usr/bin/g++ 30
    sudo update-alternatives --set c++ /usr/bin/g++


#### 3. Compile VTK 6.3 with flags and dependencies
Install dependencies of VTK 6.3 :

`sudo apt-get install libxt-dev qt4-default libqtwebkit-dev`

Download and extract the [source code of VTK 6.3](https://www.vtk.org/files/release/6.3/VTK-6.3.0.tar.gz) on a folder to store *enGrid* dependencies. Then proceed with the build:

    cd your/path/to/VTK-6.3.0
    mkdir build
    cd build
    ccmake ../

On the cmake interface press `[c]` to configure and `[e]` to exit if any message is shown. Then look for the tag `VTK_Group_Qt` and press enter to change it from `OFF` to `ON`, then press `[c]` and `[c]` again to complete the setup and `[g]` to generate the files for compiling. Finally build VTK with the command `make` with the `-j` option followed by the number of cores of the computer:

`make -j4`


#### 4. Compile libboost 1.58.0

Download and extract the [source code of Boost C++ Libraries 1.58.0](https://sourceforge.net/projects/boost/files/boost/1.58.0/boost_1_58_0.tar.gz/download). Then proceed with the build:

    cd your/path/to/boost_1_58_0
    ./bootstrap.sh --with-python=python2
    ./b2

#### 5. Compile CGAL 4.7 with flags and dependencies


Install dependencies of CGAL 4.7 :

`sudo apt-get install libgmp-dev libmpfr-dev`

Download and extract the [source code of CGAL 4.7](https://github.com/CGAL/cgal/releases/download/releases%2FCGAL-4.7/CGAL-4.7.tar.gz). Then proceed with the build:

    cd your/path/to/CGAL-4.7
    mkdir build
    cd build
    ccmake ../


On the cmake interface press `[c]` to configure and `[e]` to exit if any message is shown. Then press `[t]` to toggle advanced mode and look and set the following tags:

+ `Boost_INCLUDE_DIR` set to `/absolute/path/to/boost_1_58_0/`
+ `Boost_LIBRARY_DIR_RELEASE` set to `/absolute/path/to/boost_1_58_0/stage/lib`
+ `WITH_CGAL_Qt5` set to `OFF`

Then press `[c]` and `[e]` to complete the setup and `[g]` to generate the files for compiling. Finally build CGAL:

`make -j4`


#### 6. Compile *enGrid* with the necessary flags

Download or clone and extract the [*enGrid* source code repository](https://github.com/dinlink/engrid/archive/master.zip). Then proceed with the build:

    cd your/path/to/engrid-master/
    mkdir build
    cd build
    ccmake ../src/

On the cmake interface:

1. press `[c]`, `[e]`
2. set `VTK_DIR` to `/absolute/path/to/VTK-6.3/build/`
3. press `[c]`, `[e]`
4. set `CGAL_DIR` to `/absolute/path/to/CGAL-4.7/build/`
5. press `[c]`
6. set `CGAL_INCLUDE_PATH` to `/absolute/path/to/CGAL-4.7/build/include/`
7. press `[c]`
6. set `BOOST_INCLUDE_PATH` to `/absolute/path/to/boost_1_58_0/`
7. press `[c]`, `[g]`

Build *enGrid* 

`make -j4`
