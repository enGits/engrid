#!/bin/bash
set -ex

#set up environment
PATH=/usr/bin:/opt/shared/Qt/4.5.1/debug/bin:$PATH
source ./scripts/setup_paths.sh engits yes

echo "BUILDING TOOLS":
echo "gcc= $(which gcc)"
echo "g++= $(which g++)"
echo "qmake= $(which qmake)"
echo "make= $(which make)"
gcc -v
g++ -v
qmake -v
make -v

echo "PATHS:"
echo PATH = $PATH
echo QTDIR = $QTDIR
echo VTKLIBDIR = $VTKLIBDIR
echo VTKINCDIR = $VTKINCDIR
echo CGNSINCDIR = $CGNSINCDIR
echo CGNSLIBDIR = $CGNSLIBDIR
echo LD_LIBRARY_PATH = $LD_LIBRARY_PATH

echo "Building netgen"
./build-nglib.sh

echo "Building engrid.pro release version"
qmake && make distclean && qmake engrid.pro && make -j2 || exit 1

echo "Building engrid.pro.cgns release version"
qmake && make distclean && qmake engrid.pro.cgns && make -j2 || exit 1

echo "Building engrid.pro debug version"
qmake && make distclean && qmake engrid.pro && make -j2 debug || exit 1

echo "Building engrid.pro.cgns debug version"
qmake && make distclean && qmake engrid.pro.cgns && make -j2 debug || exit 1

echo "SUCCESS: Everything compiles."
