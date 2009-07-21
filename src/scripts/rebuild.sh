#!/bin/bash
set -ex

#set up environment
source ./scripts/setup_paths.sh engits yes

echo "BUILDING TOOLS":
echo "gcc = $(which gcc)"
echo "g++ = $(which g++)"
echo "qmake = $(which qmake)"
echo "make = $(which make)"
gcc -v
g++ -v
qmake -v
make -v

echo "PATHS:"
echo QTDIR = $QTDIR
echo VTKLIBDIR = $VTKLIBDIR
echo VTKINCDIR = $VTKINCDIR
echo CGNSINCDIR = $CGNSINCDIR
echo CGNSLIBDIR = $CGNSLIBDIR
echo PATH = $PATH
echo LD_LIBRARY_PATH = $LD_LIBRARY_PATH

echo "Building netgen"
./build-nglib.sh

MAKEOPTIONS=""

MSG="Building engrid.pro release version"
echo $MSG
qmake && make distclean && qmake engrid.pro && make $MAKEOPTIONS release || ( echo "$MSG failed." && exit 1 )

MSG="Building engrid.pro.cgns release version"
echo $MSG
qmake && make distclean && qmake engrid.pro.cgns && make $MAKEOPTIONS release || ( echo "$MSG failed." && exit 1 )

MSG="Building engrid.pro debug version"
echo $MSG
qmake && make distclean && qmake engrid.pro && make $MAKEOPTIONS debug || ( echo "$MSG failed." && exit 1 )

MSG="Building engrid.pro.cgns debug version"
echo $MSG
qmake && make distclean && qmake engrid.pro.cgns && make $MAKEOPTIONS debug || ( echo "$MSG failed." && exit 1 )

echo "SUCCESS: Everything compiles."
