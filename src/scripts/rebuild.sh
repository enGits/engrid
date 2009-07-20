#!/bin/bash
set -ex

#set up environment
PATH=/usr/bin:/opt/shared/Qt/4.5.1/debug/bin:$PATH
source ./scripts/setup_paths.cgns.sh engits

echo "PATHS:"
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
