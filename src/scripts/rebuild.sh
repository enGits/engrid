#!/bin/bash
set -eux

#set up environment
PATH=/usr/bin:/opt/shared/Qt/4.5.0/debug/bin:$PATH
source ./scripts/setup_paths.cgns.sh engits

echo "PATHS:"
echo $QTDIR
echo $VTKLIBDIR
echo $VTKINCDIR
echo $CGNSINCDIR
echo $CGNSLIBDIR
echo $LD_LIBRARY_PATH

echo "Building engrid.pro release version"
qmake && make distclean && qmake engrid.pro && make -j2 || exit 1
echo "Building engrid.pro.cgns release version"
qmake && make distclean && qmake engrid.pro.cgns && make -j2 || exit 1

echo "Building engrid.pro debug version"
qmake && make distclean && qmake engrid.pro && make -j2 debug || exit 1
echo "Building engrid.pro.cgns debug version"
qmake && make distclean && qmake engrid.pro.cgns && make -j2 debug || exit 1

echo "SUCCESS: Everything compiles."
