#!/bin/bash
set -ex

FAILURE=0

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
qmake && make distclean && qmake engrid.pro && make $MAKEOPTIONS release || FAILURE=1
if [ $FAILURE -eq 1 ]
then
  echo "$MSG failed."
  exit 1
fi

MSG="Building engrid.pro.cgns release version"
echo $MSG
qmake && make distclean && qmake engrid.pro.cgns && make $MAKEOPTIONS release || FAILURE=1
if [ $FAILURE -eq 1 ]
then
  echo "$MSG failed."
  exit 1
fi

MSG="Building engrid.pro debug version"
echo $MSG
qmake && make distclean && qmake engrid.pro && make $MAKEOPTIONS debug || FAILURE=1
if [ $FAILURE -eq 1 ]
then
  echo "$MSG failed."
  exit 1
fi

MSG="Building engrid.pro.cgns debug version"
echo $MSG
qmake && make distclean && qmake engrid.pro.cgns && make $MAKEOPTIONS debug || FAILURE=1
if [ $FAILURE -eq 1 ]
then
  echo "$MSG failed."
  exit 1
fi

echo "SUCCESS: Everything compiles."
