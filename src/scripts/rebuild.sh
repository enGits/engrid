#!/usr/bin/env bash
# 
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +                                                                      +
# + This file is part of enGrid.                                         +
# +                                                                      +
# + Copyright 2008-2011 enGits GmbH                                     +
# +                                                                      +
# + enGrid is free software: you can redistribute it and/or modify       +
# + it under the terms of the GNU General Public License as published by +
# + the Free Software Foundation, either version 3 of the License, or    +
# + (at your option) any later version.                                  +
# +                                                                      +
# + enGrid is distributed in the hope that it will be useful,            +
# + but WITHOUT ANY WARRANTY; without even the implied warranty of       +
# + MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        +
# + GNU General Public License for more details.                         +
# +                                                                      +
# + You should have received a copy of the GNU General Public License    +
# + along with enGrid. If not, see <http://www.gnu.org/licenses/>.       +
# +                                                                      +
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 

set -ex

FAILURE=0

#set up environment
# source ./scripts/setup_paths.sh engits yes
source /opt/engits/bin/engrid_environment.sh

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
./scripts/build-nglib.sh

MAKEOPTIONS=""

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

echo "SUCCESS: Everything compiles."
