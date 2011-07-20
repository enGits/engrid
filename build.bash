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

help ()
{
  echo "usage :"
  echo "`basename $0` CONFIGURATION"
  echo "CONFIGURATION = ubuntu-10.04"
  echo "                ubuntu-11.04"
  echo "                opensuse-11.4"
  exit 0
}

# Check if all parameters are present
# If no, exit
if [ $# -ne 1 ]
then
  help
fi

# Ubuntu
if [ $1 = 'ubuntu-10.04' ]
then
  export VTKLIBDIR=/usr/lib/
  export VTKINCDIR=/usr/include/vtk-5.4/
  sudo apt-get install git-core subversion libvtk5-qt4-dev qt4-dev-tools
elif [ $1 = 'ubuntu-11.04' ]
then
  export VTKLIBDIR=/usr/lib/
  export VTKINCDIR=/usr/include/vtk-5.4/
  sudo apt-get install git subversion libvtk5-qt4-dev qt4-dev-tools
elif [ $1 = 'opensuse-11.4' ]
then
#  export VTKLIBDIR=/usr/lib64
  export VTKINCDIR=/usr/include/vtk-5.6
  sudo zypper addrepo http://download.opensuse.org/repositories/science/openSUSE_11.4/ science
  sudo zypper install git-core subversion libqt4-devel make vtk-qt vtk-devel
else
  help
fi

export LD_LIBRARY_PATH=$QTDIR/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$VTKLIBDIR:$LD_LIBRARY_PATH

git clone git://engrid.git.sourceforge.net/gitroot/engrid/engrid
cd engrid
git checkout -b release-1.3 remotes/origin/release-1.3
cd src
source scripts/setup_pathes.bash $1
scripts/build-nglib.sh
cd libengrid
qmake
make -j4
cd ..
qmake
make
cd ../..
echo "You can start enGrid by typing: `pwd`/engrid/run.bash (as non-root user)"

