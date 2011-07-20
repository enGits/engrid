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
  echo "CONFIGURATION = ubuntu"
  echo "                opensuse-11.2"
  echo "                opensuse-11.3"
  echo "                opensuse-11.4"
  exit 0
}

# Check if all parameters are present
# If no, exit
if [ $# -ne 1 ]
then
  help
fi

config_name = $1

# Ubuntu
if [ $1 = 'ubuntu' ]
then
  sudo apt-get install git-core subversion libvtk5-qt4-dev qt4-dev-tools
elif [ $1 = 'opensuse-11.2' ]
then
  sudo zypper addrepo http://download.opensuse.org/repositories/science/openSUSE_11.2/ science
  sudo zypper install git-core subversion libqt4-devel make vtk-qt vtk-devel
  config_name="opensuse"
elif [ $1 = 'opensuse-11.3' ]
then
  sudo zypper addrepo http://download.opensuse.org/repositories/science/openSUSE_11.3/ science
  sudo zypper install git-core subversion libqt4-devel make vtk-qt vtk-devel
  config_name="opensuse"
elif [ $1 = 'opensuse-11.4' ]
then
  sudo zypper addrepo http://download.opensuse.org/repositories/science/openSUSE_11.4/ science
  sudo zypper install git-core subversion libqt4-devel make vtk-qt vtk-devel
  config_name="opensuse"
else
  help
fi

git clone git://engrid.git.sourceforge.net/gitroot/engrid/engrid
cd engrid
git checkout -b release-1.3 remotes/origin/release-1.3
cd src
source scripts/setup_pathes.bash $config_name
scripts/build-nglib.sh
cd libengrid
qmake
make -j4
cd ..
qmake
make
cd ../..
echo "You can start enGrid by typing: `pwd`/engrid/run.bash (as non-root user)"

