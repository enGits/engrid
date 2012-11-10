#!/usr/bin/env bash
# 
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +                                                                      +
# + This file is part of enGrid.                                         +
# +                                                                      +
# + Copyright 2008-2012 enGits GmbH                                     +
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
  echo "CONFIGURATION = fedora"
  echo "                ubuntu"
  echo "                opensuse-11.2"
  echo "                opensuse-11.3"
  echo "                opensuse-11.4"
  echo "                opensuse-12.1"
  echo "                opensuse-12.2"
}

# Check if all parameters are present
# If no, exit
if [ $# -ne 1 ]
then
  help
else
  echo ""
  echo "This script makes use of the command 'sudo' to execute"
  echo "the system's package manager in order to install all"
  echo "required dependencies"
  echo ""
  whoami=`sudo whoami`
  if [ "$whoami" != 'root' ]
  then
    echo "You seem to not be able to execute commands as root (via sudo)."
    echo "Please make sure you have sufficient permissions; alternatively"
    echo "you can directly execute this script as root."
    echo ""
  else
    config_name=$1
    if [ $1 = 'ubuntu' ]
    then
      sudo apt-get install git-core subversion g++ libvtk5-qt4-dev qt4-dev-tools patch
    elif [ $1 = 'opensuse-11.2' ]
    then
      sudo zypper addrepo http://download.opensuse.org/repositories/science/openSUSE_11.2/ science
      sudo zypper install git-core subversion libqt4-devel make vtk-qt vtk-devel patch
    elif [ $1 = 'opensuse-11.3' ]
    then
      sudo zypper addrepo http://download.opensuse.org/repositories/science/openSUSE_11.3/ science
      sudo zypper install git-core subversion libqt4-devel make vtk-qt vtk-devel patch
    elif [ $1 = 'opensuse-11.4' ]
    then
      sudo zypper addrepo http://download.opensuse.org/repositories/science/openSUSE_11.4/ science
      sudo zypper install git-core subversion libqt4-devel make vtk-qt vtk-devel patch
    elif [ $1 = 'opensuse-12.1' ]
    then
      sudo zypper addrepo http://download.opensuse.org/repositories/science/openSUSE_11.4/ science
      sudo zypper install git-core subversion libqt4-devel make vtk-qt vtk-devel patch
    elif [ $1 = 'opensuse-12.2' ]
    then
      sudo zypper addrepo http://download.opensuse.org/repositories/science/openSUSE_12.2/ science
      sudo zypper install git-core subversion libqt4-devel make vtk-qt vtk-devel patch
    elif [ $1 = 'fedora' ]
    then
      sudo yum -y install git
      sudo yum -y install subversion
      sudo yum -y install wget
      sudo yum -y install gcc-c++
      sudo yum -y install vtk-qt
      sudo yum -y install vtk-devel
      sudo yum -y install patch
    else
      help
    fi

    for url_address in git://github.com/enGits/engrid.git \
        https://github.com/enGits/engrid.git \
        git://repo.or.cz/engrid-github.git \
        http://repo.or.cz/r/engrid-github.git; do

      if git clone $url_address engrid ; then
        break;
      else
        echo "Repository $url_address failed. Trying the next one..."
      fi
    done

    cd engrid
    git checkout -b release-1.4 remotes/origin/release-1.4
    cd src
    source scripts/setup_pathes.bash
    source scripts/build-nglib.sh
    cd libengrid
    qmake
    make -j4
    cd ..
    qmake
    make
    cd ../..
    echo ""
    echo "You can start enGrid by typing: `pwd`/engrid/run.bash (as non-root user)"
    echo "If you want to install a link in '/usr/bin' as well as an entry in the"
    echo "desktop menus, please execute 'engrid/setup_generic.bash'"
    echo ""
  fi
fi
