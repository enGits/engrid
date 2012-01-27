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
  echo "CONFIGURATION = fedora-15-32"
  echo "                fedora-15-64"
  echo "                fedora-14-32"
  echo "                fedora-14-64"
  echo "                ubuntu-10.10"
  echo "                ubuntu-11.04"
  echo "                ubuntu-11.10"
  echo "                opensuse-11.2-32"
  echo "                opensuse-11.2-64"
  echo "                opensuse-11.3-32"
  echo "                opensuse-11.3-64"
  echo "                opensuse-11.4-32"
  echo "                opensuse-11.4-64"
  echo "                opensuse-12.1-32"
  echo "                opensuse-12.1-64"
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
    if [ $1 = 'ubuntu-10.10' ]
    then
      sudo apt-get install git-core subversion libvtk5-qt4-dev qt4-dev-tools
    elif [ $1 = 'ubuntu-11.04' ]
    then
      sudo apt-get install git-core subversion libvtk5-qt4-dev qt4-dev-tools
    elif [ $1 = 'ubuntu-11.10' ]
    then
      sudo apt-get install git-core subversion g++ libvtk5-qt4-dev qt4-dev-tools
    elif [ $1 = 'opensuse-11.2-32' ]
    then
      sudo zypper addrepo http://download.opensuse.org/repositories/science/openSUSE_11.2/ science
      sudo zypper install git-core subversion libqt4-devel make vtk-qt vtk-devel
      config_name="opensuse32"
    elif [ $1 = 'opensuse-11.3-32' ]
    then
      sudo zypper addrepo http://download.opensuse.org/repositories/science/openSUSE_11.3/ science
      sudo zypper install git-core subversion libqt4-devel make vtk-qt vtk-devel
      config_name="opensuse32"
    elif [ $1 = 'opensuse-11.4-32' ]
    then
      sudo zypper addrepo http://download.opensuse.org/repositories/science/openSUSE_11.4/ science
      sudo zypper install git-core subversion libqt4-devel make vtk-qt vtk-devel
      config_name="opensuse32"
    elif [ $1 = 'opensuse-12.1-32' ]
    then
      sudo zypper addrepo http://download.opensuse.org/repositories/science/openSUSE_11.4/ science
      sudo zypper install git-core subversion libqt4-devel make vtk-qt vtk-devel
      config_name="opensuse32-12"
    elif [ $1 = 'opensuse-11.2-64' ]
    then
      sudo zypper addrepo http://download.opensuse.org/repositories/science/openSUSE_11.2/ science
      sudo zypper install git-core subversion libqt4-devel make vtk-qt vtk-devel
      config_name="opensuse64"
    elif [ $1 = 'opensuse-11.3-64' ]
    then
      sudo zypper addrepo http://download.opensuse.org/repositories/science/openSUSE_11.3/ science
      sudo zypper install git-core subversion libqt4-devel make vtk-qt vtk-devel
      config_name="opensuse64"
    elif [ $1 = 'opensuse-11.4-64' ]
    then
      sudo zypper addrepo http://download.opensuse.org/repositories/science/openSUSE_11.4/ science
      sudo zypper install git-core subversion libqt4-devel make vtk-qt vtk-devel
      config_name="opensuse64"
    elif [ $1 = 'opensuse-12.1-64' ]
    then
      sudo zypper addrepo http://download.opensuse.org/repositories/science/openSUSE_11.4/ science
      sudo zypper install git-core subversion libqt4-devel make vtk-qt vtk-devel
      config_name="opensuse64-12"
    elif [ $1 = 'fedora-15-32' ]
    then
      sudo yum -y install git
      sudo yum -y install subversion
      sudo yum -y install wget
      sudo yum -y install gcc-c++
      sudo yum -y install vtk-qt
      config_name="fedora32"
    elif [ $1 = 'fedora-15-64' ]
    then
      sudo yum -y install git
      sudo yum -y install subversion
      sudo yum -y install wget
      sudo yum -y install gcc-c++
      sudo yum -y install vtk-qt
      config_name="fedora64"
    elif [ $1 = 'fedora-14-64' ]
    then
      sudo yum -y install git
      sudo yum -y install subversion
      sudo yum -y install wget
      sudo yum -y install gcc-c++
      sudo yum -y install vtk-qt
      config_name="fedora64"
    else
      help
    fi
    #git clone git://engrid.git.sourceforge.net/gitroot/engrid/engrid
    git clone http://repo.or.cz/r/engrid.git
    echo $config_name > engrid/config.txt
    cd engrid
    git checkout -b release-1.3 remotes/origin/release-1.3
    cd src
    source scripts/setup_pathes.bash $config_name
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
