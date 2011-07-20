#!/usr/bin/env bash
olddir=`pwd`
cd `dirname $0`
sudo zypper addrepo http://download.opensuse.org/repositories/science/openSUSE_11.4/ science
sudo zypper install vtk-qt
echo "You can start enGrid by typing: `pwd`/run.bash (as non-root user)"
cd $olddir

