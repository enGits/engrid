#!/usr/bin/env bash
olddir=`pwd`
cd `dirname $0`
sudo zypper addrepo http://download.opensuse.org/repositories/science/openSUSE_11.4/ science
sudo zypper install vtk-qt
echo "`pwd`/run.bash" >> engrid
sudo cp engrid /usr/bin
cd $olddir
echo "You can start enGrid by typing 'engrid' as non-root user."
