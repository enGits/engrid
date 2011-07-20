#!/usr/bin/env bash
olddir=`pwd`
cd `dirname $0`
sudo zypper addrepo http://download.opensuse.org/repositories/science/openSUSE_11.4/ science
sudo zypper install vtk-qt
cp engrid iengrid
echo "`pwd`/run.bash" >> iengrid
sudo cp iengrid /usr/bin/engrid
cp engrid.desktop iengrid.desktop
echo "Icon=`pwd`/src/libengrid/resources/icons/G.png" >> iengrid.desktop
sudo cp iengrid.desktop /usr/share/applications/engrid.desktop
cd $olddir
echo "You can start enGrid by typing 'engrid' from the command line (as non-root user) or with the help of the desktop menus."
