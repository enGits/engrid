#!/usr/bin/env bash
echo "This script makes use of the command 'sudo' to install"
echo "a script in the directory '/usr/bin' as well as a desktop"
echo "file in the directory '/usr/share/applications'."
echo ""
if [ `sudo whoami` != 'root' ]
  echo "You seem to not be able to execute commands as root (via sudo)."
  echo "Please make sure you have sufficient permissions; alternatively"
  echo "you can directly execute this script as root."
else
  olddir=`pwd`
  cd `dirname $0`
  cp engrid iengrid
  echo "`pwd`/run.bash" >> iengrid
  sudo cp iengrid /usr/bin/engrid
  cp engrid.desktop iengrid.desktop
  echo "Icon=`pwd`/src/libengrid/resources/icons/G.png" >> iengrid.desktop
  sudo cp iengrid.desktop /usr/share/applications/engrid.desktop
  cd $olddir
  echo "You can start enGrid by typing 'engrid' from the command line (as non-root user)"
  echo "or with the help of the desktop menus."
fi