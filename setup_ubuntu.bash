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
olddir=`pwd`
cd `dirname $0`
sudo apt-get install libvtk5-qt4-dev qt4-dev-tools
cp engrid iengrid
echo "`pwd`/run.bash" >> iengrid
sudo cp iengrid /usr/bin/engrid
cp engrid.desktop iengrid.desktop
echo "Icon=`pwd`/src/libengrid/resources/icons/G.png" >> iengrid.desktop
sudo cp iengrid.desktop /usr/share/applications/engrid.desktop
cd $olddir
echo "You can start enGrid by typing 'engrid' from the command line (as non-root user) or with the help of the desktop menus."
