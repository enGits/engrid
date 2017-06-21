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

# DESCRIPTION:

cd netgen_svn
wget https://sourceforge.net/projects/netgen-mesher/files/netgen-mesher/4.9.13/netgen-4.9.13.zip
unzip netgen-4.9.13.zip
rm -f netgen-4.9.13.zip
[ -e netgen-mesher ] && rm -rf netgen-mesher
mkdir netgen-mesher
mv netgen-4.9.13 netgen-mesher/netgen
cd netgen-mesher
patch -p0 < ../nglib_engrid_mods.diff
cd ..
qmake
make clean
make -j4
cd ..

# ----------------------------------------------------------------- end-of-file
