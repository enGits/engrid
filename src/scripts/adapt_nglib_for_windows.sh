#!/usr/bin/env bash
# 
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +                                                                      +
# + This file is part of enGrid.                                         +
# +                                                                      +
# + Copyright 2008-2010 enGits GmbH                                     +
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


echo ${0%/*} #directory of script
echo ${0%/*/*} # parent directory of script
# exit 0
cd ${0%/*/*} || exit 1    # run from parent directory of script

package=netgen-mesher
(
    echo "Working directory = $(pwd)"
    cd netgen_svn || exit 1

    if [ -d netgen-mesher/.svn ]
    then
        echo "updating NETGEN from SVN repository (sourceforge.net) -- please wait"
        svn up $package
    else
        echo "downloading NETGEN from SVN repository (sourceforge.net) -- please wait"
        svn co https://netgen-mesher.svn.sourceforge.net/svnroot/$package $package
    fi

    echo
    echo "adapting nglib.h for static library on windows"
    echo

    sed -i 's/__declspec(dllexport)//' ./netgen-mesher/netgen/nglib/nglib.h
    sed -i 's/__declspec(dllimport)//' ./netgen-mesher/netgen/nglib/nglib.h
)

# ----------------------------------------------------------------- end-of-file
