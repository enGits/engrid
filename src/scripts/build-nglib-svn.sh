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

# DESCRIPTION:
# This script checks out or updates the netgen source code and creates the static netgen library.

# cf parameter expansion in bash manual
SCRIPTDIR=${0%$(basename $0)}
cd $SCRIPTDIR/..

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
    echo "starting qmake for $package"
    echo

    qmake

    echo
    echo "making $package"
    echo

    make
)

# ----------------------------------------------------------------- end-of-file
