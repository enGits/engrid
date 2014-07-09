#!/usr/bin/env bash
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +                                                                      +
# + This file is part of enGrid.                                         +
# +                                                                      +
# + Copyright 2008-2014 enGits GmbH                                      +
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

cd ${0%/*} || exit 1    # run from this directory

usage() {
    while [ "$#" -ge 1 ]; do echo "$1"; shift; done
    cat<<USAGE

usage: ${0##*/} [OPTION]
options:
  -help

* compile engrid with the OpenFOAM paraview libraries

USAGE
    exit 1
}

#------------------------------------------------------------------------------


# set the major version "<digits>.<digits>"
ParaView_MAJOR_VERSION=$(echo $ParaView_VERSION | \
    sed -e 's/^\([0-9][0-9]*\.[0-9][0-9]*\).*$/\1/')

libdir="$ParaView_DIR/lib/paraview-$ParaView_MAJOR_VERSION"

[ -d "$ParaView_INST_DIR" ] || usage "ParaView_INST_DIR not found ($ParaView_INST_DIR)"
[ -d "$ParaView_DIR" ] || usage "ParaView_DIR not found ($ParaView_DIR)"

[ -d $libdir ] || usage "paraview libraries not found"

export LD_LIBRARY_PATH=$libdir:$LD_LIBRARY_PATH

qmake engrid.pro.OpenFOAM-paraview

make
make install

# ----------------------------------------------------------------- end-of-file
