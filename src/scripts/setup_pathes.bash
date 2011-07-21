#!/usr/bin/env bash
# 
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +                                                                      +
# + This file is part of enGrid.                                         +
# +                                                                      +
# + Copyright 2008-2011 enGits GmbH                                     +
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
  echo "source `basename $0` CONFIGURATION"
  echo "CONFIGURATION = ubuntu"
  echo "                opensuse32"
  echo "                opensuse64"
  exit 0
}

# Check if all parameters are present
# If no, exit
if [ $# -ne 1 ]
then
  help
fi

if [ $1 = 'ubuntu' ]
then
  export VTKINCDIR=/usr/include/vtk-5.4/
  export VTKLIBDIR=/usr/lib
elif [ $1 = 'opensuse32' ]
then
  export VTKINCDIR=/usr/include/vtk-5.6
  export VTKLIBDIR=/usr/lib
elif [ $1 = 'opensuse64' ]
then
  export VTKINCDIR=/usr/include/vtk-5.6
  export VTKLIBDIR=/usr/lib64
else
  help
fi

export LD_LIBRARY_PATH=$QTDIR/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$VTKLIBDIR:$LD_LIBRARY_PATH
