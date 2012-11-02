#!/usr/bin/env bash
# 
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +                                                                      +
# + This file is part of enGrid.                                         +
# +                                                                      +
# + Copyright 2008-2012 enGits GmbH                                      +
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
  echo "The script is unable to detect the location of your VTK installation."
  echo "Please set the following variables manually"
  echo " - VTKINCDIR"
  echo " - VTKLIBDIR"
}

if [ -d /usr/include/vtk-5.10 ]
then
  export VTKINCDIR=/usr/include/vtk-5.10
elif [ -d /usr/include/vtk-5.8 ]
then
  export VTKINCDIR=/usr/include/vtk-5.8
elif [ -d /usr/include/vtk-5.6 ]
then
  export VTKINCDIR=/usr/include/vtk-5.6
elif [ -d /usr/include/vtk-5.4 ]
then
  export VTKINCDIR=/usr/include/vtk-5.4
elif [ -d /usr/include/vtk ]
then
  export VTKINCDIR=/usr/include/vtk
else
  help  
fi

if [ -f /usr/lib/libvtkCommon.so ]
then
  export VTKLIBDIR=/usr/lib
elif [ -f /usr/lib64/libvtkCommon.so ]
then
  export VTKLIBDIR=/usr/lib64
else
  help
fi


export LD_LIBRARY_PATH=$QTDIR/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$VTKLIBDIR:$LD_LIBRARY_PATH
