#!/bin/bash

# Check if all parameters are present
# If no, exit
if [ $# -ne 2 ]
then
        echo "usage :"
        echo "source `basename $0` CONFIGURATION ENABLE_CGNS"
        echo "CONFIGURATION = engits / debian / ubuntu / opensuse / debian_etch"
        echo "ENABLE_CGNS = yes/no"
        exit 0
fi

#Debian
if [ $1 = 'debian' ]
then
  export VTKLIBDIR=/usr/lib/
  export VTKINCDIR=/usr/include/vtk-5.0/
fi

#Debian
if [ $1 = 'debian_etch' ]
then
  export PATH=/usr/local/Trolltech/Qt-4.5.1/bin:$PATH
  export QTDIR=/usr/local/Trolltech/Qt-4.5.1/
  export LD_LIBRARY_PATH=/usr/local/Trolltech/Qt-4.5.1/lib:$LD_LIBRARY_PATH
  export VTKLIBDIR=/usr/lib/
  export VTKINCDIR=/usr/include/vtk-5.0/
fi

#Ubuntu
if [ $1 = 'ubuntu' ]
then
  export VTKLIBDIR=/usr/lib/
  export VTKINCDIR=/usr/include/vtk-5.0/
fi

#OpenSUSE
if [ $1 = 'opensuse' ]
then
  export VTKLIBDIR=/usr/lib64
  export VTKINCDIR=/usr/include/vtk
fi

#enGits local system
if [ $1 = 'engits' ]
then
  export QTDIR=/opt/shared/Qt/4.5.1/debug
  export VTKLIBDIR=/opt/shared/VTK/lib/vtk-5.4/
  export VTKINCDIR=/opt/shared/VTK/include/vtk-5.4/
  if [ $2 = 'yes' ]
  then
    export CGNSINCDIR=/opt/shared/cgns/include/
    export CGNSLIBDIR=/opt/shared/cgns/lib/
  fi
fi

export LD_LIBRARY_PATH=$QTDIR/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$VTKLIBDIR:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$CGNSLIBDIR:$LD_LIBRARY_PATH
