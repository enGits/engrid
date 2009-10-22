#!/usr/bin/env bash

#for debugging
set -eux

#variables to set
ARCHIVE_QT="qt-x11-opensource-src-4.5.2.tar.gz"
URL_QT="http://get.qt.nokia.com/qt/source/$ARCHIVE_QT"

ARCHIVE_VTK="vtk-5.4.2.tar.gz"
URL_VTK="http://www.vtk.org/files/release/5.4/$ARCHIVE_VTK"

ARCHIVE_CGNS="cgnslib_2.5-4.tar.gz"
URL_CGNS="http://prdownloads.sourceforge.net/cgns/$ARCHIVE_CGNS"

#HTTP
# URL_ENGRID="http://engits.com/git/engrid.git"
#SSH
URL_ENGRID="ssh://swordfish/srv/www/htdocs/git/engrid.git"

PREFIX="$(readlink -f .)/usr/engrid_libs"
mkdir -p $PREFIX

VTKPREFIX=$PREFIX/VTK
VTKVERSION="5.4"
QTPREFIX=$PREFIX/QT
CGNSPREFIX=$PREFIX/CGNS
# VTK_WRAP_JAVA
# VTK_WRAP_PYTHON
# VTK_WRAP_TCL

ENV_SETUP="engrid_environment.sh"

install_QT()
{
  echo "Install QT"
#   wget $URL_QT
  tar -xzvf ./$ARCHIVE_QT
  cd $(basename $ARCHIVE_QT .tar.gz)
  mkdir -p $QTPREFIX
  echo yes | ./configure --prefix=$QTPREFIX -opensource
  make && make install
  cd -
}

install_VTK()
{
  echo "Install VTK"
#   wget $URL_VTK
  tar -xzvf ./$ARCHIVE_VTK
  cd ./VTK
  mkdir -p $VTKPREFIX
  cmake -DCMAKE_INSTALL_PREFIX:PATH=$VTKPREFIX -DBUILD_SHARED_LIBS:BOOL=ON -DVTK_USE_GUISUPPORT:BOOL=ON -DVTK_USE_QVTK:BOOL=ON -DDESIRED_QT_VERSION:STRING=4  .
  chmod 644 Utilities/vtktiff/tif_fax3sm.c
  make && make install
  cd -
}

install_CGNS()
{
  echo "Install CGNS"
  wget $URL_CGNS
  tar -xzvf ./$ARCHIVE_CGNS
  cd ./cgnslib_2.5/
  mkdir -p $CGNSPREFIX
  mkdir -p $CGNSPREFIX/include
  mkdir -p $CGNSPREFIX/lib
  ./configure --prefix=$CGNSPREFIX && make && make install
  cd -
}

build_engrid()
{
  git clone $URL_ENGRID
  cd engrid/src
  echo "Build netgen"
  ./scripts/build-nglib.sh
  echo "Build enGrid"
  qmake && make distclean && qmake && make release
  cd -
}

update_netgen()
{
  echo "Update netgen"
  ./scripts/build-nglib.sh
}

update_engrid()
{
  echo "Update enGrid"
  git pull
  qmake && make distclean && qmake && make release
}

create_bash_engrid()
{
  echo "Create bash_engrid"

  echo "#!/usr/bin/env bash" > $ENV_SETUP
  echo "export VTKINCDIR=$VTKPREFIX/include/vtk-$VTKVERSION" >> $ENV_SETUP
  echo "export VTKLIBDIR=$VTKPREFIX/lib/vtk-$VTKVERSION" >> $ENV_SETUP
  echo "export LD_LIBRARY_PATH=$VTKLIBDIR:\$LD_LIBRARY_PATH" >> $ENV_SETUP
  
  echo "export CGNSINCDIR=/opt/shared/cgns/include" >> $ENV_SETUP
  echo "export CGNSLIBDIR=/opt/shared/cgns/lib" >> $ENV_SETUP
  echo "export LD_LIBRARY_PATH=$CGNSLIBDIR:\$LD_LIBRARY_PATH" >> $ENV_SETUP
  
  echo "export PATH=$QTPREFIX/bin:\$PATH" >> $ENV_SETUP
  echo "export QTDIR=$QTPREFIX" >> $ENV_SETUP
  echo "export LD_LIBRARY_PATH=$QTPREFIX/lib:\$LD_LIBRARY_PATH" >> $ENV_SETUP
}

create_bash_engrid
source $ENV_SETUP
# install_QT
# install_VTK
# install_CGNS
build_engrid
cd engrid/src
update_netgen
update_engrid
cd -

exit 0
