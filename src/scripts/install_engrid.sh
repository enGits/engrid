#!/usr/bin/env bash

#for debugging
set -eux

source "${0%/*}/engrid_installer_updater.cfg"

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
  qmake && make release
  cd -
}

update_netgen()
{
  cd ./engrid/src
  echo "Update netgen"
  ./scripts/build-nglib.sh
  cd -
}

update_engrid()
{
  cd ./engrid/src
  echo "Update enGrid"
  git pull
  qmake && make release
  cd -
}

rebuild_engrid()
{
  cd ./engrid/src
  qmake && make distclean && qmake && make release
  cd -
}

ans=$(zenity  --list  --text "Which actions should be executed?" --checklist  --column "Run" --column "Actions" \
FALSE "create_bash_engrid" \
FALSE "install_QT" \
FALSE "install_VTK" \
FALSE "install_CGNS" \
FALSE "build_engrid" \
FALSE "update_netgen" \
FALSE "update_engrid" \
FALSE "rebuild_engrid" \
--separator=":");
echo $ans

if ( echo $ans | grep "create_bash_engrid" ) then create_bash_engrid; fi
source $ENV_SETUP
if ( echo $ans | grep "install_QT" ) then install_QT; fi;
if ( echo $ans | grep "install_VTK" ) then install_VTK; fi;
if ( echo $ans | grep "install_CGNS" ) then install_CGNS; fi;
if ( echo $ans | grep "build_engrid" ) then build_engrid; fi;
if ( echo $ans | grep "update_netgen" ) then update_netgen; fi;
if ( echo $ans | grep "update_engrid" ) then update_engrid; fi;
if ( echo $ans | grep "rebuild_engrid" ) then rebuild_engrid; fi;

echo "SUCCESS"
exit 0
