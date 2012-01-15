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

#for debugging
set -eux

source "${0%/*}/engrid_installer_updater.cfg"

VTKINCDIR=$VTKPREFIX/include/vtk-$VTKVERSION
VTKLIBDIR=$VTKPREFIX/lib/vtk-$VTKVERSION

CGNSINCDIR=$CGNSPREFIX/include
CGNSLIBDIR=$CGNSPREFIX/lib

create_bash_engrid()
{
  echo "Create bash_engrid"
  mkdir -p $BINPREFIX

  echo "#!/usr/bin/env bash" > $BINPREFIX/$ENV_SETUP
  echo "export VTKINCDIR=$VTKINCDIR" >> $BINPREFIX/$ENV_SETUP
  echo "export VTKLIBDIR=$VTKLIBDIR" >> $BINPREFIX/$ENV_SETUP
  echo "export LD_LIBRARY_PATH=$VTKLIBDIR:\$LD_LIBRARY_PATH" >> $BINPREFIX/$ENV_SETUP
  
  echo "export CGNSINCDIR=$CGNSINCDIR" >> $BINPREFIX/$ENV_SETUP
  echo "export CGNSLIBDIR=$CGNSLIBDIR" >> $BINPREFIX/$ENV_SETUP
  echo "export LD_LIBRARY_PATH=$CGNSLIBDIR:\$LD_LIBRARY_PATH" >> $BINPREFIX/$ENV_SETUP
  
  echo "export PATH=$QTPREFIX/bin:\$PATH" >> $BINPREFIX/$ENV_SETUP
  echo "export QTDIR=$QTPREFIX" >> $BINPREFIX/$ENV_SETUP
  echo "export LD_LIBRARY_PATH=$QTPREFIX/lib:\$LD_LIBRARY_PATH" >> $BINPREFIX/$ENV_SETUP

  chmod 755 $BINPREFIX/$ENV_SETUP
}

create_start_engrid()
{
  echo "Create start_engrid"
  mkdir -p $BINPREFIX

  echo "#!/usr/bin/env bash" > $BINPREFIX/$START_ENGRID
  echo "source $BINPREFIX/$ENV_SETUP" >> $BINPREFIX/$START_ENGRID
  echo "$SRCPREFIX/engrid/src/engrid" >> $BINPREFIX/$START_ENGRID

  chmod 755 $BINPREFIX/$START_ENGRID
}

install_QT()
{
  echo "Install QT"
  if [ $DOWNLOAD_QT = 1 ]; then wget $URL_QT; fi
  tar -xzvf ./$ARCHIVE_QT
  cd $(basename $ARCHIVE_QT .tar.gz)
  mkdir -p $QTPREFIX
#   echo yes | ./configure --prefix=$QTPREFIX -opensource
  echo yes | ./configure --prefix=$QTPREFIX -opensource -nomake examples -nomake demos -nomake docs -no-webkit -no-phonon -no-phonon-backend -no-qt3support -no-accessibility -silent

  make $MAKEOPTIONS_ALL && make install
  cd -
}

install_VTK()
{
  echo "Install VTK"
  if [ $DOWNLOAD_VTK = 1 ]; then wget $URL_VTK; fi
  tar -xzvf ./$ARCHIVE_VTK
  cd ./VTK
  mkdir -p $VTKPREFIX
  cmake -DCMAKE_INSTALL_PREFIX:PATH=$VTKPREFIX -DBUILD_SHARED_LIBS:BOOL=ON -DVTK_USE_GUISUPPORT:BOOL=ON -DVTK_USE_QVTK:BOOL=ON -DVTK_USE_QT:BOOL=ON -DDESIRED_QT_VERSION:STRING=4  .
  chmod 644 Utilities/vtktiff/tif_fax3sm.c
  make $MAKEOPTIONS_ALL && make install
  cd -
}

install_CGNS()
{
  echo "Install CGNS"
  if [ $DOWNLOAD_CGNS = 1 ]; then wget $URL_CGNS; fi
  tar -xzvf ./$ARCHIVE_CGNS
  cd ./cgnslib_2.5/
  mkdir -p $CGNSPREFIX
  mkdir -p $CGNSPREFIX/include
  mkdir -p $CGNSPREFIX/lib
  ./configure --prefix=$CGNSPREFIX && make $MAKEOPTIONS_ALL && make install
  cd -
}

build_engrid()
{
  ORIG_WD=$(pwd)

  mkdir -p $SRCPREFIX
  cd $SRCPREFIX

  if [ $DOWNLOAD_ENGRID = 1 ];
  then
    wget $URL_ENGRID;
    tar -xzvf ./$ARCHIVE_ENGRID
  else
    git clone $GIT_URL_ENGRID
    cd $SRCPREFIX/engrid/src
    if [ $BRANCH != "master" ]; then git checkout -b $BRANCH origin/$BRANCH; fi;
  fi

  cd $SRCPREFIX/engrid/src

  echo "Build netgen"

  if [ $DOWNLOAD_NETGEN = 1 ];
  then
    cd $SRCPREFIX/engrid/src/netgen_svn
    wget $URL_NETGEN;
    tar -xzvf ./$ARCHIVE_NETGEN
    qmake && make $MAKEOPTIONS_ALL
    cd $SRCPREFIX/engrid/src
  else
    ./scripts/build-nglib.sh
  fi

  echo "Build enGrid"
  qmake $PROJECT_FILE && make $MAKEOPTIONS_ALL $MAKEOPTIONS_ENGRID

  cd $ORIG_WD
}

update_netgen()
{
  ORIG_WD=$(pwd)

  cd $SRCPREFIX/engrid/src
  echo "Update netgen"

  if [ $DOWNLOAD_NETGEN = 1 ];
  then
    cd $SRCPREFIX/engrid/src/netgen_svn
    wget $URL_NETGEN;
    tar -xzvf ./$ARCHIVE_NETGEN
    qmake && make $MAKEOPTIONS_ALL
    cd $SRCPREFIX/engrid/src
  else
    ./scripts/build-nglib.sh
  fi

  cd $ORIG_WD
}

update_engrid()
{
  ORIG_WD=$(pwd)

  cd $SRCPREFIX/engrid/src
  echo "Update enGrid"

  if [ $DOWNLOAD_ENGRID = 1 ];
  then
    cd $SRCPREFIX
    wget $URL_ENGRID;
    tar -xzvf ./$ARCHIVE_ENGRID
    cd $SRCPREFIX/engrid/src
  else
    git pull
  fi

  qmake $PROJECT_FILE && make $MAKEOPTIONS_ALL $MAKEOPTIONS_ENGRID

  cd $ORIG_WD
}

rebuild_engrid()
{
  cd $SRCPREFIX/engrid/src
  qmake && make distclean && qmake $PROJECT_FILE && make $MAKEOPTIONS_ALL $MAKEOPTIONS_ENGRID
  cd -
}

ans=$(zenity  --height=350 --list  --text "Which actions should be executed?" --checklist  --column "Run" --column "Actions" \
FALSE "create_bash_engrid" \
FALSE "install_QT" \
FALSE "install_VTK" \
FALSE "install_CGNS" \
FALSE "build_engrid" \
FALSE "update_netgen" \
FALSE "update_engrid" \
FALSE "rebuild_engrid" \
FALSE "create_start_engrid" \
--separator=":");
echo $ans

if ( echo $ans | grep -w "create_bash_engrid" ) then create_bash_engrid; fi

set +u
source $BINPREFIX/$ENV_SETUP
set -u

if ( echo $ans | grep -w "install_QT" ) then install_QT; fi;
if ( echo $ans | grep -w "install_VTK" ) then install_VTK; fi;
if ( echo $ans | grep -w "install_CGNS" ) then install_CGNS; fi;
if ( echo $ans | grep -w "build_engrid" ) then build_engrid; fi;
if ( echo $ans | grep -w "update_netgen" ) then update_netgen; fi;
if ( echo $ans | grep -w "update_engrid" ) then update_engrid; fi;
if ( echo $ans | grep -w "rebuild_engrid" ) then rebuild_engrid; fi;
if ( echo $ans | grep -w "create_start_engrid" ) then create_start_engrid; fi;

echo "SUCCESS"
exit 0
