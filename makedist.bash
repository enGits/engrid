#!/bin/bash
cd src
qmake
make clean
cd ..
export version=$1
export tarname="enGrid_"$version".tar"
export gzname="enGrid_"$version".tar.gz"
export dirname="enGrid_"$version
tar cf $tarname src/*.h
tar -f $tarname -r src/*.cpp
tar -f $tarname -r src/*.cxx
tar -f $tarname -r src/*.ui
tar -f $tarname -r src/licence.txt
tar -f $tarname -r src/license.txt
tar -f $tarname -r src/resources/icons
tar -f $tarname -r src/resources/kde_icons
tar -f $tarname -r src/engrid.pro
tar -f $tarname -r src/engrid.qrc
tar -f $tarname -r src/math/*.h
tar -f $tarname -r src/netgen_svn/ng.pro
tar -f $tarname -r src/netgen_svn/config.h
tar -f $tarname -r src/build-nglib.sh
mkdir tmp
cd tmp
tar xf ../$tarname
mv src $dirname
tar czf $gzname $dirname
mv $gzname ..
cd ..
rm -rf tmp
rm $tarname


