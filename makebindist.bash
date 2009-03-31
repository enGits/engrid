#!/bin/bash
export version=$1
export tarname="enGrid_linux64bit_"$version".tar"
export gzname="enGrid_linux64bit_"$version".tar.gz"
export dirname="enGrid_linux64bit_"$version
tar cf $tarname src/setup
tar -f $tarname -r src/start_engrid
tar -f $tarname -r src/enGrid_bin.tar.gz
cd src
qmake
make clean
cd netgen_svn
qmake
make clean
cd ../..
tar -f $tarname -r src/*.h
tar -f $tarname -r src/*.cpp
tar -f $tarname -r src/*.cxx
tar -f $tarname -r src/*.ui
tar -f $tarname -r src/licence.txt
tar -f $tarname -r src/license.txt
tar -f $tarname -r src/resources
tar -f $tarname -r src/engrid.pro
tar -f $tarname -r src/engrid.qrc
tar -f $tarname -r src/math/*.h
tar -f $tarname -r src/netgen_svn/netgen-mesher/netgen
tar -f $tarname -r src/netgen_svn/ng.pro
mkdir tmp
cd tmp
tar xf ../$tarname
mv src $dirname
tar czf $gzname $dirname
mv $gzname ..
cd ..
rm -rf tmp
rm $tarname


