#!/bin/bash
set -eux

echo "Building engrid.pro release version"
qmake && make distclean && qmake engrid.pro && make -j2
echo "Building engrid.pro.cgns release version"
qmake && make distclean && qmake engrid.pro.cgns && make -j2

echo "Building engrid.pro debug version"
qmake && make distclean && qmake engrid.pro && make -j2 debug
echo "Building engrid.pro.cgns debug version"
qmake && make distclean && qmake engrid.pro.cgns && make -j2 debug

echo "SUCCESS: Everything compiles."
