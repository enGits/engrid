#!/bin/bash
set -eux

echo "Building engrid.pro release version"
make distclean && qmake engrid.pro && make -j2
echo "Building engrid.pro.cgns release version"
make distclean && qmake engrid.pro.cgns && make -j2

echo "Building engrid.pro debug version"
make distclean && qmake engrid.pro && make -j2 debug
echo "Building engrid.pro.cgns debug version"
make distclean && qmake engrid.pro.cgns && make -j2 debug

echo "SUCCESS: Everything compiles."
