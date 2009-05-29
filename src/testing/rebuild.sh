#!/bin/bash
set -eux

make distclean && qmake engrid.pro && make
make distclean && qmake engrid.pro.cgns && make

make distclean && qmake engrid.pro && make debug
make distclean && qmake engrid.pro.cgns && make debug

echo "SUCCESS: Both project files work."
