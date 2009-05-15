#!/bin/bash
set -eux

make distclean && qmake engrid.pro && make
make distclean && qmake engrid.pro.cgns && make

echo "SUCCESS: Both project files work."

