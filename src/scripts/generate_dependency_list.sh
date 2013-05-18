#!/bin/bash
# 
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +                                                                      +
# + This file is part of enGrid.                                         +
# +                                                                      +
# + Copyright 2008-2013 enGits GmbH                                      +
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
set -ux

# Check if all parameters are present
# If no, exit
if [ $# -ne 4 ]
then
        echo "usage :"
        echo "`basename $0` INPUT OUTPUT INSTALL UNINSTALL"
        exit 0
fi

INPUT=$1
OUTPUT=$2
INSTALL=$3
UNINSTALL=$4

grep DLL $1 > $OUTPUT
sed -i '/windows/d' $OUTPUT
sed -i '/Warning/d' $OUTPUT
sed -i '/KnownDLLs/d' $OUTPUT
sed -i '/IESHIMS\.DLL/d' $OUTPUT
sed -i '/WER\.DLL/d' $OUTPUT
sed -i 's/DLL\(.*\)$/DLL/g' $OUTPUT
sed -i 's/\[.*\]//g' $OUTPUT
sed -i 's/^ *//g' $OUTPUT

~/bin/remove_duplicates.pl $OUTPUT
OUTPUT_2=no_dupes_$OUTPUT
sort $OUTPUT_2 > $OUTPUT

cp $OUTPUT $INSTALL
sed -i 's/c:\\qt\\2009\.03\\mingw\\bin/${MINGWBINDIR}/g' $INSTALL
sed -i 's/c:\\qt\\2009\.03\\qt\\bin/${QTBINDIR}/g' $INSTALL
sed -i 's/c:\\qwt-5\.3\.0-svn\\lib/${QWTLIBDIR}/g' $INSTALL
sed -i 's/c:\\qwtplot3d\\lib/${QWTPLOT3DLIBDIR}/g' $INSTALL
sed -i 's/c:\\vtk-5\.4\.2-install\\bin/${VTKBINDIR}/g' $INSTALL
sed -i 's/c:\\libraries\\vtk-5\.4\.2-install\\bin/${VTKBINDIR}/g' $INSTALL

sed -i 's/^/  File "/' $INSTALL
sed -i 's/$/"/' $INSTALL

cp $OUTPUT $UNINSTALL
sed -i 's/.*\\\([^\\]*\)\.DLL/\1.DLL/g' $UNINSTALL
sed -i 's/^/  Delete "$INSTDIR\\/' $UNINSTALL
sed -i 's/$/"/' $UNINSTALL

echo "=== OUTPUT ==="
cat $OUTPUT

echo "=== INSTALL ==="
cat $INSTALL

echo "=== UNINSTALL ==="
cat $UNINSTALL
