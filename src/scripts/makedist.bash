#!/usr/bin/env bash
#
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +                                                                      +
# + This file is part of enGrid.                                         +
# +                                                                      +
# + Copyright 2008,2009 Oliver Gloth                                     +
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
# + along with enGrid. If not, see <http:#www.gnu.org/licenses/>.        +
# +                                                                      +
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# DESCRIPTION:
# This script creates a simple source tarball from the current git repository.

# set -x

# Check if all parameters are present
# If no, exit
if [ $# -ne 2 ]
then
	echo
	echo 'archives a git repository as <dir>_$(date +%Y%m%d_%H%M%S).tar.gz'
	echo "usage :"
	echo "$0 <GITDIR> <DESTDIR>"
        echo "GITDIR = directory to archive. Usually 'engrid'"
        echo "DESTDIR = Where the archive should be placed."
        echo
        exit 0
fi

DIR=$(readlink -f $1)
BASE="engrid"
DESTDIR=$(readlink -f $2)

ORIG=`pwd`

DATE=$(date +%Y%m%d_%H%M%S)
ARCHIVE=$DESTDIR/$BASE\_$DATE.tar.gz

# exit 0

cd $DIR
if [ $? -ne 0 ]
then
	echo "ERROR: Could not change directory"
	exit 2
fi

git status
if [ $? -eq 128 ]
then
	echo "ERROR: No git repository found."
	exit 2
fi

if [ $( git status | wc -l ) -ne 2 ]
then
	echo "Please commit latest changes before archiving. ;)"
        echo "If you want to force archival, use:"
        echo " cd $DIR"
        echo " git archive --format=tar --prefix=$BASE/ HEAD | gzip >$ARCHIVE"
	exit 2
fi

#hg archive -v -t"tgz" $ARCHIVE
git archive --format=tar --prefix=$BASE/ HEAD | gzip >$ARCHIVE

if [ $? -eq 0 ]
then
	echo "All archived in $ARCHIVE"
	exit 0
else
	echo "ERROR: Archiving failed"
	exit 1
fi

cd $ORIG

#####################################################
#deprecated code used when CVS was still used:
# cd src
# qmake
# make clean
# cd ..
# export version=$1
# export tarname="enGrid_"$version".tar"
# export gzname="enGrid_"$version".tar.gz"
# export dirname="enGrid_"$version
# tar cf $tarname src/*.h
# tar -f $tarname -r src/*.cpp
# tar -f $tarname -r src/*.cxx
# tar -f $tarname -r src/*.ui
# tar -f $tarname -r src/licence.txt
# tar -f $tarname -r src/license.txt
# tar -f $tarname -r src/resources/icons
# tar -f $tarname -r src/resources/kde_icons
# tar -f $tarname -r src/engrid.pro
# tar -f $tarname -r src/engrid.qrc
# tar -f $tarname -r src/math/*.h
# tar -f $tarname -r src/netgen_svn/ng.pro
# tar -f $tarname -r src/netgen_svn/config.h
# tar -f $tarname -r src/build-nglib.sh
# mkdir tmp
# cd tmp
# tar xf ../$tarname
# mv src $dirname
# tar czf $gzname $dirname
# mv $gzname ..
# cd ..
# rm -rf tmp
# rm $tarname
