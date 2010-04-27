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
DESTDIR=$2

ORIG=`pwd`

DATE=$(date +%Y%m%d_%H%M%S)
ARCHIVE=$DESTDIR/$BASE\_$DATE.tar.gz

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
