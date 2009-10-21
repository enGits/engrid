#!/bin/bash
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
# This is a script run every night for several tasks: Documentation update, todo list generation, nightly builds, etc
# USAGE:
# This script must be run from the "engrid/src" directory.

set -x

RECIPIENTS='mtaverne@engits.com ogloth@engits.com'

#Update online documentation
/usr/bin/doxygen Doxyfile

#Generate TODO lists
./checkcomments.py *.h *.cxx *.cpp math/*.h > comments.mail
mailx -s "ENGRID: comments" $RECIPIENTS < comments.mail

#test build
touch build.log
pwd
./scripts/rebuild.sh 1>build.log 2>&1
if [ $? -ne 0 ]
then
  echo "BUILD FAILED"
  mailx -s "ENGRID: build test failed" $RECIPIENTS < ./build.log
else
  echo "BUILD SUCCESSFUL"
  mailx -s "ENGRID: build test successful" $RECIPIENTS < ./build.log
fi

# copy nightly build into /opt/shared/bin/
cp -v ./engrid /opt/shared/bin/ || (echo mailx -s "failed to copy engrid into /opt/shared/bin/" $RECIPIENTS)
