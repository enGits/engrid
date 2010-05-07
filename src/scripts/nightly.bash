#!/usr/bin/env bash
# 
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +                                                                      +
# + This file is part of enGrid.                                         +
# +                                                                      +
# + Copyright 2008-2010 enGits GmbH                                     +
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

set -x

RECIPIENTS='mtaverne@engits.com ogloth@engits.com'

NIGHTLYDIR=/var/www/ftp/nightly

#Create a nightly source tarball and put it on the FTP server
./scripts/makedist.bash .. $NIGHTLYDIR

#Update online documentation
/usr/bin/doxygen Doxyfile

#Generate TODO lists
./scripts/checkcomments.py *.h *.cxx *.cpp math/*.h > comments.mail
if [ -s comments.mail ]
then
	mailx -s "ENGRID: comments" $RECIPIENTS < comments.mail
fi

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

# copy nightly build into nightly build directory
DATE=$(date +%Y%m%d_%H%M%S)
cp -v ./engrid "$NIGHTLYDIR/engrid_$DATE" || (echo mailx -s "failed to copy engrid into $NIGHTLYDIR" $RECIPIENTS)
