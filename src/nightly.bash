#!/bin/bash
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
