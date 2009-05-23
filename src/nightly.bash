#!/bin/bash
/usr/bin/doxygen Doxyfile
./checkcomments.py *.h *.cxx *.cpp math/*.h > comments.mail
mailx -s "ENGRID: comments" ogloth@engits.com < comments.mail
mailx -s "ENGRID: comments" mtaverne@engits.com < comments.mail
