#!/bin/bash
/usr/bin/doxygen Doxyfile
./checkcomments.py > comments.mail
mailx -s "ENGRID: comments" ogloth@engits.com < comments.mail
rm -f comments.mail
