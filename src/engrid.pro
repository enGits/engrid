TEMPLATE  = subdirs
LANGUAGE  = C++
CONFIG   += ordered

CONFIG += debug_and_release

SUBDIRS   = netgen
SUBDIRS  += libengrid
SUBDIRS  += engrid

netgen.file = netgen_svn/ng.pro

libengrid.file    = libengrid/libengrid.pro
libengrid.depends = netgen

engrid.file    = engrid_app.pro
engrid.depends = libengrid 
