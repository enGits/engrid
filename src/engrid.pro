TEMPLATE  = subdirs
LANGUAGE  = C++
CONFIG   += ordered recursive

CONFIG += debug_and_release

SUBDIRS   = netgen
SUBDIRS  += libengrid
SUBDIRS  += engrid

netgen.file = netgen_svn/ng.pro
libengrid.file    = libengrid/libengrid.pro
libengrid.depends = netgen

engrid.file    = engrid.pro.app
engrid.depends = libengrid 
