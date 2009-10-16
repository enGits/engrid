TEMPLATE = app
LANGUAGE = C++
TARGET = engrid

# install
target.path = /usr/bin

# target.path = $$PREFIX/bin
INSTALLS += target

# CONFIG += qt release thread
# CONFIG += qt debug thread
CONFIG += qt \
    debug_and_release \
    thread

# DEFINES += QT_NO_DEBUG
# DEFINES += QT_DEBUG
# QMAKE_CXXFLAGS += -DAPP_VERSION=\\\"`date +'\"%a_%b_%d,_%Y\"'`\\\"
# ###############################
# VERSION INFO
# get "git revision number"
win32:QMAKE_CXXFLAGS += -DENGRID_VERSION=\\\"alpha\\\"
else:QMAKE_CXXFLAGS += -DENGRID_VERSION=\\\"`git \
    describe`\\\"

# ###############################
# to get rid of deprecated header warnings caused by including QVTKwidget.h
# DEFINES += VTK_EXCLUDE_STRSTREAM_HEADERS
# DEFINES += VTK_LEGACY_REMOVE
QMAKE_CXXFLAGS += -Wall

# QMAKE_CXXFLAGS += -pg
# QMAKE_CXXFLAGS += -O3
# QMAKE_LFLAGS += -pg
QT += xml \
    network \
    opengl

# #######################
# VTK
INCLUDEPATH += $(VTKINCDIR)
LIBS += -L$(VTKLIBDIR)

# #######################
# #######################
# NETGEN
INCLUDEPATH += ./netgen_svn/netgen-mesher/netgen/nglib
INCLUDEPATH += ./netgen_svn/netgen-mesher/netgen/libsrc/general

# #######################
include(engrid-standard.pri)
