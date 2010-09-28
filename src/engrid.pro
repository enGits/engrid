TEMPLATE = app
LANGUAGE = C++
TARGET   = engrid
SUBDIRS  = libengrid

CONFIG         += qt debug_and_release thread
QT             += xml network opengl
QMAKE_CXXFLAGS += -Wall
QMAKE_CXXFLAGS += -Wno-deprecated
QMAKE_CXXFLAGS += -DGIT_VERSION=\\\"`git describe`"\\\"

INCLUDEPATH += ./libengrid
INCLUDEPATH += ./netgen_svn/netgen-mesher/netgen/nglib
INCLUDEPATH += ./netgen_svn/netgen-mesher/netgen/libsrc/general

LIBS += -lm
LIBS += -L./netgen_svn -lng
LIBS += -L./libengrid -lengrid

INCLUDEPATH += $(VTKINCDIR)
LIBS        += -L$(VTKLIBDIR)
LIBS        += -lQVTK
LIBS        += -lvtkCommon
LIBS        += -lvtkDICOMParser
LIBS        += -lvtkexoIIc
LIBS        += -lvtkFiltering
LIBS        += -lvtkftgl
LIBS        += -lvtkGenericFiltering
LIBS        += -lvtkGraphics
LIBS        += -lvtkHybrid
LIBS        += -lvtkImaging
LIBS        += -lvtkIO
LIBS        += -lvtkNetCDF
LIBS        += -lvtkRendering
LIBS        += -lvtksys
LIBS        += -lvtkVolumeRendering
LIBS        += -lvtkWidgets

OTHER_FILES += checkcomments.py todo.txt
RESOURCES   += engrid.qrc

SOURCES = main.cpp 
