include (engrid.pri)

TEMPLATE = app
LANGUAGE = C++
#TARGET   = engrid

# Enable this if the VTK from the ParaView sources and 
# installation want to be used
# Note: Currently only for Windows Compiles with MSVC
Use_VTK_Win_ParaView = yes

CONFIG += qt \
          debug_and_release \
          thread
                  
QT     += xml \
          network \
          opengl

win32-msvc* {
  QMAKE_CXXFLAGS  += -W3
  DEFINES         += LIBENGRID_EXPORTS
  DEFINES         += DLL_EXPORT
  DEFINES         += _USE_MATH_DEFINES
  INCLUDEPATH     += ../../VTK/include/vtk-5.10
  LIBS            += -L../../VTK/lib/vtk-5.10
  LIBS            += -L../build-engrid-Desktop-Release/netgen_svn/release -lnglib
  LIBS            += -L../build-engrid-Desktop-Release/libengrid/release -lengrid
  brlcad {
    INCLUDEPATH += ../../BRL-CAD/include
    INCLUDEPATH += ../../BRL-CAD/include/openNURBS
    LIBS        += ../../BRL-CAD/lib/librt.lib
    LIBS        += ../../BRL-CAD/lib/libbu.lib
    DEFINES     += BRLCAD_SUPPORT
  }
} else {
  INCLUDEPATH     += $(VTKINCDIR)
  LIBS            += -L./netgen_svn -lng
  LIBS            += -L./libengrid -lengrid
  brlcad {
    INCLUDEPATH += $(BRLCADINCDIR)
    INCLUDEPATH += $(BRLCADINCDIR)/openNURBS
    LIBS        += $(BRLCADLIBDIR)/librt.so
    DEFINES     += BRLCAD_SUPPORT
  }
}

LIBS += -lQVTK
LIBS += -lvtkCommon
LIBS += -lvtkDICOMParser
LIBS += -lvtkexoIIc
LIBS += -lvtkFiltering
LIBS += -lvtkftgl
LIBS += -lvtkGenericFiltering
LIBS += -lvtkGraphics
LIBS += -lvtkHybrid
LIBS += -lvtkImaging
LIBS += -lvtkIO
LIBS += -lvtkRendering
LIBS += -lvtksys
LIBS += -lvtkVolumeRendering
LIBS += -lvtkWidgets

INCLUDEPATH += ./libengrid
INCLUDEPATH += ./libengrid-build
INCLUDEPATH += ../engrid-build
INCLUDEPATH += ./netgen_svn/netgen-mesher/netgen/nglib
INCLUDEPATH += ./netgen_svn/netgen-mesher/netgen/libsrc/general

OTHER_FILES += checkcomments.py todo.txt
RESOURCES   += libengrid/engrid.qrc

SOURCES = main.cpp


