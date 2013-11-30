include (engrid.pri)

TEMPLATE = app
LANGUAGE = C++
TARGET   = engrid

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

#first: the essential libraries
win32-msvc* {
  QMAKE_CXXFLAGS  += -W3
  DEFINES         += _USE_MATH_DEFINES
  INCLUDEPATH     += $(VTKINCDIR)
  LIBS            += -L$(VTKLIBDIR)
  LIBS            += ./netgen_svn/release/nglib.lib
  LIBS            += ./libengrid/release/engrid.lib
  LIBS            += -lvtkzlib
  LIBS            += -lvtkexpat
  LIBS            += -lvtkfreetype
  LIBS            += -lvtkjpeg
  LIBS            += -lvtkpng
  LIBS            += -lvtktiff
} else {
  mydebug {
    QMAKE_CXXFLAGS  += -Wno-deprecated -g
  }

  INCLUDEPATH     += $(VTKINCDIR)
  
  LIBS            += -L./libengrid -lengrid

  debian {
    LIBS            += -lnglib
  } else {
    LIBS            += -L./netgen_svn -lng
  }
}

#second: the VTK libraries
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

#third: the 3rd party libraries
win32-msvc* {
  netcdf {
    DEFINES     += TAU_SUPPORT
    INCLUDEPATH += ../../netCDF/include
    LIBS        += ../../netCDF/lib/netcdf.lib
    LIBS        += ../../netCDF/lib/netcdfcxx.lib
  }
  brlcad {
    INCLUDEPATH += ../../BRL-CAD/include
    INCLUDEPATH += ../../BRL-CAD/include/openNURBS
    LIBS        += ../../BRL-CAD/lib/librt.lib
    LIBS        += ../../BRL-CAD/lib/libbu.lib
    DEFINES     += BRLCAD_SUPPORT
  }
} else {
  cgns {
    LIBS        += -lcgns
  }
  brlcad {
    INCLUDEPATH += $(BRLCADINCDIR)
    INCLUDEPATH += $(BRLCADINCDIR)/openNURBS
    DEFINES     += BRLCAD_SUPPORT

    #This is a necessary hack, so that this library does not conflict
    #with Qt's dependencies on the system's own librt.so
    LIBS += -lz $(BRLCADLIBDIR)/librt.so
  }
}

INCLUDEPATH += ./libengrid
INCLUDEPATH += ./libengrid-build
INCLUDEPATH += ../engrid-build
INCLUDEPATH += ./netgen_svn/netgen-mesher/netgen/nglib
INCLUDEPATH += ./netgen_svn/netgen-mesher/netgen/libsrc/general

OTHER_FILES += checkcomments.py todo.txt
RESOURCES   += libengrid/engrid.qrc

SOURCES = main.cpp
