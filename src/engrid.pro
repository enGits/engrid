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
# get "git revision number"
QMAKE_CXXFLAGS += -DENGRID_VERSION=\\\"`git \
    describe`\\\"
QMAKE_CXXFLAGS += -Wall
QMAKE_CXXFLAGS += -pg
QMAKE_LFLAGS += -pg
QT += xml \
    network \
    opengl
!win32 {
    LIBS += -L./netgen_svn
    LIBS += -L$(VTKLIBDIR)

    # LIBS += -Wl,-rpath
    QMAKE_CXXFLAGS += -Wno-deprecated
    INCLUDEPATH += $(VTKINCDIR)
    INCLUDEPATH += ./netgen_svn/netgen-mesher/netgen/nglib
    INCLUDEPATH += ./netgen_svn/netgen-mesher/netgen/libsrc/general
}
win32 {
    VTK_DIR = C:\VTK
    VTK_SRCDIR = C:\VTK\5.0.4
    LIBS += -L$$VTK_DIR\bin\release
    LIBS += -lvtkRendering
    LIBS += -lvtkFiltering
    LIBS += -lvtkIO
    LIBS += -lvtkfreetype
    LIBS += -lvtkftgl
    LIBS += -lvtkexpat
    LIBS += -lvtkzlib
    INCLUDEPATH += $$VTK_SRCDIR\COMMON
    INCLUDEPATH += $$VTK_SRCDIR\FILTER~1
    INCLUDEPATH += $$VTK_SRCDIR\GUISUP~1\QT
    INCLUDEPATH += $$VTK_SRCDIR\GENERI~1
    INCLUDEPATH += $$VTK_SRCDIR\GRAPHICS
    INCLUDEPATH += $$VTK_SRCDIR\HYBRID
    INCLUDEPATH += $$VTK_SRCDIR\IMAGING
    INCLUDEPATH += $$VTK_SRCDIR\IO
    INCLUDEPATH += $$VTK_SRCDIR\RENDER~1
    INCLUDEPATH += $$VTK_DIR
    INCLUDEPATH += netgen_svn\netgen-mesher\netgen\nglib
    INCLUDEPATH += netgen_svn\netgen-mesher\netgen\libsrc\general
    LIBS += -Lnetgen_svn\release
    DEFINES += _USE_MATH_DEFINES
}
include(engrid-standard.pri)
