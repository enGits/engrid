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

win32-msvc* {
    QMAKE_CXXFLAGS += -W3
} win32-g++* {
    CONFIG += console
    QMAKE_CXXFLAGS += -Wall
    QMAKE_CXXFLAGS += -Wno-deprecated
    QMAKE_CXXFLAGS += -Wl,--no-undefined
    QMAKE_CXXFLAGS += -Wl,--enable-runtime-pseudo-reloc
} else {
    QMAKE_CXXFLAGS += -Wall
    QMAKE_CXXFLAGS += -Wno-deprecated
}

INCLUDEPATH += netgen_svn/netgen-mesher/netgen/nglib
INCLUDEPATH += netgen_svn/netgen-mesher/netgen/libsrc/general

#various paths for the same thing, due to some crazy bugs between versions on qt/qmake/gcc
INCLUDEPATH += libengrid
INCLUDEPATH += $${OUT_PWD}/libengrid
INCLUDEPATH += $${OUT_PWD}/.

win32-msvc* {
    LIBS += -L./libengrid/release -lengrid
} win32-g++* {
    LIBS += -L./libengrid/release -lengrid
} else {
    LIBS += -lm
    LIBS += -L./libengrid -lengrid
    LIBS += -L./netgen_svn -lng
}

win32-msvc* {
    DEFINES += _USE_MATH_DEFINES
        
    !isEmpty(Use_VTK_Win_ParaView) {
        include(misc/engrid-vtk-win_paraview.pri)
    } else {
        INCLUDEPATH += $(VTKINCDIR)
    }
} win32-g++* {
    INCLUDEPATH += $(VTKINCDIR)
} else {
    INCLUDEPATH += $(VTKINCDIR)
}


OTHER_FILES += checkcomments.py todo.txt
RESOURCES   += engrid.qrc

SOURCES = main.cpp

 
