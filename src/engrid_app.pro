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
} else {
    QMAKE_CXXFLAGS += -Wall
    QMAKE_CXXFLAGS += -Wno-deprecated
}

INCLUDEPATH += netgen_svn/netgen-mesher/netgen/nglib
INCLUDEPATH += netgen_svn/netgen-mesher/netgen/libsrc/general

INCLUDEPATH += libengrid
INCLUDEPATH += $${OUT_PWD}/libengrid

win32-msvc* {
    LIBS += -L./libengrid/release -lengrid
} else {
    LIBS += -lm
    LIBS += -L./libengrid -lengrid
}

win32-msvc* {
    DEFINES += _USE_MATH_DEFINES
        
    !isEmpty(Use_VTK_Win_ParaView) {
        include(misc/engrid-vtk-win_paraview.pri)
    } else {
        INCLUDEPATH += $(VTKINCDIR)
    }
} else {
    INCLUDEPATH += $(VTKINCDIR)
}


OTHER_FILES += checkcomments.py todo.txt
RESOURCES   += engrid.qrc

SOURCES = main.cpp

 
