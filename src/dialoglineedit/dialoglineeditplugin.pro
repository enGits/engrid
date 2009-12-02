CONFIG      += designer plugin
TARGET      = $$qtLibraryTarget($$TARGET)
TEMPLATE    = lib
QTDIR_build:DESTDIR     = $$QT_BUILD_TREE/plugins/designer

#input
################################
#UTILITIES
INCLUDEPATH += ..
HEADERS += ../utilities.h ../error.h
SOURCES += ../utilities.cpp ../error.cpp
################################

HEADERS     = dialoglineedit.h \
              dialoglineeditplugin.h
SOURCES     = dialoglineedit.cpp \
              dialoglineeditplugin.cpp

# install
target.path = $$[QT_INSTALL_PLUGINS]/designer
sources.files = $$SOURCES $$HEADERS *.pro
sources.path = $$[QT_INSTALL_EXAMPLES]/designer/dialoglineeditplugin
INSTALLS += target sources
