TEMPLATE = app
LANGUAGE = C++
TARGET   = engrid

# CONFIG += qt release thread
# CONFIG += qt debug thread
CONFIG += qt \
    debug_and_release \
    thread

# DEFINES += QT_NO_DEBUG
# DEFINES += QT_DEBUG

include(engrid-version.pri)
include(engrid-standard.pri)

# install
target.path = /usr/bin
# target.path = $$PREFIX/bin
INSTALLS += target

# #######################
# VTK
INCLUDEPATH += $(VTKINCDIR)
LIBS += -L$(VTKLIBDIR)
# #######################
