TEMPLATE      = lib
CONFIG       += plugin
INCLUDEPATH  += ../..
INCLUDEPATH  += ../../libengrid
INCLUDEPATH  += $(VTKINCDIR)
HEADERS       = exampleplugin.h
SOURCES       = exampleplugin.cpp
FORMS         = exampleplugin.ui

QMAKE_CXXFLAGS += -Wno-deprecated
