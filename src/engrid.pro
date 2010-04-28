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
