# settings for enabling cgns support
DEFINES += CGNS_SUPPORT
!win32:!debian { 
    LIBS += -L$(CGNSLIBDIR)
    INCLUDEPATH += $(CGNSINCDIR)
}
LIBS += -lcgns
debian:LIBS += -lhdf5
