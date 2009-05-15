#
# settings for enabling cgns support
#

DEFINES += CGNS_SUPPORT

!win32 {
    LIBS += -L$(CGNSLIBDIR)
    INCLUDEPATH += $(CGNSINCDIR)
}

LIBS    += -lcgns

#
# end
#
