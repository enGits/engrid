#
# settings for enabling cgns support
#

DEFINES += CGNS_SUPPORT

!win32 {
    LIBS += -L$(CGNSLIBDIR)
}

LIBS    += -lcgns

#
# end
#
