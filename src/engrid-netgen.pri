# #######################
# NETGEN
# #######################

!win32 {
    INCLUDEPATH += ./netgen_svn/netgen-mesher/netgen/nglib
    INCLUDEPATH += ./netgen_svn/netgen-mesher/netgen/libsrc/general
    LIBS += -L./netgen_svn
    LIBS += -lng
}

win32 {
    INCLUDEPATH += ./netgen_svn/netgen-mesher/netgen/nglib
    INCLUDEPATH += ./netgen_svn/netgen-mesher/netgen/libsrc/general
    LIBS += -Lnetgen_svn/release
    LIBS += -lng
}
