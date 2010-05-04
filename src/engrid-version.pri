# ###############################
# VERSION INFO
# get "git revision number"
ENGRID_VERSION = \\\"1.1-pre-release\\\"
win32:GIT_DESCRIBE = \\\"\\\"
else:GIT_DESCRIBE = \\\"$$system(git describe)\\\"
message(GIT_DESCRIBE : $${GIT_DESCRIBE} )
message(ENGRID_VERSION : $${ENGRID_VERSION})
message(Qt version : $$[QT_VERSION])
QMAKE_CXXFLAGS += -DENGRID_VERSION=$${ENGRID_VERSION} \
    -DGIT_DESCRIBE=$${GIT_DESCRIBE}
message(QMAKE_CXXFLAGS : $${QMAKE_CXXFLAGS} )

# deprecated date/time define
# QMAKE_CXXFLAGS += -DAPP_VERSION=\\\"`date +'\"%a_%b_%d,_%Y\"'`\\\"
