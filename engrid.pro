TEMPLATE = app
LANGUAGE = C++
TARGET = engrid

# install
target.path = /usr/bin

# target.path = $$PREFIX/bin
INSTALLS += target
CONFIG += qt release thread
#CONFIG += qt debug thread
DEFINES += QT_NO_DEBUG

DEFINES += CGNS_SUPPORT
LIBS    += -lcgns

QT += xml network opengl

#VTK libs
LIBS += -lQVTK
LIBS += -lvtkCommon
LIBS += -lvtkDICOMParser
LIBS += -lvtkexoIIc
#LIBS += -lvtkexpat
LIBS += -lvtkFiltering
#LIBS += -lvtkfreetype
LIBS += -lvtkftgl
LIBS += -lvtkGenericFiltering
LIBS += -lvtkGraphics
LIBS += -lvtkHybrid
LIBS += -lvtkImaging
#LIBS += -lvtkInfovis
LIBS += -lvtkIO
#LIBS += -lvtkjpeg
#LIBS += -lvtklibxml2
#LIBS += -lvtkmetaio
LIBS += -lvtkNetCDF
#LIBS += -lvtkpng
LIBS += -lvtkRendering
#LIBS += -lvtksqlite
LIBS += -lvtksys
#LIBS += -lvtktiff
#LIBS += -lvtkViews
LIBS += -lvtkVolumeRendering
LIBS += -lvtkWidgets
#LIBS += -lvtkzlib

#netgen lib
LIBS += -lng

!win32 { 
    LIBS += -L./netgen_svn
    LIBS += -L$(VTKLIBDIR)
    LIBS += -Wl,-rpath
    QMAKE_CXXFLAGS += -Wno-deprecated
    INCLUDEPATH += $(VTKINCDIR)
    INCLUDEPATH += ./netgen_svn/netgen-mesher/netgen/nglib
    INCLUDEPATH += ./netgen_svn/netgen-mesher/netgen/libsrc/general
}
win32 { 
    VTK_DIR = C:\VTK
    VTK_SRCDIR = C:\VTK\5.0.4
    LIBS += -L$$VTK_DIR\bin\release
    LIBS += -lvtkRendering
    LIBS += -lvtkFiltering
    LIBS += -lvtkIO
    LIBS += -lvtkfreetype
    LIBS += -lvtkftgl
    LIBS += -lvtkexpat
    LIBS += -lvtkzlib
    INCLUDEPATH += $$VTK_SRCDIR\COMMON
    INCLUDEPATH += $$VTK_SRCDIR\FILTER~1
    INCLUDEPATH += $$VTK_SRCDIR\GUISUP~1\QT
    INCLUDEPATH += $$VTK_SRCDIR\GENERI~1
    INCLUDEPATH += $$VTK_SRCDIR\GRAPHICS
    INCLUDEPATH += $$VTK_SRCDIR\HYBRID
    INCLUDEPATH += $$VTK_SRCDIR\IMAGING
    INCLUDEPATH += $$VTK_SRCDIR\IO
    INCLUDEPATH += $$VTK_SRCDIR\RENDER~1
    INCLUDEPATH += $$VTK_DIR
    INCLUDEPATH += netgen_svn\netgen-mesher\netgen\nglib
    INCLUDEPATH += netgen_svn\netgen-mesher\netgen\libsrc\general
    LIBS += -Lnetgen_svn\release
    DEFINES += _USE_MATH_DEFINES
}
RESOURCES += engrid.qrc

INCLUDEPATH += ./
INCLUDEPATH += ./nglib

HEADERS = \
\
boundarycondition.h \
celllayeriterator.h \
cellneighbouriterator.h \
cgnswriter.h \
containertricks.h \
correctsurfaceorientation.h \
createvolumemesh.h \
deletecells.h \
deletetetras.h \
deletepickedcell.h \
deletevolumegrid.h \
dialogoperation.h \
egvtkobject.h \
elements.h \
engrid.h \
error.h \
fixstl.h \
foamreader.h \
foamwriter.h \
geometrytools.h \
gmshiooperation.h \
gmshreader.h \
gmshwriter.h \
gridsmoother.h \
iooperation.h \
iterator.h \
layeriterator.h \
neutralwriter.h \
nodelayeriterator.h \
operation.h \
optimisation.h \
polydatareader.h \
polymesh.h \
seedsimpleprismaticlayer.h \
setboundarycode.h \
simplefoamwriter.h \
sortablevector.h \
std_connections.h \
std_includes.h \
stlreader.h \
stlwriter.h \
uniquevector.h \
swaptriangles.h \
tvtkoperation.h \
vtkreader.h \
\
vtkEgBoundaryCodesFilter.h \
vtkEgEliminateShortEdges.h \
vtkEgExtractVolumeCells.h \
vtkEgGridFilter.h \
vtkEgNormalExtrusion.h \
vtkEgPolyDataToUnstructuredGridFilter.h \
\
guicreateboundarylayer.h \
guideletebadaspecttris.h \
guidivideboundarylayer.h \
guieditboundaryconditions.h \
guiimproveaspectratio.h \
guimainwindow.h \
guinormalextrusion.h \
guiselectboundarycodes.h \
guisetboundarycode.h \
guismoothsurface.h \
guisettingstab.h \
guisettingsviewer.h \
 \
 guitransform.h \
 vtksmoothpolydatafilter2.h \
 createspecialmapping.h \
 vertexdelegate.h \
 vertexmeshdensity.h \
 smoothingutilities.h \
 settingssheet.h \
 vtkeggridsmoother.h \
 vtkeggridsmoothpolydatafilter.h \
 vtkeggridwindowedsincpolydatafilter.h \
 laplacesmoother.h \
 deletepickedpoint.h \
 text3d.h \
 guipick.h \
 egvtkinteractorstyle.h

SOURCES = \
main.cpp \
\
boundarycondition.cpp \
celllayeriterator.cpp \
cellneighbouriterator.cpp \
cgnswriter.cpp \
correctsurfaceorientation.cpp \
createvolumemesh.cpp \
deletecells.cpp \
deletepickedcell.cpp \
deletetetras.cpp \
deletevolumegrid.cpp \
egvtkobject.cpp \
elements.cpp \
error.cpp \
fixstl.cpp \
foamreader.cpp \
foamwriter.cpp \
geometrytools.cpp \
gmshiooperation.cpp \
gmshreader.cpp \
gmshwriter.cpp \
gridsmoother.cpp \
iooperation.cpp \
iterator.cpp \
layeriterator.cpp \
neutralwriter.cpp \
nodelayeriterator.cpp \
operation.cpp \
optimisation.cpp \
polydatareader.cpp \
polymesh.cpp \
seedsimpleprismaticlayer.cpp \
setboundarycode.cpp \
simplefoamwriter.cpp \
stlreader.cpp \
stlwriter.cpp \
swaptriangles.cpp \
vtkreader.cpp \
\
vtkEgBoundaryCodesFilter.cxx \
vtkEgEliminateShortEdges.cxx \
vtkEgExtractVolumeCells.cxx \
vtkEgGridFilter.cxx \
vtkEgNormalExtrusion.cxx \
vtkEgPolyDataToUnstructuredGridFilter.cxx \
\
guicreateboundarylayer.cpp \
guideletebadaspecttris.cpp \
guidivideboundarylayer.cpp \
guieditboundaryconditions.cpp \
guiimproveaspectratio.cpp \
guimainwindow.cpp \
guinormalextrusion.cpp \
guiselectboundarycodes.cpp \
guisetboundarycode.cpp \
guismoothsurface.cpp \
guisettingstab.cpp \
guisettingsviewer.cpp \
 \
 guitransform.cpp \
 vtksmoothpolydatafilter2.cpp \
 createspecialmapping.cpp \
 vertexdelegate.cpp \
 vertexmeshdensity.cpp \
 smoothingutilities.cpp \
 settingssheet.cpp \
 vtkeggridsmoother.cpp \
 vtkeggridsmoothpolydatafilter.cpp \
 vtkeggridwindowedsincpolydatafilter.cpp \
 laplacesmoother.cpp \
 deletepickedpoint.cpp \
 text3d.cpp \
 guipick.cpp \
 egvtkinteractorstyle.cpp

FORMS = \
guicreateboundarylayer.ui \
guideletebadaspecttris.ui \
guidivideboundarylayer.ui \
guieditboundaryconditions.ui \
guimainwindow.ui \
guiimproveaspectratio.ui \
guinormalextrusion.ui \
guiselectboundarycodes.ui \
guisetboundarycode.ui \
guismoothsurface.ui \
guitransform.ui \
guipick.ui
