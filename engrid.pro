TEMPLATE = app
LANGUAGE = C++

#CONFIG += qt release thread
CONFIG += qt debug thread
QT += xml network opengl

LIBS += -lvtkCommon 
LIBS += -lvtkGraphics 
LIBS += -lvtkImaging 
LIBS += -lvtkHybrid 
LIBS += -lQVTK


LIBS += -lng

#DEFINES += QT_NO_DEBUG

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
 vtkpolydataalgorithm2.h \
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
 pick_cell_point.h \
 guipick.h

SOURCES = \
main.cpp \
\
boundarycondition.cpp \
celllayeriterator.cpp \
cellneighbouriterator.cpp \
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
 vtkpolydataalgorithm2.cpp \
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
 pick_cell_point.cpp \
 guipick.cpp

FORMS = \
guicreateboundarylayer.ui \
guideletebadaspecttris.ui \
guidivideboundarylayer.ui \
guieditboundaryconditions.ui \
guimainwindow.ui \
guiimproveaspectratio.ui \
guinormalextrusion.ui \
guioutputwindow.ui \
guiselectboundarycodes.ui \
guisetboundarycode.ui \
guismoothsurface.ui \
 \
 guitransform.ui \
 guipick.ui


SOURCES -= settingstab.cpp \
settingsviewer.cpp \
 vtkpolydataalgorithm2.cpp \
 pick_cell_point.cpp
HEADERS -= settingstab.h \
settingsviewer.h \
 vtkpolydataalgorithm2.h \
 pick_cell_point.h
