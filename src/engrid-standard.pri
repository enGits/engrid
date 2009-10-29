# common libraries, includes, source files for the engrid*.pro files
QT += xml \
    network \
    opengl
RESOURCES += engrid.qrc

################
# netgen lib
################
!win32 {
    LIBS += -L./netgen_svn
    LIBS += -lng
#   LIBS += -Wl,-rpath
    QMAKE_CXXFLAGS += -Wno-deprecated
}

win32 {
	INCLUDEPATH += ./netgen_svn/netgen-mesher/netgen/nglib
	INCLUDEPATH += ./netgen_svn/netgen-mesher/netgen/libsrc/general
        LIBS += -Lnetgen_svn/release
        LIBS += -lsuperman

	#Z:\mtaverne\Development\engrid\src\netgen_svn\release
}
################

# VTK libs
LIBS += -lQVTK
LIBS += -lvtkCommon
LIBS += -lvtkDICOMParser
LIBS += -lvtkexoIIc

# LIBS += -lvtkexpat
LIBS += -lvtkFiltering

# LIBS += -lvtkfreetype
LIBS += -lvtkftgl
LIBS += -lvtkGenericFiltering
LIBS += -lvtkGraphics
LIBS += -lvtkHybrid
LIBS += -lvtkImaging

# LIBS += -lvtkInfovis
LIBS += -lvtkIO

# LIBS += -lvtkjpeg
# LIBS += -lvtklibxml2
# LIBS += -lvtkmetaio
LIBS += -lvtkNetCDF

# LIBS += -lvtkpng
LIBS += -lvtkRendering

# LIBS += -lvtksqlite
LIBS += -lvtksys

# LIBS += -lvtktiff
# LIBS += -lvtkViews
LIBS += -lvtkVolumeRendering
LIBS += -lvtkWidgets
HEADERS = boundarycondition.h \
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
    meshpartition.h \
    neutralwriter.h \
    nodelayeriterator.h \
    operation.h \
    optimisation.h \
    physicalboundarycondition.h \
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
    plywriter.h \
    uniquevector.h \
    swaptriangles.h \
    tvtkoperation.h \
    volumedefinition.h \
    vtkreader.h \
    vtkEgBoundaryCodesFilter.h \
    vtkEgEliminateShortEdges.h \
    vtkEgExtractVolumeCells.h \
    vtkEgGridFilter.h \
    vtkEgNormalExtrusion.h \
    vtkEgPolyDataToUnstructuredGridFilter.h \
    vtkImplicitPolyData.h \
    guicreateboundarylayer.h \
    guicreatevolumemesh.h \
    guideletebadaspecttris.h \
    guidivideboundarylayer.h \
    guieditboundaryconditions.h \
    guiimproveaspectratio.h \
    guimainwindow.h \
    guinormalextrusion.h \
    guiselectboundarycodes.h \
    guisetboundarycode.h \
    guicreatesurfacemesh.h \
    guisettingstab.h \
    guisettingsviewer.h \
    guivolumedelegate.h \
    guitransform.h \
    vertexdelegate.h \
    vertexmeshdensity.h \
    smoothingutilities.h \
    settingssheet.h \
    laplacesmoother.h \
    deletepickedpoint.h \
    text3d.h \
    guipick.h \
    egvtkinteractorstyle.h \
    insertpoints.h \
    removepoints.h \
    reducedpolydatareader.h \
    showinfo.h \
    surfacemesher.h \
    updatedesiredmeshdensity.h \
    boxselect.h \
    checksurfaceintegrity.h \
    surfaceoperation.h \
    surfaceprojection.h \
    octree.h \
    filetemplate.h \
    openfoamcase.h \
    multipagewidget.h \
    tricoord.h \
    updatesurfproj.h \
    foamobject.h \
    multipagewidgetpage.h \
    xmlhandler.h \
    openfoamtools.h

SOURCES = main.cpp \
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
    meshpartition.cpp \
    neutralwriter.cpp \
    nodelayeriterator.cpp \
    operation.cpp \
    optimisation.cpp \
    physicalboundarycondition.cpp \
    polydatareader.cpp \
    polymesh.cpp \
    seedsimpleprismaticlayer.cpp \
    setboundarycode.cpp \
    simplefoamwriter.cpp \
    stlreader.cpp \
    stlwriter.cpp \
    plywriter.cpp \
    swaptriangles.cpp \
    volumedefinition.cpp \
    vtkreader.cpp \
    vtkEgBoundaryCodesFilter.cxx \
    vtkEgEliminateShortEdges.cxx \
    vtkEgExtractVolumeCells.cxx \
    vtkEgGridFilter.cxx \
    vtkEgNormalExtrusion.cxx \
    vtkEgPolyDataToUnstructuredGridFilter.cxx \
    vtkImplicitPolyData.cpp \
    guicreateboundarylayer.cpp \
    guicreatevolumemesh.cpp \
    guideletebadaspecttris.cpp \
    guidivideboundarylayer.cpp \
    guieditboundaryconditions.cpp \
    guiimproveaspectratio.cpp \
    guimainwindow.cpp \
    guinormalextrusion.cpp \
    guiselectboundarycodes.cpp \
    guisetboundarycode.cpp \
    guicreatesurfacemesh.cpp \
    guisettingstab.cpp \
    guisettingsviewer.cpp \
    guivolumedelegate.cpp \
    guitransform.cpp \
    vertexdelegate.cpp \
    vertexmeshdensity.cpp \
    smoothingutilities.cpp \
    settingssheet.cpp \
    laplacesmoother.cpp \
    deletepickedpoint.cpp \
    text3d.cpp \
    guipick.cpp \
    egvtkinteractorstyle.cpp \
    insertpoints.cpp \
    removepoints.cpp \
    showinfo.cpp \
    surfacemesher.cpp \
    updatedesiredmeshdensity.cpp \
    boxselect.cpp \
    checksurfaceintegrity.cpp \
    surfaceoperation.cpp \
    surfaceprojection.cpp \
    octree.cpp \
    filetemplate.cpp \
    openfoamcase.cpp \
    multipagewidget.cpp \
    tricoord.cpp \
    updatesurfproj.cpp \
    foamobject.cpp \
    multipagewidgetpage.cpp \
    xmlhandler.cpp \
    reducedpolydatareader.cpp \
    openfoamtools.cpp

FORMS = guicreateboundarylayer.ui \
    guideletebadaspecttris.ui \
    guidivideboundarylayer.ui \
    guieditboundaryconditions.ui \
    guimainwindow.ui \
    guiimproveaspectratio.ui \
    guinormalextrusion.ui \
    guiselectboundarycodes.ui \
    guisetboundarycode.ui \
    guicreatesurfacemesh.ui \
    guitransform.ui \
    guipick.ui \
    guicreatevolumemesh.ui

OTHER_FILES += checkcomments.py \
    todo.txt

HEADERS += surfacealgorithm.h
SOURCES += surfacealgorithm.cpp
HEADERS += reducesurfacetriangulation.h
SOURCES += reducesurfacetriangulation.cpp
HEADERS += eliminatesmallbranches.h
SOURCES += eliminatesmallbranches.cpp
HEADERS += smoothandswapsurface.h
SOURCES += smoothandswapsurface.cpp
HEADERS += seligairfoilreader.h
SOURCES += seligairfoilreader.cpp
HEADERS += fixcadgeometry.h
SOURCES += fixcadgeometry.cpp
HEADERS += blenderreader.h
SOURCES += blenderreader.cpp
HEADERS += beziertriangle.h
SOURCES += beziertriangle.cpp
HEADERS += projection_test.h
SOURCES += projection_test.cpp
