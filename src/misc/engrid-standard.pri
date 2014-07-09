# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +                                                                      +
# + This file is part of enGrid.                                         +
# +                                                                      +
# + Copyright 2008-2014 enGits GmbH                                      +
# +                                                                      +
# + enGrid is free software: you can redistribute it and/or modify       +
# + it under the terms of the GNU General Public License as published by +
# + the Free Software Foundation, either version 3 of the License, or    +
# + (at your option) any later version.                                  +
# +                                                                      +
# + enGrid is distributed in the hope that it will be useful,            +
# + but WITHOUT ANY WARRANTY; without even the implied warranty of       +
# + MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        +
# + GNU General Public License for more details.                         +
# +                                                                      +
# + You should have received a copy of the GNU General Public License    +
# + along with enGrid. If not, see <http://www.gnu.org/licenses/>.       +
# +                                                                      +
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
TEMPLATE = app
LANGUAGE = C++
TARGET = engrid

CONFIG         += qt debug_and_release thread

QMAKE_CXXFLAGS += -Wall
QT             += xml network opengl
QMAKE_CXXFLAGS += -Wno-deprecated

LIBS += -lm

INCLUDEPATH += ./netgen_svn/netgen-mesher/netgen/nglib
INCLUDEPATH += ./netgen_svn/netgen-mesher/netgen/libsrc/general
LIBS        += -L./netgen_svn
LIBS        += -lng

INCLUDEPATH += $(VTKINCDIR)
LIBS        += -L$(VTKLIBDIR)
LIBS        += -lQVTK
LIBS        += -lvtkCommon
LIBS        += -lvtkDICOMParser
LIBS        += -lvtkexoIIc
LIBS        += -lvtkFiltering
LIBS        += -lvtkftgl
LIBS        += -lvtkGenericFiltering
LIBS        += -lvtkGraphics
LIBS        += -lvtkHybrid
LIBS        += -lvtkImaging
LIBS        += -lvtkIO
LIBS        += -lvtkNetCDF
LIBS        += -lvtkRendering
LIBS        += -lvtksys
LIBS        += -lvtkVolumeRendering
LIBS        += -lvtkWidgets

OTHER_FILES += checkcomments.py todo.txt
RESOURCES   += engrid.qrc

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
    openfoamtools.h \
    sharpenedges.h \
    checkforoverlap.h \
    timer.h \
    facefinder.h \
    math/linsolve.h \
    math/mathvector.h \
    math/mathvector_methods.h \
    math/mathvector_operators.h \
    math/mathvector_structs.h \
    math/smallsquarematrix.h \
    pointfinder.h
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
    openfoamtools.cpp \
    sharpenedges.cpp \
    checkforoverlap.cpp \
    timer.cpp \
    facefinder.cpp \
    pointfinder.cpp
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
HEADERS += blenderwriter.h
SOURCES += blenderwriter.cpp
HEADERS += dialoglineedit/dialoglineedit.h
SOURCES += dialoglineedit/dialoglineedit.cpp
HEADERS += utilities.h
SOURCES += utilities.cpp
HEADERS += edgelengthsourcemanager.h edgelengthsource.h
SOURCES += edgelengthsourcemanager.cpp
FORMS += guiedgelengthsourcesphere.ui
HEADERS += guiedgelengthsourcesphere.h
SOURCES += guiedgelengthsourcesphere.cpp
HEADERS += triangle.h
SOURCES += triangle.cpp
HEADERS += mergenodes.h
SOURCES += mergenodes.cpp
FORMS += guiedgelengthsourcecone.ui
HEADERS += guiedgelengthsourcecone.h
SOURCES += guiedgelengthsourcecone.cpp
FORMS += guimergevolumes.ui
HEADERS += guimergevolumes.h
SOURCES += guimergevolumes.cpp
HEADERS += deletestraynodes.h
SOURCES += deletestraynodes.cpp
HEADERS += guimirrormesh.h
SOURCES += guimirrormesh.cpp
FORMS += guimirrormesh.ui
