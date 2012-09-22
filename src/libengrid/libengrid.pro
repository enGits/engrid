TEMPLATE = lib
LANGUAGE = C++
TARGET   = engrid

# Enable this if the VTK from the ParaView sources and 
# installation want to be used
# Note: Currently only for Windows Compiles with MSVC
Use_VTK_Win_ParaView = yes


CONFIG += qt \
          debug_and_release \
          thread

QT     += xml \
          network \
          opengl

win32-msvc* {
    QMAKE_CXXFLAGS += -W3
    DEFINES += LIBENGRID_EXPORTS
    DEFINES += DLL_EXPORT
} win32-g++* {
    CONFIG += console
    DEFINES += LIBENGRID_EXPORTS
    DEFINES += DLL_EXPORT
    QMAKE_CXXFLAGS += -Wall
    QMAKE_CXXFLAGS += -Wno-deprecated
    QMAKE_CXXFLAGS += -Wl,--no-undefined
    QMAKE_CXXFLAGS += -Wl,--enable-runtime-pseudo-reloc
} else {
    QMAKE_CXXFLAGS += -Wall
    QMAKE_CXXFLAGS += -Wno-deprecated
    QMAKE_CXXFLAGS += -fno-omit-frame-pointer
    QMAKE_CXXFLAGS += -g
}


INCLUDEPATH += ..
INCLUDEPATH += ./libengrid-build
!debian {
    INCLUDEPATH += ../netgen_svn/netgen-mesher/netgen/nglib
    INCLUDEPATH += ../netgen_svn/netgen-mesher/netgen/libsrc/general
}

#INCLUDEPATH for VTK depends on the compiler
win32-msvc* {
    DEFINES += _USE_MATH_DEFINES

    !isEmpty(Use_VTK_Win_ParaView) {
        include(../misc/engrid-vtk-win_paraview.pri)
    } else {
        INCLUDEPATH += $(VTKINCDIR)
    }
} win32-g++* {
    INCLUDEPATH += $(VTKINCDIR)
} else {
    INCLUDEPATH += $(VTKINCDIR)
}

RESOURCES += engrid.qrc

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
    ../math/linsolve.h \
    ../math/mathvector.h \
    ../math/mathvector_methods.h \
    ../math/mathvector_operators.h \
    ../math/mathvector_structs.h \
    ../math/smallsquarematrix.h \
    pointfinder.h \
    createboundarylayer.h \
    guisurfacemesher.h \
    guicreatehexcore.h \
    createhexcore.h \
    orthogonalityoptimiser.h \
    optimisenormalvector.h \
    brlcadreader.h \
    eghashset.h \
    polymolecule.h \
    su2writer.h \
    booleangeometryoperation.h \
    guibooleangeometryoperation.h \
    dolfynwriter.h

SOURCES = boundarycondition.cpp \
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
    pointfinder.cpp \
    createboundarylayer.cpp \
    guisurfacemesher.cpp \
    guicreatehexcore.cpp \
    createhexcore.cpp \
    orthogonalityoptimiser.cpp \
    optimisenormalvector.cpp \
    brlcadreader.cpp \
    polymolecule.cpp \
    su2writer.cpp \
    booleangeometryoperation.cpp \
    guibooleangeometryoperation.cpp \
    dolfynwriter.cpp

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
    guicreatevolumemesh.ui \
    guisurfacemesher.ui \
    guicreatehexcore.ui \
    guibooleangeometryoperation.ui
    
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
HEADERS += dialoglineedit.h
SOURCES += dialoglineedit.cpp
HEADERS += utilities.h
SOURCES += utilities.cpp
HEADERS += edgelengthsourcemanager.h \
           edgelengthsource.h
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
FORMS += guiedgelengthsourcebox.ui
HEADERS += guiedgelengthsourcebox.h
SOURCES += guiedgelengthsourcebox.cpp
FORMS += guimergevolumes.ui
HEADERS += guimergevolumes.h
SOURCES += guimergevolumes.cpp
HEADERS += deletestraynodes.h
SOURCES += deletestraynodes.cpp
HEADERS += guimirrormesh.h
SOURCES += guimirrormesh.cpp
FORMS += guimirrormesh.ui
