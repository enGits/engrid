//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008,2009 Oliver Gloth                                     +
// +                                                                      +
// + enGrid is free software: you can redistribute it and/or modify       +
// + it under the terms of the GNU General Public License as published by +
// + the Free Software Foundation, either version 3 of the License, or    +
// + (at your option) any later version.                                  +
// +                                                                      +
// + enGrid is distributed in the hope that it will be useful,            +
// + but WITHOUT ANY WARRANTY; without even the implied warranty of       +
// + MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        +
// + GNU General Public License for more details.                         +
// +                                                                      +
// + You should have received a copy of the GNU General Public License    +
// + along with enGrid. If not, see <http://www.gnu.org/licenses/>.       +
// +                                                                      +
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//

//connect(ui.action,SIGNAL(triggered()),this,SLOT(call()));

connect(ui.actionImportSTL,              SIGNAL(triggered()),       this, SLOT(importSTL()));
connect(ui.actionImportGmsh1Ascii,       SIGNAL(triggered()),       this, SLOT(importGmsh1Ascii()));
connect(ui.actionImportGmsh2Ascii,       SIGNAL(triggered()),       this, SLOT(importGmsh2Ascii()));
connect(ui.actionExportGmsh1Ascii,       SIGNAL(triggered()),       this, SLOT(exportGmsh1Ascii()));
connect(ui.actionExportGmsh2Ascii,       SIGNAL(triggered()),       this, SLOT(exportGmsh2Ascii()));
connect(ui.actionExportNeutral,          SIGNAL(triggered()),       this, SLOT(exportNeutral()));
connect(ui.actionExportAsciiStl,         SIGNAL(triggered()),       this, SLOT(exportAsciiStl()));
connect(ui.actionExportBinaryStl,        SIGNAL(triggered()),       this, SLOT(exportBinaryStl()));
connect(ui.actionExportAsciiPly,         SIGNAL(triggered()),       this, SLOT(exportAsciiPly()));
connect(ui.actionExportBinaryPly,        SIGNAL(triggered()),       this, SLOT(exportBinaryPly()));
connect(ui.actionExit,                   SIGNAL(triggered()),       this, SLOT(exit()));
connect(ui.actionZoomAll,                SIGNAL(triggered()),       this, SLOT(zoomAll()));
connect(ui.actionZoomOnPickedObject,     SIGNAL(triggered()),       this, SLOT(zoomOnPickedObject()));
connect(ui.actionPrintGrid,              SIGNAL(triggered()),       this, SLOT(printGrid()));
connect(ui.actionShowInfo,               SIGNAL(triggered()),       this, SLOT(info()));
connect(ui.actionDeselectAll,            SIGNAL(triggered()),       this, SLOT(deselectAll()));
connect(ui.actionOpen,                   SIGNAL(triggered()),       this, SLOT(open()));
connect(ui.actionSave,                   SIGNAL(triggered()),       this, SLOT(save()));
connect(ui.actionSaveAs,                 SIGNAL(triggered()),       this, SLOT(saveAs()));
connect(ui.actionBoundaryCodes,          SIGNAL(triggered()),       this, SLOT(selectBoundaryCodes()));
connect(ui.actionNormalExtrusion,        SIGNAL(triggered()),       this, SLOT(normalExtrusion()));
connect(ui.actionViewAxes,               SIGNAL(changed()),         this, SLOT(setAxesVisibility()));
connect(ui.actionViewOrthogonal,         SIGNAL(changed()),         this, SLOT(setViewingMode()));
connect(ui.actionViewNodeIDs,            SIGNAL(changed()),         this, SLOT(viewNodeIDs()));
connect(ui.actionViewCellIDs,            SIGNAL(changed()),         this, SLOT(viewCellIDs()));
connect(ui.actionChangeOrientation,      SIGNAL(triggered()),       this, SLOT(changeSurfaceOrientation()));
connect(ui.actionCheckOrientation,       SIGNAL(triggered()),       this, SLOT(checkSurfaceOrientation()));
connect(ui.actionImproveAspectRatio,     SIGNAL(triggered()),       this, SLOT(improveAspectRatio()));
connect(ui.actionRedraw,                 SIGNAL(triggered()),       this, SLOT(updateActors()));
connect(ui.actionForcedRedraw,           SIGNAL(triggered()),       this, SLOT(forceUpdateActors()));
connect(ui.actionScaleToData,            SIGNAL(triggered()),       this, SLOT(scaleToData()));
connect(ui.actionClearOutputWindow,      SIGNAL(triggered()),       this, SLOT(clearOutput()));
connect(ui.actionEditBoundaryConditions, SIGNAL(triggered()),       this, SLOT(editBoundaryConditions()));
connect(ui.actionConfigure,              SIGNAL(triggered()),       this, SLOT(configure()));
connect(ui.actionAbout,                  SIGNAL(triggered()),       this, SLOT(about()));
connect(ui.actionStoreGeometry,          SIGNAL(triggered()),       this, SLOT(callUpdateSurfProj()));

connect(ui.checkBox_UseVTKInteractor,    SIGNAL(stateChanged(int)), this, SLOT(setUseVTKInteractor(int)));

connect(ui.actionViewXP, SIGNAL(triggered()), this, SLOT(viewXP()));
connect(ui.actionViewXM, SIGNAL(triggered()), this, SLOT(viewXM()));
connect(ui.actionViewYP, SIGNAL(triggered()), this, SLOT(viewYP()));
connect(ui.actionViewYM, SIGNAL(triggered()), this, SLOT(viewYM()));
connect(ui.actionViewZP, SIGNAL(triggered()), this, SLOT(viewZP()));
connect(ui.actionViewZM, SIGNAL(triggered()), this, SLOT(viewZM()));

connect(ui.lineEditClipX, SIGNAL(textChanged(QString)), this, SLOT(setClipX(QString)));
connect(ui.lineEditClipY, SIGNAL(textChanged(QString)), this, SLOT(setClipY(QString)));
connect(ui.lineEditClipZ, SIGNAL(textChanged(QString)), this, SLOT(setClipZ(QString)));
connect(ui.lineEditClipNX, SIGNAL(textChanged(QString)), this, SLOT(setClipNX(QString)));
connect(ui.lineEditClipNY, SIGNAL(textChanged(QString)), this, SLOT(setClipNY(QString)));
connect(ui.lineEditClipNZ, SIGNAL(textChanged(QString)), this, SLOT(setClipNZ(QString)));

connect(ui.pushButtonMarkPosition, SIGNAL(clicked()), this, SLOT(markOutputLine()));

connect(ui.actionCreateSurfaceMesh,SIGNAL(triggered()),this,SLOT(callCreateSurfaceMesh()));
connect(ui.actionCreateBoundaryLayer,SIGNAL(triggered()),this,SLOT(callCreateBoundaryLayer()));
connect(ui.actionDivideBoundaryLayer,SIGNAL(triggered()),this,SLOT(callDivideBoundaryLayer()));
connect(ui.actionDeleteVolumeGrid,SIGNAL(triggered()),this,SLOT(callDeleteTetras()));
connect(ui.actionFixSTL,SIGNAL(triggered()),this,SLOT(callFixSTL()));
connect(ui.actionCreateVolumeMesh,SIGNAL(triggered()),this,SLOT(callCreateVolumeMesh()));
connect(ui.actionSmoothVolumeGrid,SIGNAL(triggered()),this,SLOT(callSmoothVolumeGrid()));
connect(ui.actionVtkReader,SIGNAL(triggered()),this,SLOT(callVtkReader()));
connect(ui.actionPolyDataReader,SIGNAL(triggered()),this,SLOT(callPolyDataReader()));
connect(ui.actionSetBoundaryCode,SIGNAL(triggered()),this,SLOT(callSetBoundaryCode()));
connect(ui.actionFoamWriter,SIGNAL(triggered()),this,SLOT(callFoamWriter()));
connect(ui.actionSimpleFoamWriter,SIGNAL(triggered()),this,SLOT(callSimpleFoamWriter()));
connect(ui.actionFoamCaseWriter, SIGNAL(triggered()), this, SLOT(callFoamCaseWriter()));
connect(ui.actionDeleteBadAspectTris,SIGNAL(triggered()),this,SLOT(callDeleteBadAspectTris()));
connect(ui.actionDeletePickedCell,SIGNAL(triggered()),this,SLOT(callDeletePickedCell()));
connect(ui.actionDeletePickedPoint,SIGNAL(triggered()),this,SLOT(callDeletePickedPoint()));
connect(ui.actionBoxSelect,SIGNAL(triggered()),this,SLOT(callBoxSelect()));
connect(ui.actionCheck_surface_integrity,SIGNAL(triggered()),this,SLOT(callCheckSurfaceIntegrity()));
connect(ui.actionPick_cell_point,SIGNAL(triggered()),this,SLOT(callPick_cell_point()));
connect(ui.actionTransform, SIGNAL(triggered()), this, SLOT(callTransform()));
connect(ui.actionExportCGNS, SIGNAL(triggered()), this, SLOT(callCgnsWriter()));
connect(ui.actionUndo, SIGNAL(triggered()), this, SLOT(undo()));
connect(ui.actionRedo, SIGNAL(triggered()), this, SLOT(redo()));
connect(ui.actionImportOpenFoamCase, SIGNAL(triggered()), this, SLOT(callImportOpenFoamCase()));
connect(ui.actionReducedPolyDataReader, SIGNAL(triggered()), this, SLOT(callReducedPolyDataReader()));
connect(ui.actionSurfaceMesher, SIGNAL(triggered()), this, SLOT(callSurfaceMesher()));
connect(ui.actionReduceSurfaceTriangulation, SIGNAL(triggered()), this, SLOT(callReduceSurfaceTriangulation()));
connect(ui.actionEliminateSmallBranches, SIGNAL(triggered()), this, SLOT(callEliminateSmallBranches()));
connect(ui.actionSmoothAndSwapSurface, SIGNAL(triggered()), this, SLOT(callSmoothAndSwapSurface()));
connect(ui.actionImportSeligAirfoil, SIGNAL(triggered()), this, SLOT(callSeligAirfoilReader()));
connect(ui.actionImportBlenderFile, SIGNAL(triggered()), this, SLOT(callBlenderReader()));

connect(ui.actionFixCADgeometry, SIGNAL(triggered()), this, SLOT(callFixCADGeometry()));

// OpenFOAMtools
connect(ui.actionRunSolver,             SIGNAL(triggered()), &m_OpenFoamTools, SLOT(runSolver()));
// connect(ui.actionRunFoamToVTK,          SIGNAL(triggered()), &m_OpenFoamTools, SLOT(runFoamToVTK()));
connect(ui.actionPreparePostProcessing, SIGNAL(triggered()), &m_OpenFoamTools, SLOT(runPostProcessingTools()));
connect(ui.actionStopProcesses,         SIGNAL(triggered()), &m_OpenFoamTools, SLOT(stopSolverProcess()));
connect(ui.actionImportFluentCase,      SIGNAL(triggered()), &m_OpenFoamTools, SLOT(runImportFluentCase()));

connect(ui.menuOpen_recent,SIGNAL(triggered(QAction*)),this,SLOT(openRecent(QAction*)));
// -------------------------------------------
