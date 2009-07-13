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

//connect(ui.action,SIGNAL(activated()),this,SLOT(call()));

connect(ui.actionImportSTL,              SIGNAL(activated()),       this, SLOT(importSTL()));
connect(ui.actionImportGmsh1Ascii,       SIGNAL(activated()),       this, SLOT(importGmsh1Ascii()));
connect(ui.actionImportGmsh2Ascii,       SIGNAL(activated()),       this, SLOT(importGmsh2Ascii()));
connect(ui.actionExportGmsh1Ascii,       SIGNAL(activated()),       this, SLOT(exportGmsh1Ascii()));
connect(ui.actionExportGmsh2Ascii,       SIGNAL(activated()),       this, SLOT(exportGmsh2Ascii()));
connect(ui.actionExportNeutral,          SIGNAL(activated()),       this, SLOT(exportNeutral()));
connect(ui.actionExportAsciiStl,         SIGNAL(activated()),       this, SLOT(exportAsciiStl()));
connect(ui.actionExportBinaryStl,        SIGNAL(activated()),       this, SLOT(exportBinaryStl()));
connect(ui.actionExit,                   SIGNAL(activated()),       this, SLOT(exit()));
connect(ui.actionZoomAll,                SIGNAL(activated()),       this, SLOT(zoomAll()));
connect(ui.actionZoomOnPickedObject,     SIGNAL(activated()),       this, SLOT(zoomOnPickedObject()));
connect(ui.actionPrintGrid,              SIGNAL(activated()),       this, SLOT(printGrid()));
connect(ui.actionShowInfo,               SIGNAL(activated()),       this, SLOT(info()));
connect(ui.actionDeselectAll,            SIGNAL(activated()),       this, SLOT(deselectAll()));
connect(ui.actionOpen,                   SIGNAL(activated()),       this, SLOT(open()));
connect(ui.actionSave,                   SIGNAL(activated()),       this, SLOT(save()));
connect(ui.actionSaveAs,                 SIGNAL(activated()),       this, SLOT(saveAs()));
connect(ui.actionBoundaryCodes,          SIGNAL(activated()),       this, SLOT(selectBoundaryCodes()));
connect(ui.actionNormalExtrusion,        SIGNAL(activated()),       this, SLOT(normalExtrusion()));
connect(ui.actionViewAxes,               SIGNAL(changed()),         this, SLOT(setAxesVisibility()));
connect(ui.actionViewOrthogonal,         SIGNAL(changed()),         this, SLOT(setViewingMode()));
connect(ui.actionViewNodeIDs,            SIGNAL(changed()),         this, SLOT(viewNodeIDs()));
connect(ui.actionViewCellIDs,            SIGNAL(changed()),         this, SLOT(viewCellIDs()));
connect(ui.actionChangeOrientation,      SIGNAL(activated()),       this, SLOT(changeSurfaceOrientation()));
connect(ui.actionCheckOrientation,       SIGNAL(activated()),       this, SLOT(checkSurfaceOrientation()));
connect(ui.actionImproveAspectRatio,     SIGNAL(activated()),       this, SLOT(improveAspectRatio()));
connect(ui.actionRedraw,                 SIGNAL(activated()),       this, SLOT(updateActors()));
connect(ui.actionForcedRedraw,           SIGNAL(activated()),       this, SLOT(forceUpdateActors()));
connect(ui.actionScaleToData,            SIGNAL(activated()),       this, SLOT(scaleToData()));
connect(ui.actionClearOutputWindow,      SIGNAL(activated()),       this, SLOT(clearOutput()));
connect(ui.actionEditBoundaryConditions, SIGNAL(activated()),       this, SLOT(editBoundaryConditions()));
connect(ui.actionConfigure,              SIGNAL(activated()),       this, SLOT(configure()));
connect(ui.actionAbout,                  SIGNAL(activated()),       this, SLOT(about()));
connect(ui.actionStoreGeometry,          SIGNAL(activated()),       this, SLOT(callUpdateSurfProj()));

connect(ui.checkBox_UseVTKInteractor,    SIGNAL(stateChanged(int)), this, SLOT(setUseVTKInteractor(int)));

connect(ui.actionViewXP, SIGNAL(activated()), this, SLOT(viewXP()));
connect(ui.actionViewXM, SIGNAL(activated()), this, SLOT(viewXM()));
connect(ui.actionViewYP, SIGNAL(activated()), this, SLOT(viewYP()));
connect(ui.actionViewYM, SIGNAL(activated()), this, SLOT(viewYM()));
connect(ui.actionViewZP, SIGNAL(activated()), this, SLOT(viewZP()));
connect(ui.actionViewZM, SIGNAL(activated()), this, SLOT(viewZM()));

connect(ui.lineEditClipX, SIGNAL(textChanged(QString)), this, SLOT(setClipX(QString)));
connect(ui.lineEditClipY, SIGNAL(textChanged(QString)), this, SLOT(setClipY(QString)));
connect(ui.lineEditClipZ, SIGNAL(textChanged(QString)), this, SLOT(setClipZ(QString)));
connect(ui.lineEditClipNX, SIGNAL(textChanged(QString)), this, SLOT(setClipNX(QString)));
connect(ui.lineEditClipNY, SIGNAL(textChanged(QString)), this, SLOT(setClipNY(QString)));
connect(ui.lineEditClipNZ, SIGNAL(textChanged(QString)), this, SLOT(setClipNZ(QString)));

connect(ui.pushButtonMarkPosition, SIGNAL(clicked()), this, SLOT(markOutputLine()));

connect(ui.actionCreateSurfaceMesh,SIGNAL(activated()),this,SLOT(callCreateSurfaceMesh()));
connect(ui.actionCreateBoundaryLayer,SIGNAL(activated()),this,SLOT(callCreateBoundaryLayer()));
connect(ui.actionDivideBoundaryLayer,SIGNAL(activated()),this,SLOT(callDivideBoundaryLayer()));
connect(ui.actionDeleteVolumeGrid,SIGNAL(activated()),this,SLOT(callDeleteTetras()));
connect(ui.actionFixSTL,SIGNAL(activated()),this,SLOT(callFixSTL()));
connect(ui.actionCreateVolumeMesh,SIGNAL(activated()),this,SLOT(callCreateVolumeMesh()));
connect(ui.actionSmoothVolumeGrid,SIGNAL(activated()),this,SLOT(callSmoothVolumeGrid()));
connect(ui.actionVtkReader,SIGNAL(activated()),this,SLOT(callVtkReader()));
connect(ui.actionPolyDataReader,SIGNAL(activated()),this,SLOT(callPolyDataReader()));
connect(ui.actionSetBoundaryCode,SIGNAL(activated()),this,SLOT(callSetBoundaryCode()));
connect(ui.actionFoamWriter,SIGNAL(activated()),this,SLOT(callFoamWriter()));
connect(ui.actionSimpleFoamWriter,SIGNAL(activated()),this,SLOT(callSimpleFoamWriter()));
connect(ui.actionFoamCaseWriter, SIGNAL(activated()), this, SLOT(callFoamCaseWriter()));
connect(ui.actionDeleteBadAspectTris,SIGNAL(activated()),this,SLOT(callDeleteBadAspectTris()));
connect(ui.actionDeletePickedCell,SIGNAL(activated()),this,SLOT(callDeletePickedCell()));
connect(ui.actionDeletePickedPoint,SIGNAL(activated()),this,SLOT(callDeletePickedPoint()));
connect(ui.actionBoxSelect,SIGNAL(activated()),this,SLOT(callBoxSelect()));
connect(ui.actionCheck_surface_integrity,SIGNAL(activated()),this,SLOT(callCheckSurfaceIntegrity()));
connect(ui.actionPick_cell_point,SIGNAL(activated()),this,SLOT(callPick_cell_point()));
connect(ui.actionTransform, SIGNAL(activated()), this, SLOT(callTransform()));
connect(ui.actionExportCGNS, SIGNAL(activated()), this, SLOT(callCgnsWriter()));
connect(ui.actionUndo, SIGNAL(activated()), this, SLOT(undo()));
connect(ui.actionRedo, SIGNAL(activated()), this, SLOT(redo()));
connect(ui.actionImportOpenFoamCase, SIGNAL(activated()), this, SLOT(callImportOpenFoamCase()));

// OpenFOAMtools
connect(ui.actionRunSolver,             SIGNAL(activated()), &m_OpenFoamTools, SLOT(runSolver()));
// connect(ui.actionRunFoamToVTK,          SIGNAL(activated()), &m_OpenFoamTools, SLOT(runFoamToVTK()));
connect(ui.actionPreparePostProcessing, SIGNAL(activated()), &m_OpenFoamTools, SLOT(runPostProcessingTools()));
connect(ui.actionStopProcesses,         SIGNAL(activated()), &m_OpenFoamTools, SLOT(stopSolverProcess()));
connect(ui.actionImportFluentCase,      SIGNAL(activated()), &m_OpenFoamTools, SLOT(runImportFluentCase()));
// -------------------------------------------
