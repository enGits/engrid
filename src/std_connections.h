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

connect(ui.actionSmoothSurface,SIGNAL(activated()),this,SLOT(callSmoothSurface()));
connect(ui.actionCreateBoundaryLayer,SIGNAL(activated()),this,SLOT(callCreateBoundaryLayer()));
connect(ui.actionDivideBoundaryLayer,SIGNAL(activated()),this,SLOT(callDivideBoundaryLayer()));
connect(ui.actionDeleteVolumeGrid,SIGNAL(activated()),this,SLOT(callDeleteTetras()));
connect(ui.actionFixSTL,SIGNAL(activated()),this,SLOT(callFixSTL()));
connect(ui.actionCreateVolumeMesh,SIGNAL(activated()),this,SLOT(callCreateVolumeMesh()));
connect(ui.actionSmoothVolumeGrid,SIGNAL(activated()),this,SLOT(callSmoothVolumeGrid()));
connect(ui.actionFoamReader,SIGNAL(activated()),this,SLOT(callFoamReader()));
connect(ui.actionVtkReader,SIGNAL(activated()),this,SLOT(callVtkReader()));
connect(ui.actionPolyDataReader,SIGNAL(activated()),this,SLOT(callPolyDataReader()));
connect(ui.actionSetBoundaryCode,SIGNAL(activated()),this,SLOT(callSetBoundaryCode()));
connect(ui.actionFoamWriter,SIGNAL(activated()),this,SLOT(callFoamWriter()));
connect(ui.actionSimpleFoamWriter,SIGNAL(activated()),this,SLOT(callSimpleFoamWriter()));
connect(ui.actionDeleteBadAspectTris,SIGNAL(activated()),this,SLOT(callDeleteBadAspectTris()));
connect(ui.actionDeletePickedCell,SIGNAL(activated()),this,SLOT(callDeletePickedCell()));
connect(ui.actionDeletePickedPoint,SIGNAL(activated()),this,SLOT(callDeletePickedPoint()));
connect(ui.actionBoxSelect,SIGNAL(activated()),this,SLOT(callBoxSelect()));
connect(ui.actionCheck_surface_integrity,SIGNAL(activated()),this,SLOT(callCheckSurfaceIntegrity()));
connect(ui.actionPick_cell_point,SIGNAL(activated()),this,SLOT(callPick_cell_point()));
connect(ui.actionTransform, SIGNAL(activated()), this, SLOT(callTransform()));
connect(ui.actionExportCGNS, SIGNAL(activated()), this, SLOT(callCgnsWriter()));
connect(ui.actionUndo, SIGNAL(activated()), this, SLOT(Undo()));
connect(ui.actionRedo, SIGNAL(activated()), this, SLOT(Redo()));

// -------------------------------------------
