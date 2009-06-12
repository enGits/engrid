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
#include "guimainwindow.h"
#include "guiselectboundarycodes.h"
#include "guiimproveaspectratio.h"
#include "guinormalextrusion.h"
#include "guisetboundarycode.h"
#include "guipick.h"

#include "vtkEgPolyDataToUnstructuredGridFilter.h"
#include "stlreader.h"
#include "gmshreader.h"
#include "gmshwriter.h"
#include "neutralwriter.h"
#include "stlwriter.h"
#include "correctsurfaceorientation.h"
#include "guieditboundaryconditions.h"

#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkProperty.h>
#include <vtkCallbackCommand.h>
#include <vtkCamera.h>
#include <vtkCharArray.h>
#include <vtkTextActor.h>
#include <vtkVectorText.h>
#include <vtkFollower.h>
#include <vtkLookupTable.h>
#include <vtkScalarBarActor.h>

#include <QFileDialog>
#include <QFileSystemWatcher>
#include <stdlib.h>
#include <stdio.h>

#include "geometrytools.h"
using namespace GeometryTools;

#include "guisettingsviewer.h"
#include "guitransform.h"
#include "egvtkinteractorstyle.h"
#include "showinfo.h"

QString GuiMainWindow::cwd = ".";
QSettings GuiMainWindow::qset("enGits","enGrid");
GuiMainWindow* GuiMainWindow::THIS = NULL;
QMutex GuiMainWindow::mutex;
vtkIdType GuiMainWindow::PickedPoint;
vtkIdType GuiMainWindow::PickedCell;
bool GuiMainWindow::m_UseVTKInteractor;

GuiMainWindow::GuiMainWindow() : QMainWindow(NULL)
{
  ui.setupUi(this);
  THIS = this;
  
  setGeometry(qset.value("GuiMainWindow", QRect(200,200,400,400)).toRect());
  restoreState(qset.value("dockWidget_states").toByteArray());
  
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
  
#include "std_connections.h"
  
  if (qset.contains("working_directory")) {
    cwd = qset.value("working_directory").toString();
  }

  setupVtk();

  current_filename = "untitled.vtu";
  ResetOperationCounter();//clears undo/redo list and disables undo/redo
  setWindowTitle(current_filename + " - enGrid - " + QString("%1").arg(current_operation) );
  
  status_bar = new QStatusBar(this);
  setStatusBar(status_bar);
  status_label = new QLabel(this);
  status_bar->addWidget(status_label);

  QString txt = "0 volume cells (0 tetras, 0 hexas, 0 pyramids, 0 prisms), ";
  txt += "0 surface cells (0 triangles, 0 quads), 0 nodes";
  status_label->setText(txt);
  ui.label_node_cell_info->setText(txt);

  QString user = QString(getenv("USER"));
  QString basename="enGrid_output.txt";
  
  // define temporary path
  QDir dir("/");
  if (qset.contains("tmp_directory")) {
    m_LogDir = qset.value("tmp_directory").toString();
  } else {
    m_LogDir = dir.tempPath();
  }
  QDateTime now = QDateTime::currentDateTime();
  m_LogDir = m_LogDir + "/" + "enGrid_" + QDateTime::currentDateTime().toString("yyyyMMddhhmmsszzz") + "/";
  dir.mkpath(m_LogDir);
  
  log_file_name = m_LogDir + basename;
  cout << "log_file_name=" << log_file_name.toLatin1().data() << endl;

  system_stdout = stdout;
  freopen (log_file_name.toAscii().data(), "w", stdout);
  
  busy = false;
  
  setPickMode(true,true);
  PickedPoint=-1;
  PickedCell=-1;
  
  updateStatusBar();
  
  connect(&garbage_timer, SIGNAL(timeout()), this, SLOT(periodicUpdate()));
  garbage_timer.start(1000);
  
  connect(&log_timer, SIGNAL(timeout()), this, SLOT(updateOutput()));
  log_timer.start(1000);
  
  N_chars = 0;
  
  bool exp_features=false;
  getSet("","enable experimental features",false,exp_features);
  bool undo_redo;
  getSet("","enable undo/redo",false,undo_redo);
  bool undo_redo_mode;
  getSet("","use RAM for undo/redo oprations",false,undo_redo_mode);
  
  ui.actionFoamWriter->setEnabled(exp_features);
  
  m_ReferenceSize=0.2;
  
  ui.doubleSpinBox_HueMin->setValue(0.667);
  ui.doubleSpinBox_HueMax->setValue(0);
  
  egvtkInteractorStyle *style = egvtkInteractorStyle::New();
  getInteractor()->SetInteractorStyle(style);
  style->Delete();
}

//end of GuiMainWindow::GuiMainWindow() : QMainWindow(NULL)

GuiMainWindow::~GuiMainWindow()
{
  qset.setValue("GuiMainWindow", this->geometry());
  qset.setValue("dockWidget_states", this->saveState());
  
#ifndef QT_DEBUG
  QDirIterator it(m_LogDir);
  while (it.hasNext()) {
    QString str = it.next();
    QFileInfo fileinfo(str);
    if(fileinfo.isFile()) {
      QFile file(str);
      if(!file.remove()) qDebug() << "Failed to remove " << file.fileName();
    }
  }
  QDir dir(m_LogDir);
  dir.rmdir(m_LogDir);
#endif
  
}

void GuiMainWindow::setupVtk()
{
  grid = vtkUnstructuredGrid::New();
  renderer = vtkRenderer::New();
  getRenderWindow()->AddRenderer(renderer);

  // coordinate axes
  axes = vtkCubeAxesActor2D::New();
  //
  axes->SetCamera(getRenderer()->GetActiveCamera());
  getRenderer()->AddActor(axes);
  axes->SetVisibility(0);

  // surface pipelines
  backface_property   = vtkProperty::New();
  surface_filter      = vtkGeometryFilter::New();
  m_SurfaceMapper     = vtkPolyDataMapper::New();
  m_SurfaceWireMapper = vtkPolyDataMapper::New();
  m_BCodesFilter      = vtkEgBoundaryCodesFilter::New();
  lut                 = vtkLookupTable::New();
  m_SurfaceActor      = vtkActor::New();
  m_SurfaceWireActor  = vtkActor::New();
  m_LegendActor       = vtkScalarBarActor::New();
  //
  m_BCodesFilter->SetBoundaryCodes(m_DisplayBoundaryCodes);
  m_BCodesFilter->SetInput(grid);
  surface_filter->SetInput(m_BCodesFilter->GetOutput());
  m_SurfaceMapper->SetInput(surface_filter->GetOutput());
  m_SurfaceWireMapper->SetInput(surface_filter->GetOutput());
  m_SurfaceMapper->SetLookupTable(lut);
  m_SurfaceActor->GetProperty()->SetRepresentationToSurface();
  m_SurfaceActor->GetProperty()->SetColor(0.5,1,0.5);
  m_SurfaceActor->SetBackfaceProperty(backface_property);
  m_SurfaceActor->GetBackfaceProperty()->SetColor(1,1,0.5);
  m_SurfaceActor->SetMapper(m_SurfaceMapper);
  getRenderer()->AddActor(m_SurfaceActor);
  m_SurfaceActor->SetVisibility(1);
  m_LegendActor->SetLookupTable(lut);
  getRenderer()->AddActor(m_LegendActor);
  m_LegendActor->SetVisibility(0);
  m_SurfaceWireActor->GetProperty()->SetRepresentationToWireframe();
  m_SurfaceWireActor->GetProperty()->SetColor(0,0,1);
  m_SurfaceWireActor->SetMapper(m_SurfaceWireMapper);
  getRenderer()->AddActor(m_SurfaceWireActor);
  m_SurfaceWireActor->SetVisibility(1);

  // tetra pipline
  m_ExtrTetras   = vtkEgExtractVolumeCells::New();
  m_TetraActor   = vtkActor::New();
  tetra_geometry = vtkGeometryFilter::New();
  m_TetraMapper  = vtkPolyDataMapper::New();
  //
  m_ExtrTetras->SetInput(grid);
  m_ExtrTetras->SetAllOff();
  m_ExtrTetras->SetTetrasOn();;
  tetra_geometry->SetInput(m_ExtrTetras->GetOutput());
  m_TetraMapper->SetInput(tetra_geometry->GetOutput());
  m_TetraActor->SetMapper(m_TetraMapper);
  m_TetraActor->GetProperty()->SetColor(1,0,0);
  getRenderer()->AddActor(m_TetraActor);
  m_TetraActor->SetVisibility(0);

  // pyramid pipeline
  m_PyramidActor   = vtkActor::New();
  m_ExtrPyramids   = vtkEgExtractVolumeCells::New();
  pyramid_geometry = vtkGeometryFilter::New();
  m_PyramidMapper  = vtkPolyDataMapper::New();
  //
  m_ExtrPyramids->SetInput(grid);
  m_ExtrPyramids->SetAllOff();
  m_ExtrPyramids->SetPyramidsOn();
  pyramid_geometry->SetInput(m_ExtrPyramids->GetOutput());
  m_PyramidMapper->SetInput(pyramid_geometry->GetOutput());
  m_PyramidActor->SetMapper(m_PyramidMapper);
  m_PyramidActor->GetProperty()->SetColor(1,1,0);
  getRenderer()->AddActor(m_PyramidActor);
  m_PyramidActor->SetVisibility(0);

  // wedge pipeline
  m_WedgeActor   = vtkActor::New();
  m_ExtrWedges   = vtkEgExtractVolumeCells::New();
  wedge_geometry = vtkGeometryFilter::New();
  m_WedgeMapper  = vtkPolyDataMapper::New();
  //
  m_ExtrWedges->SetInput(grid);
  m_ExtrWedges->SetAllOff();
  m_ExtrWedges->SetWedgesOn();
  wedge_geometry->SetInput(m_ExtrWedges->GetOutput());
  m_WedgeMapper->SetInput(wedge_geometry->GetOutput());
  m_WedgeActor->SetMapper(m_WedgeMapper);
  m_WedgeActor->GetProperty()->SetColor(0,1,0);
  getRenderer()->AddActor(m_WedgeActor);
  m_WedgeActor->SetVisibility(0);

  // hexa pipeline
  m_HexaActor   = vtkActor::New();
  m_ExtrHexes   = vtkEgExtractVolumeCells::New();
  hexa_geometry = vtkGeometryFilter::New();
  m_HexaMapper  = vtkPolyDataMapper::New();
  //
  m_ExtrHexes->SetInput(grid);
  m_ExtrHexes->SetAllOff();
  m_ExtrHexes->SetHexesOn();
  hexa_geometry->SetInput(m_ExtrHexes->GetOutput());
  m_HexaMapper->SetInput(hexa_geometry->GetOutput());
  m_HexaActor->SetMapper(m_HexaMapper);
  m_HexaActor->GetProperty()->SetColor(0,0.7,1);
  getRenderer()->AddActor(m_HexaActor);
  m_HexaActor->SetVisibility(0);

  // volume wire pipeline
  m_VolumeWireActor  = vtkActor::New();
  m_ExtrVol          = vtkEgExtractVolumeCells::New();
  volume_geometry    = vtkGeometryFilter::New();
  m_VolumeWireMapper = vtkPolyDataMapper::New();
  //
  m_ExtrVol->SetInput(grid);
  m_ExtrVol->SetAllOn();
  volume_geometry->SetInput(m_ExtrVol->GetOutput());
  m_VolumeWireMapper->SetInput(volume_geometry->GetOutput());
  m_VolumeWireActor->SetMapper(m_VolumeWireMapper);
  m_VolumeWireActor->GetProperty()->SetRepresentationToWireframe();
  m_VolumeWireActor->GetProperty()->SetColor(0,0,1);
  getRenderer()->AddActor(m_VolumeWireActor);
  m_VolumeWireActor->SetVisibility(0);

  // picker stuff
  pick_sphere = vtkSphereSource::New();
  pick_mapper = vtkPolyDataMapper::New();
  pick_actor  = vtkActor::New();
  CellPicker  = vtkCellPicker::New();
  PointPicker = vtkPointPicker::New();
  //
  pick_sphere->SetRadius(0.25); //in case the user starts picking points instead of cells
  vtkCallbackCommand *cbc = vtkCallbackCommand::New();
  cbc->SetCallback(pickCellCallBack);
  CellPicker->AddObserver(vtkCommand::EndPickEvent, cbc);
  cbc->SetCallback(pickPointCallBack);
  PointPicker->AddObserver(vtkCommand::EndPickEvent, cbc);
  cbc->Delete();

}

void GuiMainWindow::updateOutput()
{
  QFile log_file(log_file_name);
  log_file.open(QIODevice::ReadOnly);
  QByteArray buffer = log_file.readAll();
  if (buffer.size() > N_chars) {
    QByteArray newchars = buffer.right(buffer.size() - N_chars);
    N_chars = buffer.size();
    QString txt(newchars);
    if (txt.right(1) == "\n") {
      txt = txt.left(txt.size()-1);
    }
    ui.textEditOutput->append(txt);
  }
}

void GuiMainWindow::exit()
{
  QCoreApplication::exit();
}

vtkRenderWindow* GuiMainWindow::getRenderWindow() 
{
  return ui.qvtkWidget->GetRenderWindow();
}

vtkRenderer* GuiMainWindow::getRenderer()
{
  return renderer;
}

QVTKInteractor* GuiMainWindow::getInteractor()
{
  return ui.qvtkWidget->GetInteractor();
}

QString GuiMainWindow::getCwd()
{
  return cwd;
}

void GuiMainWindow::setCwd(QString dir)
{
  cwd = dir;
  qset.setValue("working_directory",dir);
}

void GuiMainWindow::scaleToData()
{
  int current_field=ui.comboBox_Field->currentIndex();
  if(current_field>0)
  {
    double range[2];

    surface_filter->GetOutput()->GetPointData()->GetArray(current_field-1)->GetRange(range);
    //boundary_pd->GetPointData()->GetArray(current_field-1)->GetRange(range);
    cout<<"current_field="<<current_field<<endl;
    cout<<"range[0]="<<range[0]<<endl;
    cout<<"range[1]="<<range[1]<<endl;
    ui.doubleSpinBox_FieldMin->setRange(range[0],range[1]);
    ui.doubleSpinBox_FieldMax->setRange(range[0],range[1]);
    ui.doubleSpinBox_FieldMin->setValue(range[0]);
    ui.doubleSpinBox_FieldMax->setValue(range[1]);
  }
}

void GuiMainWindow::setClipX(const QString &txt)
{
  m_ExtrVol->Setx(txt.toDouble());
  m_ExtrTetras->Setx(txt.toDouble());
  m_ExtrPyramids->Setx(txt.toDouble());
  m_ExtrWedges->Setx(txt.toDouble());
  m_ExtrHexes->Setx(txt.toDouble());
}

void GuiMainWindow::setClipY(const QString &txt)
{
  m_ExtrVol->Sety(txt.toDouble());
  m_ExtrTetras->Sety(txt.toDouble());
  m_ExtrPyramids->Sety(txt.toDouble());
  m_ExtrWedges->Sety(txt.toDouble());
  m_ExtrHexes->Sety(txt.toDouble());
}

void GuiMainWindow::setClipZ(const QString &txt)
{
  m_ExtrVol->Setz(txt.toDouble());
  m_ExtrTetras->Setz(txt.toDouble());
  m_ExtrPyramids->Setz(txt.toDouble());
  m_ExtrWedges->Setz(txt.toDouble());
  m_ExtrHexes->Setz(txt.toDouble());
}

void GuiMainWindow::setClipNX(const QString &txt)
{
  m_ExtrVol->Setnx(txt.toDouble());
  m_ExtrTetras->Setnx(txt.toDouble());
  m_ExtrPyramids->Setnx(txt.toDouble());
  m_ExtrWedges->Setnx(txt.toDouble());
  m_ExtrHexes->Setnx(txt.toDouble());
}

void GuiMainWindow::setClipNY(const QString &txt)
{
  m_ExtrVol->Setny(txt.toDouble());
  m_ExtrTetras->Setny(txt.toDouble());
  m_ExtrPyramids->Setny(txt.toDouble());
  m_ExtrWedges->Setny(txt.toDouble());
  m_ExtrHexes->Setny(txt.toDouble());
}

void GuiMainWindow::setClipNZ(const QString &txt)
{
  m_ExtrVol->Setnz(txt.toDouble());
  m_ExtrTetras->Setnz(txt.toDouble());
  m_ExtrPyramids->Setnz(txt.toDouble());
  m_ExtrWedges->Setnz(txt.toDouble());
  m_ExtrHexes->Setnz(txt.toDouble());
}

void GuiMainWindow::updateSurfaceActors(bool forced)
{
  if (ui.checkBoxSurface->isChecked()) {
    if (forced) {
      surface_filter->Update();
    }

    // fill node field combobox
    int current_field=ui.comboBox_Field->currentIndex();
    ui.comboBox_Field->clear();
    ui.comboBox_Field->addItem("None");
    for (int i = 0; i < grid->GetPointData()->GetNumberOfArrays(); ++i) {
      ui.comboBox_Field->addItem(grid->GetPointData()->GetArrayName(i));
    }
    if(current_field == -1) {
      ui.comboBox_Field->setCurrentIndex(0);
    } else {
      ui.comboBox_Field->setCurrentIndex(current_field);
    }

    // fill cell field combobox
    int current_cell_field = ui.comboBox_CellTextField->currentIndex();
    ui.comboBox_CellTextField->clear();
    ui.comboBox_CellTextField->addItem("Cell ID");
    for (int i = 0; i < surface_filter->GetOutput()->GetCellData()->GetNumberOfArrays(); ++i) {
      ui.comboBox_CellTextField->addItem(grid->GetCellData()->GetArrayName(i));
    }
    if(current_cell_field == -1) {
      ui.comboBox_CellTextField->setCurrentIndex(0);
    } else {
      ui.comboBox_CellTextField->setCurrentIndex(current_cell_field);
    }
    current_field = ui.comboBox_Field->currentIndex();
    if(current_field > 0) {
      double range[2];
      surface_filter->GetOutput()->GetPointData()->GetArray(current_field-1)->GetRange(range);
      ui.doubleSpinBox_FieldMin->setRange(range[0],range[1]);
      ui.doubleSpinBox_FieldMax->setRange(range[0],range[1]);
    }

    if(ui.comboBox_Field->currentIndex() > 0) {
      m_SurfaceMapper->SetColorModeToMapScalars();
      lut->SetNumberOfColors(ui.spinBox_Color->value());
      lut->SetHueRange(ui.doubleSpinBox_HueMin->value(),ui.doubleSpinBox_HueMax->value());
      lut->Build();
      m_SurfaceMapper->SetScalarModeToUsePointFieldData();
      m_SurfaceMapper->ColorByArrayComponent(ui.comboBox_Field->currentText().toLatin1().data(),0);
      m_SurfaceMapper->SetScalarRange(ui.doubleSpinBox_FieldMin->value(),ui.doubleSpinBox_FieldMax->value());
      m_SurfaceMapper->ScalarVisibilityOn();
      if(ui.checkBox_Legend->checkState()) {
        m_LegendActor->SetVisibility(1);
      } else {
        m_LegendActor->SetVisibility(0);
      }
    } else {
      m_SurfaceMapper->SetColorModeToDefault();
      m_SurfaceMapper->ScalarVisibilityOff();
    }
    if (forced) {
      m_BCodesFilter->Update();
    }
    if(ui.checkBox_ShowPickSphere->checkState()) {
      if(m_UseVTKInteractor) {
        if(ui.radioButton_CellPicker->isChecked()) {
          getInteractor()->SetPicker(CellPicker);
          vtkIdType cellId = getPickedCell();
          pickCell(cellId);
        } else {
          getInteractor()->SetPicker(PointPicker);
          vtkIdType nodeId = getPickedPoint();
          pickPoint(nodeId);
        }
      } else {
        if(ui.radioButton_CellPicker->isChecked()) pickCell(PickedCell);
        else pickPoint(PickedPoint);
      }
    }
    m_SurfaceActor->SetVisibility(1);
    m_SurfaceWireActor->SetVisibility(1);
  } else {
    m_SurfaceActor->SetVisibility(0);
    m_SurfaceWireActor->SetVisibility(0);
  }
}

void GuiMainWindow::updateVolumeActors(bool forced)
{
  if (ui.checkBoxVolume->isChecked()) {
    if (ui.checkBoxTetra->isChecked()) {
      m_ExtrVol->SetTetrasOn();
      if (ui.checkBoxClip->isChecked()) {
        m_ExtrTetras->SetClippingOn();
      } else {
        m_ExtrTetras->SetClippingOff();
      }
      if (forced) {
        tetra_geometry->Update();
      }
      m_TetraActor->SetVisibility(1);
    } else {
      m_ExtrVol->SetTetrasOff();
      m_TetraActor->SetVisibility(0);
    }
    if (ui.checkBoxPyramid->isChecked()) {
      m_ExtrVol->SetPyramidsOn();
      if (ui.checkBoxClip->isChecked()) {
        m_ExtrPyramids->SetClippingOn();
      } else {
        m_ExtrPyramids->SetClippingOff();
      }
      if (forced) {
        pyramid_geometry->Update();
      }
      m_PyramidActor->SetVisibility(1);
    } else {
      m_ExtrVol->SetPyramidsOff();
      m_PyramidActor->SetVisibility(0);
    }
    if (ui.checkBoxWedge->isChecked()) {
      m_ExtrVol->SetWedgesOn();
      if (ui.checkBoxClip->isChecked()) {
        m_ExtrWedges->SetClippingOn();
      } else {
        m_ExtrWedges->SetClippingOff();
      }
      if (forced) {
        wedge_geometry->Update();
      }
      m_WedgeActor->SetVisibility(1);
    } else {
      m_ExtrVol->SetWedgesOff();
      m_WedgeActor->SetVisibility(0);
    }
    if (ui.checkBoxHexa->isChecked()) {
      m_ExtrVol->SetHexesOn();
      if (ui.checkBoxClip->isChecked()) {
        m_ExtrHexes->SetClippingOn();
      } else {
        m_ExtrHexes->SetClippingOff();
      }
      if (forced) {
        hexa_geometry->Update();
      }
      m_HexaActor->SetVisibility(1);
    } else {
      m_ExtrVol->SetHexesOff();
      m_HexaActor->SetVisibility(0);
    }

    // wireframe
    if (ui.checkBoxClip->isChecked()) {
      m_ExtrVol->SetClippingOn();
    } else {
      m_ExtrVol->SetClippingOff();
    }
    if (forced) {
      volume_geometry->Update();
    }
    m_VolumeWireActor->SetVisibility(1);
  } else {
    m_TetraActor->VisibilityOff();
    m_PyramidActor->VisibilityOff();
    m_WedgeActor->VisibilityOff();
    m_HexaActor->VisibilityOff();
    m_VolumeWireActor->VisibilityOff();
  }
}

void GuiMainWindow::updateActors(bool forced)
{
  if (!tryLock()) return;
  try {
    axes->SetInput(grid);
    updateSurfaceActors(forced);
    updateVolumeActors(forced);
    updateStatusBar();
    getRenderWindow()->Render();
  } catch (Error err) {
    err.display();
  }
  unlock();
}



void GuiMainWindow::forceUpdateActors()
{
  updateActors(true);
}

void GuiMainWindow::setPickMode(bool a_UseVTKInteractor,bool a_CellPickerMode)
{
  m_UseVTKInteractor=a_UseVTKInteractor;
  if (a_UseVTKInteractor) {
    ui.checkBox_UseVTKInteractor->setCheckState(Qt::Checked);
  } else {
    ui.checkBox_UseVTKInteractor->setCheckState(Qt::Unchecked);
  }
  if (a_CellPickerMode) {
    ui.radioButton_CellPicker->toggle();
  } else {
    ui.radioButton_PointPicker->toggle();
  }
}

void GuiMainWindow::setUseVTKInteractor(int a_UseVTKInteractor)
{
  m_UseVTKInteractor = a_UseVTKInteractor;
}

bool GuiMainWindow::pickPoint(vtkIdType nodeId)
{
  if (nodeId >= 0 && nodeId<grid->GetNumberOfPoints()) {
    vec3_t x(0,0,0);
    grid->GetPoints()->GetPoint(nodeId, x.data());
    pick_sphere->SetCenter(x.data());
    pick_mapper->SetInput(pick_sphere->GetOutput());
    
    if (pick_actor) {
      getRenderer()->RemoveActor(pick_actor);
      pick_actor->Delete();cout<<"Deleting pick_actor="<<pick_actor<<endl;
      pick_actor = NULL;
    }
    pick_actor = vtkActor::New();cout<<"New pick_actor="<<pick_actor<<endl;
    pick_actor->SetMapper(pick_mapper);
    pick_actor->GetProperty()->SetRepresentationToSurface();
    pick_actor->GetProperty()->SetColor(0,0,1);
    getRenderer()->AddActor(pick_actor);
    PickedPoint=nodeId;
    return(true);
  } else {
    return(false);
  }
}

bool GuiMainWindow::pickCell(vtkIdType cellId)
{
  if (cellId >= 0 && cellId<grid->GetNumberOfCells()) {
    vtkIdType *pts, Npts;
    grid->GetCellPoints(cellId, Npts, pts);
    vec3_t x(0,0,0);
    for (vtkIdType i = 0; i < Npts; ++i) {
      vec3_t xp;
      grid->GetPoints()->GetPoint(pts[i], xp.data());
      x += double(1)/Npts * xp;
    }
    pick_sphere->SetCenter(x.data());
    double R = 1e99;
    for (vtkIdType i = 0; i < Npts; ++i) {
      vec3_t xp;
      grid->GetPoints()->GetPoint(pts[i], xp.data());
      R = min(R, 0.25*(xp-x).abs());
    }
    m_ReferenceSize=R;//Used for text annotations too!
    pick_sphere->SetRadius(R);
    pick_mapper->SetInput(pick_sphere->GetOutput());
    
    if (pick_actor) {
      getRenderer()->RemoveActor(pick_actor);
      pick_actor->Delete();
      pick_actor = NULL;
    }
    pick_actor = vtkActor::New();
    pick_actor->SetMapper(pick_mapper);
    pick_actor->GetProperty()->SetRepresentationToSurface();
    pick_actor->GetProperty()->SetColor(1,0,0);
    getRenderer()->AddActor(pick_actor);
    PickedCell = cellId;
    return(true);
  } else {
    return(false);
  }
}

void GuiMainWindow::importSTL()
{
  StlReader stl;
  stl();
  updateBoundaryCodes(true);
  updateActors();
  updateStatusBar();
  zoomAll();
}

void GuiMainWindow::importGmsh1Ascii()
{
  GmshReader gmsh;
  gmsh.setV1Ascii();
  gmsh();
  updateBoundaryCodes(true);
  updateActors();
  updateStatusBar();
  zoomAll();
}

void GuiMainWindow::exportGmsh1Ascii()
{
  GmshWriter gmsh;
  gmsh.setV1Ascii();
  gmsh();
}

void GuiMainWindow::importGmsh2Ascii()
{
  GmshReader gmsh;
  gmsh.setV2Ascii();
  gmsh();
  updateBoundaryCodes(true);
  updateActors();
  updateStatusBar();
  zoomAll();
}

void GuiMainWindow::exportGmsh2Ascii()
{
  GmshWriter gmsh;
  gmsh.setV2Ascii();
  gmsh();
}

void GuiMainWindow::exportNeutral()
{
  NeutralWriter neutral;
  neutral();
}

void GuiMainWindow::zoomAll()
{
  getRenderer()->ResetCamera();
  getRenderWindow()->Render();
}

void GuiMainWindow::zoomOnPickedObject()
{
  if(pick_actor!=NULL)
  {
    getRenderer()->ResetCamera(pick_actor->GetBounds());
    getRenderWindow()->Render();
  }
}

void GuiMainWindow::deselectAll()
{
  cout<<"void GuiMainWindow::deselectAll()"<<endl;
  pick_actor->VisibilityOff();
  updateActors();
}

///@@@  TODO: Should display a window
void GuiMainWindow::info()
{
  ShowInfo info(ui.radioButton_CellPicker->isChecked(),PickedPoint,PickedCell);
  info();
}

int GuiMainWindow::quickSave()
{
  ///@@@ might be re-activated with RAM support
  /*
  if(grid->GetNumberOfPoints()>0)
  {
    current_operation++;
    QFileInfo fileinfo(current_filename);
    QString l_filename = m_LogDir + fileinfo.completeBaseName() + "_" + QString("%1").arg(current_operation);
    last_operation=current_operation;
    cout<<"Operation "<<current_operation<<endl;//" : Saving as l_filename="<<l_filename.toLatin1().data()<<endl;
    QuickSave(l_filename);
    setWindowTitle(current_filename + " - enGrid - " + QString("%1").arg(current_operation) );
    if(current_operation>0) ui.actionUndo->setEnabled(true);
    ui.actionRedo->setEnabled(false);
  }
  else cout<<"No grid to save!"<<endl;
  return(current_operation);
  */
  return 0;
}

void GuiMainWindow::quickLoad(int a_operation)
{
  a_operation = a_operation;
  ///@@@ might be re-activated with RAM support
  /*
  QFileInfo fileinfo(current_filename);
  QString l_filename = m_LogDir + fileinfo.completeBaseName() + "_" + QString("%1").arg(a_operation) + ".vtu";
  QuickLoad(l_filename);
  setWindowTitle(current_filename + " - enGrid - " + QString("%1").arg(current_operation) );
  */
}

void GuiMainWindow::Undo()
{
  QMessageBox::critical(this, "de-activated", "Undo is not doing anything at the moment!");
  /*
  cout << "Undoing operation " << current_operation << endl;
  current_operation--;
  QuickLoad(current_operation);
  ui.actionRedo->setEnabled(true);
  if(current_operation<=0) ui.actionUndo->setEnabled(false);
  */
}

void GuiMainWindow::Redo()
{
  QMessageBox::critical(this, "de-activated", "Redo is not doing anything at the moment!");
  /*
  current_operation++;
  cout << "Redoing operation " << current_operation << endl;
  QuickLoad(current_operation);
  ui.actionUndo->setEnabled(true);
  if(current_operation>=last_operation) ui.actionRedo->setEnabled(false);
  */
}

void GuiMainWindow::ResetOperationCounter()
{
  current_operation=-1;
  last_operation=current_operation;
  ui.actionUndo->setEnabled(false);
  ui.actionRedo->setEnabled(false);
}

void GuiMainWindow::openBC()
{
  openBC(current_filename);
}

void GuiMainWindow::saveBC()
{
  saveBC(current_filename);
}

void GuiMainWindow::openBC(QString a_file)
{
  QString bc_file = a_file + ".bcs";
  QFile file(bc_file);
  bcmap.clear();
  volmap.clear();
  if (file.exists()) {
    file.open(QIODevice::ReadOnly | QIODevice::Text);
    QTextStream f(&file);
    while (!f.atEnd()) {
      QString name, type;
      int i;
      f >> i >> name >> type;
      if (i >= 0) {
        bcmap[i] = BoundaryCondition(name,type);
      } else {
        VolumeDefinition V(name, -i);
        QString text = type.replace(",", " ").replace(":", " ");
        QTextStream s(&text);
        while (!s.atEnd()) {
          QString bc_txt, sign_txt;
          s >> bc_txt >> sign_txt;
          V.addBC(bc_txt.toInt(), sign_txt.toInt());
        }
        volmap[name] = V;
      }
    }
  }
}

void GuiMainWindow::saveBC(QString a_file)
{
  QString bc_file = a_file + ".bcs";
  QFile file(bc_file);
  file.open(QIODevice::WriteOnly | QIODevice::Text);
  QTextStream f(&file);
  foreach (int i, m_AllBoundaryCodes) {
    BoundaryCondition bc = bcmap[i];
    f << i << " " << bc.getName() << " " << bc.getType() << "\n";
  }
  foreach (VolumeDefinition V, volmap) {
    QString dirs = "";
    bool first = true;
    foreach (int i, m_AllBoundaryCodes) {
      BoundaryCondition bc = bcmap[i];
      if (!first) {
        dirs += ",";
      } else {
        first = false;
      }
      QString num;
      num.setNum(i);
      dirs += num + ":";
      num.setNum(V.getSign(i));
      dirs += num;
    }
    f << "-" << V.getVC() << " " << V.getName() << " " << dirs << "\n";
  }
}

///@@@  TODO: I think this should also be a done by a subclass of IOOperation just like for import operations
void GuiMainWindow::open()
{
  current_filename = QFileDialog::getOpenFileName(NULL, "open grid from file", getCwd(), "VTK unstructured grid files (*.vtu *.VTU)");
  if (!current_filename.isNull()) {
    GuiMainWindow::setCwd(QFileInfo(current_filename).absolutePath());
    EG_VTKSP(vtkXMLUnstructuredGridReader,vtu);
    vtu->SetFileName(current_filename.toAscii().data());
    vtu->Update();
    grid->DeepCopy(vtu->GetOutput());
    createBasicFields(grid, grid->GetNumberOfCells(), grid->GetNumberOfPoints());
    setWindowTitle(current_filename + " - enGrid - " + QString("%1").arg(current_operation) );
    openBC();
    updateBoundaryCodes(true);
    createIndices(grid);
    updateActors();
    updateStatusBar();
    zoomAll();
    ResetOperationCounter();
    quickSave();
  }
}

void GuiMainWindow::save()
{
  cout << current_filename.toAscii().data() << endl;
  if (current_filename == "untitled.vtu") {
    saveAs();
  } else {
    EG_VTKDCC(vtkDoubleArray, cell_VA, grid, "cell_VA");
    for (vtkIdType cellId = 0; cellId < grid->GetNumberOfCells(); ++cellId) {
      cell_VA->SetValue(cellId, GeometryTools::cellVA(grid, cellId, true));
    }
    EG_VTKSP(vtkXMLUnstructuredGridWriter,vtu);
    addVtkTypeInfo();
    createIndices(grid);
    vtu->SetFileName(current_filename.toAscii().data());
    vtu->SetDataModeToBinary();
    vtu->SetInput(grid);
    vtu->Write();
    if(vtu->GetErrorCode()) {
      QMessageBox::critical(this, tr("Save failed"), tr("The grid could not be saved as:\n%1").arg(current_filename));
    }
    else{
      saveBC();
    }
  }
}

void GuiMainWindow::saveAs()
{
  current_filename = QFileDialog::getSaveFileName
    (
      NULL,
      "write grid to file",
      getCwd(),
      "VTK unstructured grid files (*.vtu *.VTU)"
    );
  if (!current_filename.isNull()) {
    if (current_filename.right(4) != ".vtu") {
      if (current_filename.right(4) != ".VTU") {
        current_filename += ".vtu";
      }
    }
    EG_VTKDCC(vtkDoubleArray, cell_VA, grid, "cell_VA");
    for (vtkIdType cellId = 0; cellId < grid->GetNumberOfCells(); ++cellId) {
      cell_VA->SetValue(cellId, GeometryTools::cellVA(grid, cellId, true));
    }
    GuiMainWindow::setCwd(QFileInfo(current_filename).absolutePath());
    EG_VTKSP(vtkXMLUnstructuredGridWriter,vtu);
    addVtkTypeInfo();
    createIndices(grid);
    vtu->SetFileName(current_filename.toAscii().data());
    vtu->SetDataModeToBinary();
    vtu->SetInput(grid);
    vtu->Write();
    
    if(vtu->GetErrorCode()) {
      QMessageBox::critical(this, tr("Save failed"), tr("The grid could not be saved as:\n%1").arg(current_filename));
    }
    else{
      saveBC();
      //for the undo/redo operations
      setWindowTitle(current_filename + " - enGrid - " + QString("%1").arg(current_operation) );
      ResetOperationCounter();
      quickSave();
    }
  }
}

void GuiMainWindow::quickSave(QString a_filename)
{
  if(grid->GetNumberOfPoints()>0)
  {
    QFileInfo fileinfo(a_filename);
    cout<<"a_filename="<<a_filename.toLatin1().data()<<endl;
    cout<<"fileinfo.suffix()="<<fileinfo.suffix().toLatin1().data()<<endl;
    if(fileinfo.suffix()!="vtu") a_filename=a_filename + ".vtu";
    
    cout << "Saving as " << a_filename.toAscii().data() << endl;
    
    EG_VTKDCC(vtkDoubleArray, cell_VA, grid, "cell_VA");
    for (vtkIdType cellId = 0; cellId < grid->GetNumberOfCells(); ++cellId) {
      cell_VA->SetValue(cellId, GeometryTools::cellVA(grid, cellId, true));
    }
    EG_VTKSP(vtkXMLUnstructuredGridWriter,vtu);
    addVtkTypeInfo();
    createIndices(grid);
    vtu->SetFileName(a_filename.toAscii().data());
    vtu->SetDataModeToBinary();
    vtu->SetInput(grid);
    vtu->Write();
    if(vtu->GetErrorCode()) {
      QMessageBox::critical(this, tr("Save failed"), tr("The grid could not be saved as:\n%1").arg(current_filename));
    }
    else{
      saveBC(a_filename);
    }
  }
  else cout<<"No grid to save!"<<endl;
}

void GuiMainWindow::quickLoad(QString a_filename)
{
    cout << "Loading " << a_filename.toAscii().data() << endl;

    if (!a_filename.isNull()) {
      // GuiMainWindow::setCwd(QFileInfo(a_filename).absolutePath());
      EG_VTKSP(vtkXMLUnstructuredGridReader,vtu);
      vtu->SetFileName(a_filename.toAscii().data());
      vtu->Update();
      grid->DeepCopy(vtu->GetOutput());
      createBasicFields(grid, grid->GetNumberOfCells(), grid->GetNumberOfPoints());
      setWindowTitle(a_filename + " - enGrid - " + QString("%1").arg(current_operation) );
      openBC(a_filename);
      updateBoundaryCodes(true);
      updateActors();
      updateStatusBar();
      zoomAll();
    }
  
}

void GuiMainWindow::updateStatusBar()
{
  QString num, txt = "enGrid is currently busy with an operation ...";
  if (!busy) {
    txt = "";
  }
  if (!tryLock()) {
    status_label->setText(txt);
    ui.label_node_cell_info->setText(txt);
    return;
  }
  vtkIdType Ncells = grid->GetNumberOfCells();
  vtkIdType Nnodes = grid->GetNumberOfPoints();
  vtkIdType Ntris  = 0;
  vtkIdType Nquads = 0;
  vtkIdType Ntets  = 0;
  vtkIdType Npyras = 0;
  vtkIdType Nprism = 0;
  vtkIdType Nhexas = 0;
  for (vtkIdType i = 0; i < Ncells; ++i) {
    int ct = grid->GetCellType(i);
    if      (ct == VTK_TRIANGLE)   ++Ntris;
    else if (ct == VTK_QUAD)       ++Nquads;
    else if (ct == VTK_TETRA)      ++Ntets;
    else if (ct == VTK_WEDGE)      ++Nprism;
    else if (ct == VTK_PYRAMID)    ++Npyras;
    else if (ct == VTK_HEXAHEDRON) ++Nhexas;
  }
  num.setNum(Ntets + Npyras + Nprism + Nhexas); txt += num + " volume cells(";
  num.setNum(Ntets);  txt += num + " tetras, ";
  num.setNum(Npyras); txt += num + " pyramids, ";
  num.setNum(Nprism); txt += num + " prisms, ";
  num.setNum(Nhexas); txt += num + " hexas), ";
  num.setNum(Ntris + Nquads); txt += num + " surface cells(";
  num.setNum(Ntris);  txt += num + " triangles, ";
  num.setNum(Nquads); txt += num + " quads), ";
  num.setNum(Nnodes); txt += num + " nodes";
  
  if(ui.radioButton_CellPicker->isChecked())
  {
    QString pick_txt = ", picked cell: ";
    vtkIdType id_cell = PickedCell;
    if (id_cell < 0) {
      pick_txt += "no cell picked";
    } else {
      vtkIdType type_cell = grid->GetCellType(id_cell);
      if      (type_cell == VTK_TRIANGLE)   pick_txt += "tri";
      else if (type_cell == VTK_QUAD)       pick_txt += "qua";
      else if (type_cell == VTK_TETRA)      pick_txt += "tet";
      else if (type_cell == VTK_PYRAMID)    pick_txt += "pyr";
      else if (type_cell == VTK_WEDGE)      pick_txt += "pri";
      else if (type_cell == VTK_HEXAHEDRON) pick_txt += "hex";
      vtkIdType N_pts, *pts;
      grid->GetCellPoints(id_cell, N_pts, pts);
      pick_txt += " [";
      for (int i_pts = 0; i_pts < N_pts; ++i_pts) {
        QString num;
        num.setNum(pts[i_pts]);
        pick_txt += num;
        if (i_pts < N_pts-1) {
          pick_txt += ",";
        }
      }
      pick_txt += "]";
      QString tmp;
      EG_VTKDCC(vtkIntArray, cell_code, grid, "cell_code");
      tmp.setNum(cell_code->GetValue(id_cell));
      pick_txt += " code=" + tmp;
      tmp.setNum(id_cell);
      pick_txt += " id=" + tmp;
    }
    txt += pick_txt;
  }
  else
  {
    QString pick_txt = ", picked node: ";
    vtkIdType id_node = PickedPoint;
    if (id_node < 0) {
      pick_txt += "no node picked";
    } else {
      vec3_t x;
      grid->GetPoints()->GetPoint(id_node,x.data());
      pick_txt += " [";
      for (int i = 0; i < 3; i++) {
        QString num;
        num.setNum(x[i]);
        pick_txt += num;
        if (i < 2) {
          pick_txt += ",";
        }
      }
      pick_txt += "]";
      QString tmp;
      EG_VTKDCN(vtkDoubleArray, node_meshdensity_desired, grid, "node_meshdensity_desired");
      tmp.setNum(node_meshdensity_desired->GetValue(id_node));
      pick_txt += " wanted density=" + tmp;
      EG_VTKDCN(vtkDoubleArray, node_meshdensity_current, grid, "node_meshdensity_current");
      tmp.setNum(node_meshdensity_current->GetValue(id_node));
      pick_txt += " current density=" + tmp;
      EG_VTKDCN(vtkIntArray, node_specified_density, grid, "node_specified_density");
      tmp.setNum(node_specified_density->GetValue(id_node));
      pick_txt += " node_specified_density=" + tmp;
      EG_VTKDCN(vtkCharArray, node_type, grid, "node_type");
      pick_txt += " type=" + QString(VertexType2Str( node_type->GetValue(id_node)));
      tmp.setNum(id_node);
      pick_txt += " id_node=" + tmp;
    }
    
    txt += pick_txt;
  }
  
  ///@@@ TODO: Reduce size of text for small screens or better: allow making the window smaller than the text
  status_label->setText(txt);
  ui.label_node_cell_info->setText(txt);
  unlock();
}

void GuiMainWindow::selectBoundaryCodes()
{
  GuiSelectBoundaryCodes bcodes;
  bcodes.setDisplayBoundaryCodes(m_DisplayBoundaryCodes);
  bcodes.setBoundaryCodes(m_AllBoundaryCodes);
  bcodes();
  bcodes.getThread().wait();
  bcodes.getSelectedBoundaryCodes(m_DisplayBoundaryCodes);
  m_BCodesFilter->SetBoundaryCodes(m_DisplayBoundaryCodes);
  updateActors();
}

void GuiMainWindow::updateBoundaryCodes(bool all_on)
{
  try {
    m_AllBoundaryCodes.clear();
    EG_VTKDCC(vtkIntArray, cell_code, grid, "cell_code");
    for (vtkIdType i = 0; i < grid->GetNumberOfCells(); ++i) {
      int ct = grid->GetCellType(i);
      if ((ct == VTK_TRIANGLE) || (ct == VTK_QUAD)) {
        m_AllBoundaryCodes.insert(cell_code->GetValue(i));
      }
    }
    if (all_on) {
      m_DisplayBoundaryCodes.clear();
      foreach (int bc, m_AllBoundaryCodes) {
        m_DisplayBoundaryCodes.insert(bc);
      }
    } else {
      QSet<int> dbcs;
      foreach (int bc, m_DisplayBoundaryCodes) {
        if (m_AllBoundaryCodes.contains(bc)) {
          dbcs.insert(bc);
        }
      }
      m_DisplayBoundaryCodes.clear();
      foreach (int bc, m_AllBoundaryCodes) {
        if (dbcs.contains(bc)) {
          m_DisplayBoundaryCodes.insert(bc);
        }
      }
    }
    m_BCodesFilter->SetBoundaryCodes(m_DisplayBoundaryCodes);
  } catch (Error err) {
    err.display();
  }
}

void GuiMainWindow::normalExtrusion()
{
  GuiNormalExtrusion extr;
  extr();
  updateBoundaryCodes(false);
  updateActors();
}

void GuiMainWindow::setAxesVisibility()
{
  if (ui.actionViewAxes->isChecked()) {
    axes->SetVisibility(1);
  } else {
    axes->SetVisibility(0);
  }
  getRenderWindow()->Render();
}

void GuiMainWindow::setViewingMode()
{
  if (ui.actionViewOrthogonal->isChecked()) getRenderer()->GetActiveCamera()->ParallelProjectionOn();
  else getRenderer()->GetActiveCamera()->ParallelProjectionOff();
  getRenderWindow()->Render();
}

void GuiMainWindow::viewNodeIDs()
{
  int N=grid->GetNumberOfPoints();
  cout<<"N="<<N<<endl;
  if (ui.actionViewNodeIDs->isChecked()) {
    cout<<"Activating node ID view"<<endl;
    m_NodeTextVectorText.resize(N);
    m_NodeTextPolyDataMapper.resize(N);
    m_NodeTextFollower.resize(N);
    for(int i = 0; i < N; ++i){
      m_NodeTextVectorText[i]=vtkVectorText::New();
      QString tmp;
      tmp.setNum(i);
      m_NodeTextVectorText[i]->SetText(tmp.toLatin1().data());
      m_NodeTextPolyDataMapper[i]=vtkPolyDataMapper::New();
      m_NodeTextPolyDataMapper[i]->SetInputConnection(m_NodeTextVectorText[i]->GetOutputPort());
      m_NodeTextFollower[i]=vtkFollower::New();
      m_NodeTextFollower[i]->SetMapper(m_NodeTextPolyDataMapper[i]);
      m_NodeTextFollower[i]->SetScale(m_ReferenceSize, m_ReferenceSize, m_ReferenceSize);
      vec3_t M;
      grid->GetPoint(i, M.data());
      vec3_t tmp_M = M;
      vec3_t OffSet = m_ReferenceSize*tmp_M.normalise();
      M = M + OffSet;
      m_NodeTextFollower[i]->AddPosition(M[0], M[1], M[2]);
      m_NodeTextFollower[i]->SetCamera(getRenderer()->GetActiveCamera());
      m_NodeTextFollower[i]->GetProperty()->SetColor(0,0,1);
      getRenderer()->AddActor(m_NodeTextFollower[i]);
    }
  }
  else {
    cout<<"Deactivating node ID view"<<endl;
    for(unsigned int i = 0; i < m_NodeTextFollower.size();i++){
      getRenderer()->RemoveActor(m_NodeTextFollower[i]);
      m_NodeTextFollower[i]->Delete();
      m_NodeTextPolyDataMapper[i]->Delete();
      m_NodeTextVectorText[i]->Delete();
    }
    m_NodeTextFollower.clear();
    m_NodeTextPolyDataMapper.clear();
    m_NodeTextVectorText.clear();
  }
  
  getRenderWindow()->Render();
}

void GuiMainWindow::viewCellIDs()
{
  vtkIdType N=grid->GetNumberOfCells();
  cout<<"N="<<N<<endl;
  if (ui.actionViewCellIDs->isChecked()) {
    cout<<"Activating cell ID view"<<endl;
    m_CellTextVectorText.resize(N);
    m_CellTextPolyDataMapper.resize(N);
    m_CellTextFollower.resize(N);
    for (vtkIdType id_cell = 0; id_cell < N; ++id_cell){
      m_CellTextVectorText[id_cell] = vtkVectorText::New();
      
      QString tmp;
      
      if(ui.comboBox_CellTextField->currentIndex()==0) {
        tmp.setNum(id_cell);
      } else if (ui.comboBox_CellTextField->currentIndex()>0) {
        EG_VTKDCC(vtkIntArray, current_cell_field, grid, ui.comboBox_CellTextField->currentText().toLatin1().data());
        tmp.setNum(current_cell_field->GetValue(id_cell));
      }
      else EG_BUG;
      
      m_CellTextVectorText[id_cell]->SetText(tmp.toLatin1().data());
      m_CellTextPolyDataMapper[id_cell]=vtkPolyDataMapper::New();
      m_CellTextPolyDataMapper[id_cell]->SetInputConnection(m_CellTextVectorText[id_cell]->GetOutputPort());
      m_CellTextFollower[id_cell]=vtkFollower::New();
      m_CellTextFollower[id_cell]->SetMapper(m_CellTextPolyDataMapper[id_cell]);
      m_CellTextFollower[id_cell]->SetScale(m_ReferenceSize, m_ReferenceSize, m_ReferenceSize);
      vtkIdType N_pts,*pts;
      grid->GetCellPoints(id_cell,N_pts,pts);
      vec3_t Center(0,0,0);
      for (int p = 0; p < N_pts; ++p) {
        vec3_t M;
        grid->GetPoint(pts[p],M.data());
        Center+=M.data();
      }
      vec3_t OffSet = m_ReferenceSize*triNormal(grid, pts[0], pts[1], pts[2]).normalise();
      Center = 1.0/(double)N_pts*Center+OffSet;
      m_CellTextFollower[id_cell]->AddPosition(Center[0], Center[1], Center[2]);
      m_CellTextFollower[id_cell]->SetCamera(getRenderer()->GetActiveCamera());
      m_CellTextFollower[id_cell]->GetProperty()->SetColor(1, 0, 0);
      getRenderer()->AddActor(m_CellTextFollower[id_cell]);
    }
  } else {
    cout<<"Deactivating cell ID view"<<endl;
    for (vtkIdType id_cell = 0; id_cell < (vtkIdType) m_CellTextFollower.size(); ++id_cell) {
      getRenderer()->RemoveActor(m_CellTextFollower[id_cell]);
      m_CellTextFollower[id_cell]->Delete();
      m_CellTextPolyDataMapper[id_cell]->Delete();
      m_CellTextVectorText[id_cell]->Delete();
    }
    m_CellTextFollower.clear();
    m_CellTextPolyDataMapper.clear();
    m_CellTextVectorText.clear();
  }
  
  getRenderWindow()->Render();
}

void GuiMainWindow::addVtkTypeInfo()
{
  EG_VTKSP(vtkIntArray, vtk_type);
  vtk_type->SetName("vtk_type");
  vtk_type->SetNumberOfValues(grid->GetNumberOfCells());
  for (vtkIdType cellId = 0; cellId < grid->GetNumberOfCells(); ++cellId) {
    vtk_type->SetValue(cellId, grid->GetCellType(cellId));
  }
  grid->GetCellData()->AddArray(vtk_type);
}

void GuiMainWindow::pickCellCallBack
(
  vtkObject *caller, 
  unsigned long int eid, 
  void *clientdata, 
  void *calldata
)
{
  caller = caller;
  eid = eid;
  clientdata = clientdata;
  calldata = calldata;
  THIS->updateActors();
  THIS->updateStatusBar();
  cout<<"pickCellCallBack"<<endl;
}

void GuiMainWindow::pickPointCallBack
(
  vtkObject *caller, 
  unsigned long int eid, 
  void *clientdata, 
  void *calldata
)
{
  caller = caller;
  eid = eid;
  clientdata = clientdata;
  calldata = calldata;
  THIS->updateActors();
  THIS->updateStatusBar();
  cout<<"pickPointCallBack"<<endl;
}

vtkIdType GuiMainWindow::getPickedCell()
{
  vtkIdType picked_cell = -1;
  if (THIS->grid->GetNumberOfCells() > 0) {
    THIS->m_BCodesFilter->Update();
    EG_VTKDCC(vtkLongArray_t, cell_index, THIS->m_BCodesFilter->GetOutput(), "cell_index");
    
    vtkIdType cellId;
    if(m_UseVTKInteractor) cellId = THIS->CellPicker->GetCellId();
    else cellId = PickedCell;
    
    if (cellId >= 0) {
      if (cellId < THIS->m_BCodesFilter->GetOutput()->GetNumberOfCells()) {
        picked_cell = cell_index->GetValue(cellId);
      }
    }
  }
  return picked_cell;
}

vtkIdType GuiMainWindow::getPickedPoint()
{
  vtkIdType picked_point = -1;
  if (THIS->grid->GetNumberOfCells() > 0) {
    THIS->m_BCodesFilter->Update();
    
    vtkIdType pointId;
    if(m_UseVTKInteractor) pointId = THIS->PointPicker->GetPointId();
    else pointId = PickedPoint;
    
    if (pointId >= 0) {
      picked_point = pointId;
    }
  }
  return picked_point;
}

void GuiMainWindow::changeSurfaceOrientation()
{
  for (vtkIdType cellId = 0; cellId < grid->GetNumberOfCells(); ++cellId) {
    vtkIdType Npts, *pts;
    grid->GetCellPoints(cellId, Npts, pts);
    QVector<vtkIdType> nodes(Npts);
    for (vtkIdType j = 0; j < Npts; ++j) nodes[j]          = pts[j];
    for (vtkIdType j = 0; j < Npts; ++j) pts[Npts - j - 1] = nodes[j];
  }
  updateActors();
}

void GuiMainWindow::checkSurfaceOrientation()
{
  CorrectSurfaceOrientation corr_surf;
  vtkIdType picked_cell = getPickedCell();
  if (picked_cell >= 0) {
    corr_surf.setStart(picked_cell);
  }
  corr_surf();
  updateActors();
}

void GuiMainWindow::improveAspectRatio()
{
  GuiImproveAspectRatio impr_ar;
  impr_ar();
  updateActors();
}

void GuiMainWindow::exportAsciiStl()
{
  StlWriter stl;
  stl();
}

void GuiMainWindow::exportBinaryStl()
{
}

void GuiMainWindow::periodicUpdate()
{
  Operation::collectGarbage(); 
  updateStatusBar();
}

void GuiMainWindow::viewXP()
{
  getRenderer()->ResetCamera();
  double x[3];
  getRenderer()->GetActiveCamera()->GetFocalPoint(x);
  x[0] += 1;
  getRenderer()->GetActiveCamera()->SetPosition(x);
  getRenderer()->GetActiveCamera()->ComputeViewPlaneNormal();
  getRenderer()->GetActiveCamera()->SetViewUp(0,1,0);
  getRenderer()->ResetCamera();
  getRenderWindow()->Render();
}

void GuiMainWindow::viewXM()
{
  getRenderer()->ResetCamera();
  double x[3];
  getRenderer()->GetActiveCamera()->GetFocalPoint(x);
  x[0] -= 1;
  getRenderer()->GetActiveCamera()->SetPosition(x);
  getRenderer()->GetActiveCamera()->ComputeViewPlaneNormal();
  getRenderer()->GetActiveCamera()->SetViewUp(0,1,0);
  getRenderer()->ResetCamera();
  getRenderWindow()->Render();
}

void GuiMainWindow::viewYP()
{
  getRenderer()->ResetCamera();
  double x[3];
  getRenderer()->GetActiveCamera()->GetFocalPoint(x);
  x[1] += 1;
  getRenderer()->GetActiveCamera()->SetPosition(x);
  getRenderer()->GetActiveCamera()->ComputeViewPlaneNormal();
  getRenderer()->GetActiveCamera()->SetViewUp(0,0,-1);
  getRenderer()->ResetCamera();
  getRenderWindow()->Render();
}

void GuiMainWindow::viewYM()
{
  getRenderer()->ResetCamera();
  double x[3];
  getRenderer()->GetActiveCamera()->GetFocalPoint(x);
  x[1] -= 1;
  getRenderer()->GetActiveCamera()->SetPosition(x);
  getRenderer()->GetActiveCamera()->ComputeViewPlaneNormal();
  getRenderer()->GetActiveCamera()->SetViewUp(0,0,-1);
  getRenderer()->ResetCamera();
  getRenderWindow()->Render();
}

void GuiMainWindow::viewZP()
{
  getRenderer()->ResetCamera();
  double x[3];
  getRenderer()->GetActiveCamera()->GetFocalPoint(x);
  x[2] += 1;
  getRenderer()->GetActiveCamera()->SetPosition(x);
  getRenderer()->GetActiveCamera()->ComputeViewPlaneNormal();
  getRenderer()->GetActiveCamera()->SetViewUp(0,1,0);
  getRenderer()->ResetCamera();
  getRenderWindow()->Render();
}

void GuiMainWindow::viewZM()
{
  getRenderer()->ResetCamera();
  double x[3];
  getRenderer()->GetActiveCamera()->GetFocalPoint(x);
  x[2] -= 1;
  getRenderer()->GetActiveCamera()->SetPosition(x);
  getRenderer()->GetActiveCamera()->ComputeViewPlaneNormal();
  getRenderer()->GetActiveCamera()->SetViewUp(0,1,0);
  getRenderer()->ResetCamera();
  getRenderWindow()->Render();
}

void GuiMainWindow::callFixSTL()
{
  FixSTL *fix;
  fix = new FixSTL();
  fix->setGui();
  (*fix)();
  updateBoundaryCodes(false);
  updateActors();
}

void GuiMainWindow::editBoundaryConditions()
{
  GuiEditBoundaryConditions editbcs;
  editbcs.setBoundaryCodes(m_AllBoundaryCodes);
  editbcs.setMap(&bcmap);
  editbcs();
}

void GuiMainWindow::configure()
{
  {
    // Just to create initial entries in the settings file 
    // so that the options menu isn't empty at first start.
    GridSmoother tmp01;
    GuiCreateBoundaryLayer tmp02;
  }
  GuiSettingsViewer settings(&qset);
  settings.exec();
}

void GuiMainWindow::about()
{
  QMessageBox box(this);
  
  QString title="ENGRID";
  QString version = QString("version ") + ENGRID_VERSION;
  
  version += " built on ";
  version += QString(__DATE__);
  version += " at ";
  version += QString(__TIME__);
  
  QString address = tr("ENGRID is being developed and maintained by:<br/>"
                       "enGits GmbH<br/>"
                       "Marie-Curie-Strasse 8<br/>"
                       "79539 Loerrach<br/>"
                       "Germany<br/>");
  
  QString mainurl="<a href=\"http://www.engits.com\">www.engits.com</a>";
  QString mail="<a href=\"mailto:info@engits.com\">info@engits.com</a>";
  QString gnuurl="<a href=\"http://www.gnu.org/licenses\">http://www.gnu.org/licenses</a>";
  QString license=tr("ENGRID is licenced under the GPL version 3.<br/>"
                     "(see ")+gnuurl+tr(" for details)<br/>");
  QString bugurl="<a href=\"http://sourceforge.net/tracker2/?func=add&group_id=245110&atid=1126548\">the bugtracker available on Sourceforge</a>";
  QString bugreporting=tr("To submit a bug report, please use ")+bugurl;
  box.setText(QString::fromLatin1("<center><img src=\":/icons/resources/icons/G.png\">"
                                  "<h3>%1</h3>"
                                  "<p>%2</p>"
                                  "<p>%3</p>"
                                  "<p>Homepage: %4</p>"
                                  "<p>E-mail: %5</p>"
                                  "<p>%6</p>"
                                  "<p>%7</p></center>")
              .arg(title).arg(version).arg(address).arg(mainurl).arg(mail).arg(license).arg(bugreporting));
  box.setWindowTitle(tr("about ENGRID"));
  box.setIcon(QMessageBox::NoIcon);
  box.exec();
  
}

void GuiMainWindow::getAllBoundaryCodes(QSet<int> &bcs)
{
  bcs.clear();
  foreach (int bc, m_AllBoundaryCodes) {
    bcs.insert(bc);
  }
}

void GuiMainWindow::getDisplayBoundaryCodes(QSet<int> &bcs)
{
  bcs.clear();
  foreach (int bc, m_DisplayBoundaryCodes) {
    bcs.insert(bc);
  }
}

QList<VolumeDefinition> GuiMainWindow::getAllVols()
{
  QList<VolumeDefinition> vols;
  foreach(VolumeDefinition vol, volmap) {
    vols.push_back(vol);
  }
  return vols;
}

void GuiMainWindow::setAllVols(QList<VolumeDefinition> vols)
{
  volmap.clear();
  foreach (VolumeDefinition V, vols) {
    volmap[V.getName()] = V;
  }
}

void GuiMainWindow::createDefaultVol()
{
  QList<VolumeDefinition> vols = getAllVols();
  if (vols.size() == 0) {
    VolumeDefinition V("default", 1);
    QSet<int> bcs;
    getAllBoundaryCodes(bcs);
    foreach (int bc, bcs) {
      V.addBC(bc, 1);
    }
    vols.append(V);
    setAllVols(vols);
  }
}

QString GuiMainWindow::getFilePath()
{
  QFileInfo fileinfo(current_filename);
  return fileinfo.absolutePath()+"/";
}

void GuiMainWindow::markOutputLine()
{
  cout << "\n****************************************\n";
  cout << QTime::currentTime().toString("hh:mm:ss").toAscii().data();
  cout << "\n****************************************\n" << endl;
}
