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
#include "plywriter.h"
#include "correctsurfaceorientation.h"
#include "guieditboundaryconditions.h"
#include "laplacesmoother.h"
#include "swaptriangles.h"

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
#include <QFileInfo>
#include <stdlib.h>
#include <stdio.h>

#include "geometrytools.h"

using namespace GeometryTools;

#include "guisettingsviewer.h"
#include "guitransform.h"
#include "egvtkinteractorstyle.h"
#include "showinfo.h"

QString GuiMainWindow::m_cwd = ".";
QSettings GuiMainWindow::m_qset("enGits","enGrid");
GuiMainWindow* GuiMainWindow::THIS = NULL;
QMutex GuiMainWindow::m_Mutex;
bool GuiMainWindow::m_UnSaved = true;

GuiMainWindow::GuiMainWindow() : QMainWindow(NULL)
{
  setupGuiMainWindow();
  if(m_open_last) {
    if(m_qset.contains("LatestFile")) {
//       qDebug()<<"Opening latest";
      open(m_qset.value("LatestFile").toString());
    }
  }
}

GuiMainWindow::GuiMainWindow(QString file_name) : QMainWindow(NULL)
{
  setupGuiMainWindow();
  open(file_name);
}

void GuiMainWindow::setupGuiMainWindow()
{
  ui.setupUi(this);
  THIS = this;
  
  // restore window size
  if(m_qset.contains("GuiMainWindow")) {
    setGeometry(m_qset.value("GuiMainWindow").toRect());
  }
  else {
    this->setWindowState(Qt::WindowMaximized);
  }
  
  // restore dockwidget positions
  if(m_qset.contains("dockWidget_states")) {
    restoreState(m_qset.value("dockWidget_states").toByteArray());
  }
  else {
    tabifyDockWidget(ui.dockWidget_output, ui.dockWidget_node_cell_info);
    tabifyDockWidget(ui.dockWidget_DisplayOptions, ui.dockWidget_DebuggingUtilities);
    ui.dockWidget_node_cell_info->hide();
    ui.dockWidget_DebuggingUtilities->hide();
  }
  
# include "std_connections.h"
  
  if (m_qset.contains("working_directory")) {
    m_cwd = m_qset.value("working_directory").toString();
  }

  setupVtk();

  resetOperationCounter();//clears undo/redo list and disables undo/redo
  m_CurrentFilename = "untitled.egc";
  setWindowTitle(m_CurrentFilename + " - enGrid - " + QString("%1").arg(m_CurrentOperation) );
  setUnsaved(true);
  
  m_StatusLabel = new QLabel(this);
  statusBar()->addWidget(m_StatusLabel);
  
  QString txt = "0 volume cells (0 tetras, 0 hexas, 0 pyramids, 0 prisms), ";
  txt += "0 surface cells (0 triangles, 0 quads), 0 nodes";
  m_StatusLabel->setText(txt);
  ui.label_node_cell_info->setText(txt);

  QString user = QString(getenv("USER"));
  QString basename="enGrid_output.txt";
  
  // define temporary path
  QDir dir("/");
  if (m_qset.contains("tmp_directory")) {
    m_LogDir = m_qset.value("tmp_directory").toString();
  } else {
    m_LogDir = dir.tempPath();
  }
  QDateTime now = QDateTime::currentDateTime();
  m_LogDir = m_LogDir + "/" + "enGrid_" + QDateTime::currentDateTime().toString("yyyyMMddhhmmsszzz") + "/";
  dir.mkpath(m_LogDir);
  
  m_LogFileName = m_LogDir + basename;
  cout << "m_LogFileName = " << qPrintable(m_LogFileName) << endl;

  m_SystemStdout = stdout;
  if(freopen (qPrintable(m_LogFileName), "w", stdout)==NULL) EG_BUG;
  
  m_Busy = false;
  
  setPickMode(true,true);
  m_PickedPoint = -1;
  m_PickedCell = -1;
  
  updateStatusBar();
  
  connect(&m_GarbageTimer, SIGNAL(timeout()), this, SLOT(periodicUpdate()));
  m_GarbageTimer.start(1000);
  
  connect(&m_LogTimer, SIGNAL(timeout()), this, SLOT(updateOutput()));
  m_LogTimer.start(1000);
  
  m_N_chars = 0;
  
  bool exp_features=false;
  getSet("General","enable experimental features",false,exp_features);
  getSet("General","enable undo+redo",false,m_undo_redo_enabled);
  bool undo_redo_mode;
  getSet("General","use RAM for undo+redo operations",false,undo_redo_mode);
  getSet("General", "open last used file on startup", false, m_open_last);
  
  ui.actionFoamWriter->setEnabled(exp_features);
  
  m_ReferenceSize=0.2;
  
  ui.doubleSpinBox_HueMin->setValue(0.667);
  ui.doubleSpinBox_HueMax->setValue(0);
  
  egvtkInteractorStyle *style = egvtkInteractorStyle::New();
  getInteractor()->SetInteractorStyle(style);
  style->Delete();

  // initialise XML document
  m_XmlHandler = new XmlHandler("engridcase");
//   this->resetXmlDoc();
  
  m_SolverIndex = 0;
  
  readRecentFiles();
  
}
//end of GuiMainWindow::GuiMainWindow() : QMainWindow(NULL)

void GuiMainWindow::resetXmlDoc()
{
  m_XmlHandler->resetXmlDoc();
/*  m_XmlDoc.clear();
  QDomElement root = m_XmlDoc.createElement("engridcase");
  m_XmlDoc.appendChild(root);*/
}

GuiMainWindow::~GuiMainWindow()
{
  writeRecentFiles();
  
  m_qset.setValue("GuiMainWindow", this->geometry());
  m_qset.setValue("dockWidget_states", this->saveState());
  
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
  
  delete m_XmlHandler;
}

void GuiMainWindow::setupVtk()
{
  m_Grid = vtkUnstructuredGrid::New();
  m_Renderer = vtkRenderer::New();
  getRenderWindow()->AddRenderer(m_Renderer);

  // coordinate axes
  m_Axes = vtkCubeAxesActor2D::New();
  //
  m_Axes->SetCamera(getRenderer()->GetActiveCamera());
  getRenderer()->AddActor(m_Axes);
  m_Axes->SetVisibility(0);

  // surface pipelines
  m_BackfaceProperty  = vtkProperty::New();
  m_SurfaceFilter     = vtkGeometryFilter::New();
  m_SurfaceMapper     = vtkPolyDataMapper::New();
  m_SurfaceWireMapper = vtkPolyDataMapper::New();
  m_BCodesFilter      = vtkEgBoundaryCodesFilter::New();
  m_LookupTable       = vtkLookupTable::New();
  m_SurfaceActor      = vtkActor::New();
  m_SurfaceWireActor  = vtkActor::New();
  m_LegendActor       = vtkScalarBarActor::New();
  //
  m_BCodesFilter->SetBoundaryCodes(m_DisplayBoundaryCodes);
  m_BCodesFilter->SetInput(m_Grid);
  m_SurfaceFilter->SetInput(m_BCodesFilter->GetOutput());
  m_SurfaceMapper->SetInput(m_SurfaceFilter->GetOutput());
  m_SurfaceWireMapper->SetInput(m_SurfaceFilter->GetOutput());
  m_SurfaceMapper->SetLookupTable(m_LookupTable);
  m_SurfaceActor->GetProperty()->SetRepresentationToSurface();
  m_SurfaceActor->GetProperty()->SetColor(0.5,1,0.5);
  m_SurfaceActor->SetBackfaceProperty(m_BackfaceProperty);
  m_SurfaceActor->GetBackfaceProperty()->SetColor(1,1,0.5);
  m_SurfaceActor->SetMapper(m_SurfaceMapper);
  getRenderer()->AddActor(m_SurfaceActor);
  m_SurfaceActor->SetVisibility(1);
  m_LegendActor->SetLookupTable(m_LookupTable);
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
  m_TetraGeometry = vtkGeometryFilter::New();
  m_TetraMapper  = vtkPolyDataMapper::New();
  //
  m_ExtrTetras->SetInput(m_Grid);
  m_ExtrTetras->SetAllOff();
  m_ExtrTetras->SetTetrasOn();;
  m_TetraGeometry->SetInput(m_ExtrTetras->GetOutput());
  m_TetraMapper->SetInput(m_TetraGeometry->GetOutput());
  m_TetraActor->SetMapper(m_TetraMapper);
  m_TetraActor->GetProperty()->SetColor(1,0,0);
  getRenderer()->AddActor(m_TetraActor);
  m_TetraActor->SetVisibility(0);

  // pyramid pipeline
  m_PyramidActor   = vtkActor::New();
  m_ExtrPyramids   = vtkEgExtractVolumeCells::New();
  m_PyramidGeometry = vtkGeometryFilter::New();
  m_PyramidMapper  = vtkPolyDataMapper::New();
  //
  m_ExtrPyramids->SetInput(m_Grid);
  m_ExtrPyramids->SetAllOff();
  m_ExtrPyramids->SetPyramidsOn();
  m_PyramidGeometry->SetInput(m_ExtrPyramids->GetOutput());
  m_PyramidMapper->SetInput(m_PyramidGeometry->GetOutput());
  m_PyramidActor->SetMapper(m_PyramidMapper);
  m_PyramidActor->GetProperty()->SetColor(1,1,0);
  getRenderer()->AddActor(m_PyramidActor);
  m_PyramidActor->SetVisibility(0);

  // wedge pipeline
  m_WedgeActor   = vtkActor::New();
  m_ExtrWedges   = vtkEgExtractVolumeCells::New();
  m_WedgeGeometry = vtkGeometryFilter::New();
  m_WedgeMapper  = vtkPolyDataMapper::New();
  //
  m_ExtrWedges->SetInput(m_Grid);
  m_ExtrWedges->SetAllOff();
  m_ExtrWedges->SetWedgesOn();
  m_WedgeGeometry->SetInput(m_ExtrWedges->GetOutput());
  m_WedgeMapper->SetInput(m_WedgeGeometry->GetOutput());
  m_WedgeActor->SetMapper(m_WedgeMapper);
  m_WedgeActor->GetProperty()->SetColor(0,1,0);
  getRenderer()->AddActor(m_WedgeActor);
  m_WedgeActor->SetVisibility(0);

  // hexa pipeline
  m_HexaActor   = vtkActor::New();
  m_ExtrHexes   = vtkEgExtractVolumeCells::New();
  m_HexaGeometry = vtkGeometryFilter::New();
  m_HexaMapper  = vtkPolyDataMapper::New();
  //
  m_ExtrHexes->SetInput(m_Grid);
  m_ExtrHexes->SetAllOff();
  m_ExtrHexes->SetHexesOn();
  m_HexaGeometry->SetInput(m_ExtrHexes->GetOutput());
  m_HexaMapper->SetInput(m_HexaGeometry->GetOutput());
  m_HexaActor->SetMapper(m_HexaMapper);
  m_HexaActor->GetProperty()->SetColor(0,0.7,1);
  getRenderer()->AddActor(m_HexaActor);
  m_HexaActor->SetVisibility(0);

  // volume wire pipeline
  m_VolumeWireActor  = vtkActor::New();
  m_ExtrVol          = vtkEgExtractVolumeCells::New();
  m_VolumeGeometry    = vtkGeometryFilter::New();
  m_VolumeWireMapper = vtkPolyDataMapper::New();
  //
  m_ExtrVol->SetInput(m_Grid);
  m_ExtrVol->SetAllOn();
  m_VolumeGeometry->SetInput(m_ExtrVol->GetOutput());
  m_VolumeWireMapper->SetInput(m_VolumeGeometry->GetOutput());
  m_VolumeWireActor->SetMapper(m_VolumeWireMapper);
  m_VolumeWireActor->GetProperty()->SetRepresentationToWireframe();
  m_VolumeWireActor->GetProperty()->SetColor(0,0,1);
  getRenderer()->AddActor(m_VolumeWireActor);
  m_VolumeWireActor->SetVisibility(0);

  // picker stuff
  m_PickSphere  = vtkSphereSource::New();
  m_PickMapper  = vtkPolyDataMapper::New();
  m_PickActor   = vtkActor::New();
  m_CellPicker  = vtkCellPicker::New();
  m_PointPicker = vtkPointPicker::New();
  
  m_PickSphere->SetRadius(0.25); //in case the user starts picking points instead of cells
  m_PickMapper->SetInput(m_PickSphere->GetOutput());
  m_PickActor->SetMapper(m_PickMapper);
  m_PickActor->GetProperty()->SetRepresentationToSurface();
  m_PickActor->GetProperty()->SetColor(0,0,1);
  m_PickActor->VisibilityOff();
  getRenderer()->AddActor(m_PickActor);
  
  vtkCallbackCommand *cbc = vtkCallbackCommand::New();
  cbc->SetCallback(pickCallBack);
  
  m_CellPicker->AddObserver(vtkCommand::EndPickEvent, cbc);
  m_PointPicker->AddObserver(vtkCommand::EndPickEvent, cbc);
  m_PickedObject = 0;
//   cbc->Delete();
}

void GuiMainWindow::updateOutput()
{
  QFile log_file(m_LogFileName);
  log_file.open(QIODevice::ReadOnly);
  QByteArray buffer = log_file.readAll();
  if (buffer.size() > m_N_chars) {
    QByteArray newchars = buffer.right(buffer.size() - m_N_chars);
    m_N_chars = buffer.size();
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
  return m_Renderer;
}

QVTKInteractor* GuiMainWindow::getInteractor()
{
  return ui.qvtkWidget->GetInteractor();
}

QString GuiMainWindow::getCwd()
{
  return m_cwd;
}

void GuiMainWindow::setCwd(QString dir)
{
  m_cwd = dir;
  m_qset.setValue("working_directory",dir);
}

void GuiMainWindow::setUnsaved(bool unsaved)
{
  m_UnSaved = unsaved;
}

void GuiMainWindow::scaleToData()
{
  int current_field=ui.comboBox_Field->currentIndex();
  if(current_field>0)
  {
    double range[2];

    m_SurfaceFilter->GetOutput()->GetPointData()->GetArray(current_field-1)->GetRange(range);
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
      m_SurfaceFilter->Update();
    }

    // fill node field combobox
    int current_field=ui.comboBox_Field->currentIndex();
    ui.comboBox_Field->clear();
    ui.comboBox_Field->addItem("None");
    for (int i = 0; i < m_Grid->GetPointData()->GetNumberOfArrays(); ++i) {
      ui.comboBox_Field->addItem(m_Grid->GetPointData()->GetArrayName(i));
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
    for (int i = 0; i < m_SurfaceFilter->GetOutput()->GetCellData()->GetNumberOfArrays(); ++i) {
      ui.comboBox_CellTextField->addItem(m_Grid->GetCellData()->GetArrayName(i));
    }
    if(current_cell_field == -1) {
      ui.comboBox_CellTextField->setCurrentIndex(0);
    } else {
      ui.comboBox_CellTextField->setCurrentIndex(current_cell_field);
    }
    current_field = ui.comboBox_Field->currentIndex();
    if(current_field > 0) {
      double range[2];
      m_SurfaceFilter->GetOutput()->GetPointData()->GetArray(current_field-1)->GetRange(range);
      ui.doubleSpinBox_FieldMin->setRange(range[0],range[1]);
      ui.doubleSpinBox_FieldMax->setRange(range[0],range[1]);
    }

    if(ui.comboBox_Field->currentIndex() > 0) {
      m_SurfaceMapper->SetColorModeToMapScalars();
      m_LookupTable->SetNumberOfColors(ui.spinBox_Color->value());
      m_LookupTable->SetHueRange(ui.doubleSpinBox_HueMin->value(),ui.doubleSpinBox_HueMax->value());
      m_LookupTable->Build();
      m_SurfaceMapper->SetScalarModeToUsePointFieldData();
      m_SurfaceMapper->ColorByArrayComponent(qPrintable(ui.comboBox_Field->currentText()),0);
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
      m_LegendActor->SetVisibility(0);
    }
    if (forced) {
      m_BCodesFilter->Update();
    }
    if(ui.checkBox_ShowPickSphere->checkState()) {
      if(m_UseVTKInteractor) {
        if(ui.radioButton_CellPicker->isChecked()) {
          getInteractor()->SetPicker(m_CellPicker);
          vtkIdType id_cell = getPickedCell();
          pickCell(id_cell);
        } else {
          getInteractor()->SetPicker(m_PointPicker);
          vtkIdType id_node = getPickedPoint();
          pickPoint(id_node);
        }
      } else {
        if (ui.radioButton_CellPicker->isChecked()) {
          pickCell(m_PickedCell);
        } else {
          pickPoint(m_PickedPoint);
        }
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
        m_TetraGeometry->Update();
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
        m_PyramidGeometry->Update();
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
        m_WedgeGeometry->Update();
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
        m_HexaGeometry->Update();
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
      m_VolumeGeometry->Update();
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
//   qDebug()<<"QApplication::setOverrideCursor(QCursor(Qt::WaitCursor)); called()";
  QApplication::setOverrideCursor(QCursor(Qt::WaitCursor));
  
  //if (!tryLock()) return;
  try {
    m_Axes->SetInput(m_Grid);
    updateSurfaceActors(forced);
    updateVolumeActors(forced);
    updateStatusBar();
  } catch (Error err) {
    err.display();
  }
  //unlock();
  
//   qDebug()<<"QApplication::restoreOverrideCursor(); called()";
  QApplication::restoreOverrideCursor();
}



void GuiMainWindow::forceUpdateActors()
{
//   qDebug()<<"void GuiMainWindow::forceUpdateActors() START";
  updateActors(true);
  getRenderWindow()->Render();
//   qDebug()<<"void GuiMainWindow::forceUpdateActors() END";
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

bool GuiMainWindow::pickPoint(vtkIdType id_node)
{
  if ((id_node >= 0) && (id_node < m_Grid->GetNumberOfPoints())) {
    vec3_t x(0,0,0);
    m_Grid->GetPoints()->GetPoint(id_node, x.data());
    m_PickSphere->SetCenter(x.data());
    m_PickedPoint = id_node;
    m_PickActor->GetProperty()->SetColor(0,0,1);
    m_PickActor->VisibilityOn();
    m_PickedObject = 1;
    return(true);
  } else {
    m_PickActor->VisibilityOff();
    m_PickedObject = 0;
    return(false);
  }
}

bool GuiMainWindow::pickCell(vtkIdType id_cell)
{
  if ((id_cell >= 0) && (id_cell < m_Grid->GetNumberOfCells())) {
    vtkIdType *pts, Npts;
    m_Grid->GetCellPoints(id_cell, Npts, pts);
    vec3_t x(0,0,0);
    for (vtkIdType i = 0; i < Npts; ++i) {
      vec3_t xp;
      m_Grid->GetPoints()->GetPoint(pts[i], xp.data());
      x += double(1)/Npts * xp;
    }
    m_PickSphere->SetCenter(x.data());
    double R = 1e99;
    for (vtkIdType i = 0; i < Npts; ++i) {
      vec3_t xp;
      m_Grid->GetPoints()->GetPoint(pts[i], xp.data());
      R = min(R, 0.25*(xp-x).abs());
    }
    m_ReferenceSize = R; //Used for text annotations too!
    m_PickSphere->SetRadius(R);
    m_PickedCell = id_cell;
    m_PickActor->GetProperty()->SetColor(1,0,0);
    m_PickActor->VisibilityOn();
    m_PickedObject = 2;
    return(true);
  } else {
    m_PickActor->VisibilityOff();
    m_PickedObject = 0;
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
  if(m_PickActor->GetVisibility()) {
    getRenderer()->ResetCamera(m_PickActor->GetBounds());
    getRenderWindow()->Render();
  }
}

void GuiMainWindow::deselectAll()
{
  cout << "void GuiMainWindow::deselectAll()" << endl;
  m_PickActor->VisibilityOff();
  updateActors();
}

///\todo Should display a window
void GuiMainWindow::info()
{
  ShowInfo info(ui.radioButton_CellPicker->isChecked(), m_PickedPoint, m_PickedCell);
  info();
}

int GuiMainWindow::quickSave()
{
  ///\todo add RAM support
  if(m_undo_redo_enabled) {
    if(m_Grid->GetNumberOfPoints()>0)
    {
      m_CurrentOperation++;
      QFileInfo fileinfo(m_CurrentFilename);
      QString l_filename = m_LogDir + fileinfo.completeBaseName() + "_" + QString("%1").arg(m_CurrentOperation);
      m_LastOperation=m_CurrentOperation;
      cout<<"Operation "<<m_CurrentOperation<<endl;
      saveAs(l_filename, false);
      if(m_CurrentOperation>0) ui.actionUndo->setEnabled(true);
      ui.actionRedo->setEnabled(false);
    }
    else cout<<"No grid to save!"<<endl;
    return(m_CurrentOperation);
  }
  return 0;
}

void GuiMainWindow::quickLoad(int a_operation)
{
  ///\todo add RAM support
  if(m_undo_redo_enabled) {
    QFileInfo fileinfo(m_CurrentFilename);
    QString l_filename = m_LogDir + fileinfo.completeBaseName() + "_" + QString("%1").arg(a_operation) + ".egc";
    open(l_filename, false);
  }
}

void GuiMainWindow::undo()
{
  if(m_undo_redo_enabled) {
    cout << "Undoing operation " << m_CurrentOperation << endl;
    m_CurrentOperation--;
    quickLoad(m_CurrentOperation);
    ui.actionRedo->setEnabled(true);
    if(m_CurrentOperation<=0) ui.actionUndo->setEnabled(false);
  }
  else {
    resetOperationCounter();
    QMessageBox::critical(this, "de-activated", "Undo is not doing anything at the moment!");
  }
}

void GuiMainWindow::redo()
{
  if(m_undo_redo_enabled) {
    m_CurrentOperation++;
    cout << "Redoing operation " << m_CurrentOperation << endl;
    quickLoad(m_CurrentOperation);
    ui.actionUndo->setEnabled(true);
    if(m_CurrentOperation>=m_LastOperation) ui.actionRedo->setEnabled(false);
  }
  else {
    resetOperationCounter();
    QMessageBox::critical(this, "de-activated", "Redo is not doing anything at the moment!");
  }
}

void GuiMainWindow::resetOperationCounter()
{
  m_CurrentOperation=-1;
  m_LastOperation=m_CurrentOperation;
  ui.actionUndo->setEnabled(false);
  ui.actionRedo->setEnabled(false);
}

QString GuiMainWindow::getXmlSection(QString name)
{
  return m_XmlHandler->getXmlSection(name);
/*  QStringList tags = name.toLower().split("/", QString::SkipEmptyParts);
  QDomElement element = m_XmlDoc.documentElement();
  bool found = true;
  QString section_text = "";
  try {
    foreach (QString tag, tags) {
      QDomNodeList nodes = element.elementsByTagName(tag);
      if (nodes.size() > 1) {
        EG_ERR_RETURN("error retrieving XML section '" + name + "'");
      }
      if (nodes.size() == 0) {
        found = false;
        break;
      }
      if (!nodes.at(0).isElement()) {
        EG_ERR_RETURN("error retrieving XML section '" + name + "'");
      }
      element = nodes.at(0).toElement();
    }
  } catch (Error err) {
    err.display();
  }
  if (found) {
    section_text = element.text();
  }
  return section_text;*/
}

void GuiMainWindow::setXmlSection(QString name, QString contents)
{
  m_XmlHandler->setXmlSection(name,contents);
/*  QStringList tags = name.toLower().split("/", QString::SkipEmptyParts);
  QDomElement element = m_XmlDoc.documentElement();
  try {
    foreach (QString tag, tags) {
      QDomNodeList nodes = element.elementsByTagName(tag);
      if (nodes.size() > 1) {
        EG_ERR_RETURN("error retrieving XML section '" + name + "'");
      }
      if (nodes.size() == 0) {
        QDomElement new_element = m_XmlDoc.createElement(tag);
        element.appendChild(new_element);
        element = new_element;
      } else if (!nodes.at(0).isElement()) {
        EG_ERR_RETURN("error retrieving XML section '" + name + "'");
      } else {
        element = nodes.at(0).toElement();
      }
    }
    while (element.hasChildNodes()) {
      element.removeChild(element.firstChild());
    }
    QDomText text_node = m_XmlDoc.createTextNode(contents);
    element.appendChild(text_node);
  } catch (Error err) {
    err.display();
  }*/
}

void GuiMainWindow::openPhysicalBoundaryConditions()
{
  m_PhysicalBoundaryConditionsMap.clear();
  QString buffer = getXmlSection("engrid/physical");
  QTextStream f(&buffer, QIODevice::ReadOnly);
  while (!f.atEnd()) {
    QString name, type;
    int index;
    f >> index >> name >> type;
    if ((name != "") && (type != "")) {
      PhysicalBoundaryCondition PBC;
      PBC.setName(name);
      PBC.setIndex(index);
      PBC.setType(type);
      for (int i = 0; i < PBC.getNumVars(); ++i) {
        double v;
        f >> v;
        PBC.setValue(i, v);
      }
      m_PhysicalBoundaryConditionsMap[name] = PBC;
    }
  }
}

void GuiMainWindow::savePhysicalBoundaryConditions()
{
  QString buffer("");
  QTextStream f(&buffer, QIODevice::WriteOnly);
  f << "\n";
  foreach (PhysicalBoundaryCondition PBC, m_PhysicalBoundaryConditionsMap) {
    f << PBC.getIndex() << " " << PBC.getName() << " " << PBC.getType();
    for (int i = 0; i < PBC.getNumVars(); ++i) {
      f << " " << PBC.getVarValue(i);
    }
    f << "\n";
  }
  setXmlSection("engrid/physical", buffer);
}

void GuiMainWindow::openBC()
{
  m_bcmap.clear();
  m_VolMap.clear();
  QString buffer = getXmlSection("engrid/bc");
  QTextStream f(&buffer, QIODevice::ReadOnly);
  while (!f.atEnd()) {
    QString name, type;
    int i;
    f >> i >> name >> type;
    if(name!="" && type!="") {
      if (i >= 0) {
        m_bcmap[i] = BoundaryCondition(name,type);
      } else {
        VolumeDefinition V(name, -i);
        QString text = type.replace(",", " ").replace(":", " ");
        QTextStream s(&text);
        while (!s.atEnd()) {
          QString bc_txt, sign_txt;
          s >> bc_txt >> sign_txt;
          V.addBC(bc_txt.toInt(), sign_txt.toInt());
        }
        m_VolMap[name] = V;
      }
    }
  }
}

void GuiMainWindow::saveBC()
{
  QString buffer("");
  QTextStream f(&buffer, QIODevice::WriteOnly);
  f << "\n";
  foreach (int i, m_AllBoundaryCodes) {
    BoundaryCondition bc = m_bcmap[i];
    f << i << " " << bc.getName() << " " << bc.getType() << "\n";
  }
  foreach (VolumeDefinition V, m_VolMap) {
    QString dirs = "";
    bool first = true;
    foreach (int i, m_AllBoundaryCodes) {
      BoundaryCondition bc = m_bcmap[i];
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
  setXmlSection("engrid/bc", buffer);
}

void GuiMainWindow::openGrid(QString file_name)
{
  file_name += ".vtu";
  EG_VTKSP(vtkXMLUnstructuredGridReader,vtu);
  vtu->SetFileName(qPrintable(file_name));
  vtu->Update();
  m_Grid->DeepCopy(vtu->GetOutput());
  createBasicFields(m_Grid, m_Grid->GetNumberOfCells(), m_Grid->GetNumberOfPoints());
  openBC();
  openPhysicalBoundaryConditions();
  updateBoundaryCodes(true);
  createIndices(m_Grid);
  updateActors();
  updateStatusBar();
  zoomAll();
}

///\todo I think this should also be a done by a subclass of IOOperation just like for import operations
void GuiMainWindow::open()
{
  QFileDialog dialog(NULL, "open grid from file", getCwd(), "enGrid case files (*.egc *.EGC);; legacy grid files(*.vtu *.VTU)");
  QFileInfo file_info(m_CurrentFilename);
//   qDebug()<<"m_CurrentFilename="<<m_CurrentFilename;
  dialog.selectFile(file_info.completeBaseName() + ".egc");
  if (dialog.exec()) {
    QStringList selected_files = dialog.selectedFiles();
    QString file_name = selected_files[0];
    if (!file_name.isNull()) {
      this->open(file_name);
    }
  }
}

void GuiMainWindow::open(QString file_name, bool update_current_filename)
{
  cout << "Opening " << qPrintable(file_name) << endl;
  
  QFileInfo file_info(file_name);
  bool no_case_file = false;
  QString file_extension = getExtension(file_name);
  QString grid_file_name = file_name;
  if (file_extension.toLower() == "vtu") {
    no_case_file = true;
    grid_file_name = stripFromExtension(file_name);
  }
  if (!no_case_file) {
    if(!openXml(file_name)) {
      QMessageBox::critical(this, tr("Open failed"), tr("Error reading enGrid case file:\n%1").arg(file_name));
      return;
    }
  }
  if(update_current_filename) {
    GuiMainWindow::setCwd(QFileInfo(file_name).absolutePath());
  }
  openGrid(grid_file_name);
  openBC();
  openPhysicalBoundaryConditions();
  // update current filename
  if(update_current_filename) m_CurrentFilename = stripFromExtension(file_name) + ".egc";
  setWindowTitle(m_CurrentFilename + " - enGrid - " + QString("%1").arg(m_CurrentOperation) );
  setUnsaved(false);

  if(update_current_filename) {
    this->addRecentFile(file_name,QDateTime::currentDateTime());
//     qDebug()<<"Setting new latest file to "<<file_name;
    m_qset.setValue("LatestFile",file_name);
    resetOperationCounter();
    quickSave();
  }
}

bool GuiMainWindow::openXml(QString file_name)
{
  return m_XmlHandler->openXml(file_name);
  /*  QFile xml_file(file_name);
  if (!xml_file.open(QIODevice::ReadOnly)) {
    QString error_message = "Failed to open xml_file " + xml_file.fileName();
    error_message += QString("\n") + tr("Error reading enGrid case file:\n%1").arg(file_name);
    error_message += QString("\n") + "QDir::current()=" + QDir::current().absolutePath();
    error_message += QString("\n") + "QDir::currentPath()=" + QDir::currentPath();
    error_message += QString("\n") + "getCwd()=" + getCwd();
    
    qWarning()<<error_message;
    QMessageBox::critical(this, tr("Open failed"), error_message);
    return(false);
  }
  if (!m_XmlDoc.setContent(&xml_file)) {
    QMessageBox::critical(this, tr("Open failed"), tr("Error reading enGrid case file:\n%1").arg(file_name));
    return(false);
  }
  xml_file.close();*/
}

bool GuiMainWindow::saveXml(QString file_name)
{
  return m_XmlHandler->saveXml(file_name);
/*  QString buffer = m_XmlDoc.toString(0);
  QFile xml_file(file_name);
  xml_file.open(QIODevice::WriteOnly | QIODevice::Text);
  QTextStream f(&xml_file);
  f << buffer << endl;
  return(true);*/
}

QString GuiMainWindow::saveAs(QString file_name, bool update_current_filename)
{
  QString buffer = m_XmlDoc.toString(0);
  
  if(update_current_filename) QApplication::setOverrideCursor(QCursor(Qt::WaitCursor));
  
  QFileInfo file_info(file_name);
  if (file_info.suffix().toLower() != "egc") {
    file_name += ".egc";
  }
  cout << "Saving as " << qPrintable(file_name) << endl;
  if(update_current_filename) {
    // update current filename
    GuiMainWindow::setCwd(file_info.absolutePath());
    m_CurrentFilename = file_name;
  }
  if(!saveGrid(m_Grid, file_name)) {
    QMessageBox::critical(this, QObject::tr("Save failed"), QObject::tr("The grid could not be saved as:\n%1").arg(file_name));
  }
  
  saveBC();
  savePhysicalBoundaryConditions();
  saveXml(file_name);
  
  setWindowTitle(m_CurrentFilename + " - enGrid - " + QString("%1").arg(m_CurrentOperation) );
  setUnsaved(false);
  
  if(update_current_filename) QApplication::restoreOverrideCursor();
  
  if(update_current_filename) {
    this->addRecentFile(file_name,QDateTime::currentDateTime());
//     qDebug()<<"Setting new latest file to "<<file_name;
    m_qset.setValue("LatestFile",file_name);
  }
  
  return(file_name);
}

void GuiMainWindow::save()
{
  if ( m_CurrentFilename == "untitled.egc" || m_UnSaved ) {
    saveAs();
  } else {
    saveAs(m_CurrentFilename);
  }
}

void GuiMainWindow::saveAs()
{
  QFileDialog dialog(NULL, "write case to file", getCwd(), "enGrid case files (*.egc)");
  QFileInfo file_info(m_CurrentFilename);
  dialog.selectFile(file_info.completeBaseName() + ".egc");
  dialog.setAcceptMode(QFileDialog::AcceptSave);
  dialog.setConfirmOverwrite(true);
  if (dialog.exec()) {
    QStringList selected_files = dialog.selectedFiles();
    QString file_name = selected_files[0];
    if (!file_name.isNull()) {
      saveAs(file_name);
      //for the undo/redo operations
      resetOperationCounter();
      quickSave();
    }
  }
}

void GuiMainWindow::updateStatusBar()
{
  QString num, txt = "enGrid is currently busy with an operation ...";
  if (!m_Busy) {
    txt = "";
  }
  if (!tryLock()) {
    m_StatusLabel->setText(txt);
    ui.label_node_cell_info->setText(txt);
    return;
  }
  vtkIdType Ncells = m_Grid->GetNumberOfCells();
  vtkIdType Nnodes = m_Grid->GetNumberOfPoints();
  vtkIdType Ntris  = 0;
  vtkIdType Nquads = 0;
  vtkIdType Ntets  = 0;
  vtkIdType Npyras = 0;
  vtkIdType Nprism = 0;
  vtkIdType Nhexas = 0;
  for (vtkIdType i = 0; i < Ncells; ++i) {
    int ct = m_Grid->GetCellType(i);
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
    vtkIdType id_cell = m_PickedCell;
    if (id_cell < 0 || id_cell>=m_Grid->GetNumberOfCells()) {
      pick_txt += "no cell picked";
    } else {
      vtkIdType type_cell = m_Grid->GetCellType(id_cell);
      if      (type_cell == VTK_TRIANGLE)   pick_txt += "tri";
      else if (type_cell == VTK_QUAD)       pick_txt += "qua";
      else if (type_cell == VTK_TETRA)      pick_txt += "tet";
      else if (type_cell == VTK_PYRAMID)    pick_txt += "pyr";
      else if (type_cell == VTK_WEDGE)      pick_txt += "pri";
      else if (type_cell == VTK_HEXAHEDRON) pick_txt += "hex";
      vtkIdType N_pts, *pts;
      m_Grid->GetCellPoints(id_cell, N_pts, pts);
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
      EG_VTKDCC(vtkIntArray, cell_code, m_Grid, "cell_code");
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
    vtkIdType id_node = m_PickedPoint;
    if (id_node < 0) {
      pick_txt += "no node picked";
    } else {
      vec3_t x;
      m_Grid->GetPoints()->GetPoint(id_node,x.data());
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
      EG_VTKDCN(vtkDoubleArray, characteristic_length_desired, m_Grid, "node_meshdensity_desired");
      tmp.setNum(characteristic_length_desired->GetValue(id_node));
      pick_txt += " wanted density=" + tmp;
      EG_VTKDCN(vtkDoubleArray, node_meshdensity_current, m_Grid, "node_meshdensity_current");
      tmp.setNum(node_meshdensity_current->GetValue(id_node));
      pick_txt += " current density=" + tmp;
      EG_VTKDCN(vtkIntArray, node_specified_density, m_Grid, "node_specified_density");
      tmp.setNum(node_specified_density->GetValue(id_node));
      pick_txt += " node_specified_density=" + tmp;
      EG_VTKDCN(vtkCharArray, node_type, m_Grid, "node_type");
      pick_txt += " type=" + QString(VertexType2Str( node_type->GetValue(id_node)));
      tmp.setNum(id_node);
      pick_txt += " id_node=" + tmp;
    }
    
    txt += pick_txt;
  }
  
  m_StatusLabel->setText(txt);
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
    EG_VTKDCC(vtkIntArray, cell_code, m_Grid, "cell_code");
    for (vtkIdType i = 0; i < m_Grid->GetNumberOfCells(); ++i) {
      int ct = m_Grid->GetCellType(i);
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
    m_Axes->VisibilityOn();
  } else {
    m_Axes->VisibilityOff();
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
  int N = m_Grid->GetNumberOfPoints();
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
      m_NodeTextVectorText[i]->SetText(qPrintable(tmp));
      m_NodeTextPolyDataMapper[i]=vtkPolyDataMapper::New();
      m_NodeTextPolyDataMapper[i]->SetInputConnection(m_NodeTextVectorText[i]->GetOutputPort());
      m_NodeTextFollower[i]=vtkFollower::New();
      m_NodeTextFollower[i]->SetMapper(m_NodeTextPolyDataMapper[i]);
      m_NodeTextFollower[i]->SetScale(m_ReferenceSize, m_ReferenceSize, m_ReferenceSize);
      vec3_t M;
      m_Grid->GetPoint(i, M.data());
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
  vtkIdType N = m_Grid->GetNumberOfCells();
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
        EG_VTKDCC(vtkIntArray, current_cell_field, m_Grid, qPrintable(ui.comboBox_CellTextField->currentText()));
        tmp.setNum(current_cell_field->GetValue(id_cell));
      }
      else EG_BUG;
      
      m_CellTextVectorText[id_cell]->SetText(qPrintable(tmp));
      m_CellTextPolyDataMapper[id_cell]=vtkPolyDataMapper::New();
      m_CellTextPolyDataMapper[id_cell]->SetInputConnection(m_CellTextVectorText[id_cell]->GetOutputPort());
      m_CellTextFollower[id_cell]=vtkFollower::New();
      m_CellTextFollower[id_cell]->SetMapper(m_CellTextPolyDataMapper[id_cell]);
      m_CellTextFollower[id_cell]->SetScale(m_ReferenceSize, m_ReferenceSize, m_ReferenceSize);
      vtkIdType N_pts,*pts;
      m_Grid->GetCellPoints(id_cell,N_pts,pts);
      vec3_t Center(0,0,0);
      for (int p = 0; p < N_pts; ++p) {
        vec3_t M;
        m_Grid->GetPoint(pts[p],M.data());
        Center+=M.data();
      }
      vec3_t OffSet = m_ReferenceSize*triNormal(m_Grid, pts[0], pts[1], pts[2]).normalise();
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

void GuiMainWindow::pickCallBack
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
}

vtkIdType GuiMainWindow::getPickedCell()
{
  if(!ui.radioButton_CellPicker->isChecked()) return(-1);
  
  vtkIdType picked_cell = -1;
  if (m_Grid->GetNumberOfCells() > 0) {
    m_BCodesFilter->Update();
    if (m_BCodesFilter->GetOutput()->GetNumberOfCells() > 0) {
      EG_VTKDCC(vtkLongArray_t, cell_index, m_BCodesFilter->GetOutput(), "cell_index");
      if (m_UseVTKInteractor) {
        picked_cell = cell_index->GetValue(m_CellPicker->GetCellId());
      } else {
        picked_cell = m_PickedCell;
      }
    }
  }
  return picked_cell;
}

vtkIdType GuiMainWindow::getPickedPoint()
{
  if(ui.radioButton_CellPicker->isChecked()) return(-1);
  
  vtkIdType picked_point = -1;
  if (m_Grid->GetNumberOfCells() > 0) {
    m_BCodesFilter->Update();
    if (m_BCodesFilter->GetOutput()->GetNumberOfCells() > 0) {
      EG_VTKDCN(vtkLongArray_t, node_index, m_BCodesFilter->GetOutput(), "node_index");
      if (m_UseVTKInteractor) {
        picked_point = node_index->GetValue(m_PointPicker->GetPointId());
      } else {
        picked_point = m_PickedPoint;
      }
    }
  }
  return picked_point;
}

void GuiMainWindow::changeSurfaceOrientation()
{
  for (vtkIdType cellId = 0; cellId < m_Grid->GetNumberOfCells(); ++cellId) {
    vtkIdType Npts, *pts;
    m_Grid->GetCellPoints(cellId, Npts, pts);
    QVector<vtkIdType> nodes(Npts);
    for (vtkIdType j = 0; j < Npts; ++j) nodes[j]          = pts[j];
    for (vtkIdType j = 0; j < Npts; ++j) pts[Npts - j - 1] = nodes[j];
  }
  updateActors();
  m_Grid->Modified();// to make sure VTK notices the changes and changes the cell colors
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
  stl.setFileTypeToASCII();
  stl();
}

void GuiMainWindow::exportBinaryStl()
{
  StlWriter stl;
  stl.setFileTypeToBinary();
  stl();
}

void GuiMainWindow::exportAsciiPly()
{
  PlyWriter ply;
  ply.setFileTypeToASCII();
  ply();
}

void GuiMainWindow::exportBinaryPly()
{
  PlyWriter ply;
  ply.setFileTypeToBinary();
  ply();
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
  fix->setLockGui();
  (*fix)();
  updateBoundaryCodes(false);
  updateActors();
}

void GuiMainWindow::callDeletePickedPoint()
{
  EG_STDINTERSLOT( DeletePickedPoint );
}

void GuiMainWindow::editBoundaryConditions()
{
  GuiEditBoundaryConditions editbcs;
  editbcs.setBoundaryCodes(m_AllBoundaryCodes);
  editbcs.setMap(&m_bcmap);
  editbcs();
}

void GuiMainWindow::configure()
{
  {
    // Just to create initial entries in the settings file 
    // so that the options menu isn't empty at first start.
    try {
      GridSmoother tmp01;
      GuiCreateBoundaryLayer tmp02;
      SurfaceProjection tmp03;
      SurfaceMesher tmp04;
      UpdateDesiredMeshDensity tmp05;
      InsertPoints tmp06;
      RemovePoints tmp07;
      LaplaceSmoother tmp08;
      SwapTriangles tmp09;
      OpenFOAMTools tmp10;
    } catch (Error err) {
      err.display();
    }
  }
  GuiSettingsViewer settings(&m_qset);
  settings.CreateViewer();
  settings.exec();
  
  getSet("General","enable undo+redo",false,m_undo_redo_enabled);
}

void GuiMainWindow::about()
{
  QMessageBox box(this);
  
  QString title="ENGRID";
  QString version = QString("version ") + ENGRID_VERSION;
  #ifdef GIT_DESCRIBE
  if(!QString(GIT_DESCRIBE).isEmpty()) {
      version +=  QString(" - ") + GIT_DESCRIBE;
  }
  #endif
  
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

///\todo Why not use bcs = m_AllBoundaryCodes; ?
void GuiMainWindow::getAllBoundaryCodes(QSet<int> &bcs)
{
  bcs = m_AllBoundaryCodes;
//   qWarning()<<"m_AllBoundaryCodes="<<m_AllBoundaryCodes;
//   bcs.clear();
//   foreach (int bc, m_AllBoundaryCodes) {
//     bcs.insert(bc);
//   }
}

QSet<int> GuiMainWindow::getAllBoundaryCodes()
{
  return m_AllBoundaryCodes;
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
  foreach(VolumeDefinition vol, m_VolMap) {
    vols.push_back(vol);
  }
  return vols;
}

void GuiMainWindow::setAllVols(QList<VolumeDefinition> vols)
{
  m_VolMap.clear();
  foreach (VolumeDefinition V, vols) {
    m_VolMap[V.getName()] = V;
  }
}

QList<PhysicalBoundaryCondition> GuiMainWindow::getAllPhysicalBoundaryConditions()
{
  QList<PhysicalBoundaryCondition> physical_boundary_conditions;
  foreach(PhysicalBoundaryCondition PBC, m_PhysicalBoundaryConditionsMap) {
    physical_boundary_conditions.push_back(PBC);
  }
  return physical_boundary_conditions;
}

void GuiMainWindow::setAllPhysicalBoundaryConditions(QList<PhysicalBoundaryCondition> physical_boundary_conditions)
{
  m_PhysicalBoundaryConditionsMap.clear();
  foreach (PhysicalBoundaryCondition PBC, physical_boundary_conditions) {
    m_PhysicalBoundaryConditionsMap[PBC.getName()] = PBC;
  }
}

void GuiMainWindow::setAllPhysicalBoundaryConditions(QMap<QString,PhysicalBoundaryCondition> physical_boundary_conditions) {
  m_PhysicalBoundaryConditionsMap = physical_boundary_conditions;
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
  QFileInfo fileinfo(m_CurrentFilename);
  return fileinfo.absolutePath()+"/";
}

void GuiMainWindow::markOutputLine()
{
  cout << "\n****************************************\n";
  cout << qPrintable(QTime::currentTime().toString("hh:mm:ss"));
  cout << "\n****************************************\n" << endl;
}

void GuiMainWindow::storeSurfaceProjection()
{
  qDebug()<<"@@@ GuiMainWindow::storeSurfaceProjection called";
  foreach (SurfaceProjection* proj, m_SurfProj) {
    delete proj;
  }
  m_SurfProj.clear();
  cout << "storing background grid for surface projection:" << endl;
//   EG_VTKSP(vtkUnstructuredGrid,new_grid);
  MeshPartition new_grid_partition;
  bool first = true;
  
  QFileInfo file_info(m_CurrentFilename);
  
  foreach (int bc, m_AllBoundaryCodes) {
    SurfaceProjection *proj = new SurfaceProjection();
    m_SurfProj[bc] = proj;
    QSet<int> bcs;
    bcs.insert(bc);
    QVector<vtkIdType> cls;
    getSurfaceCells(bcs, cls, m_Grid);
    proj->setBackgroundGrid(m_Grid, cls);
    if (proj->usesLevelSet()) {
      QString file_name;
      file_name.setNum(bc);
      file_name = "OctreeBC" + file_name;
      proj->writeOctree(file_name);
      cout << "  bc " << bc << ": " << proj->getNumOctreeCells() << endl;
    }
    QString basename = file_info.completeBaseName() + "_" + QString::number(bc);
    
    proj->m_ExactMode = 0;
    
//     DebugLevel = 100;
    if(DebugLevel>100) {
      proj->writeGridWithNormals(basename);
      proj->writeInterpolationGrid(basename);
      proj->writeTriangleGrid(basename);
      qDebug()<<"=====> bc="<<bc<<" proj->getBezierGrid()->GetNumberOfPoints()="<<proj->getBezierGrid()->GetNumberOfPoints()
        <<" proj->getBezierGrid()->GetNumberOfCells()="<<proj->getBezierGrid()->GetNumberOfCells();
      
      if(first) {
        first = false;
        new_grid_partition.setGrid(proj->getBezierGrid());
        new_grid_partition.setAllCells();
      }
      else {
        MeshPartition grid_partition(proj->getBezierGrid(), true);
        new_grid_partition.addPartition(grid_partition);
      }
    }
  }
  
  if(DebugLevel>100) writeGrid(new_grid_partition.getGrid(), file_info.completeBaseName() + "_projection_surface");
//   DebugLevel = 0;
}

SurfaceProjection* GuiMainWindow::getSurfProj(int bc)
{
  if (!m_SurfProj.contains(bc)) {
    QString bc_txt;
    bc_txt.setNum(bc);
    EG_ERR_RETURN("No surface projection found for boundary code " + bc_txt);
  }
  return m_SurfProj[bc];
}

bool GuiMainWindow::checkSurfProj()
{
  bool ok = true;
  foreach (int bc, m_AllBoundaryCodes) {
    if (!m_SurfProj.contains(bc)) {
      ok = false;
      break;
    }
  }
  return ok;
}

void GuiMainWindow::openRecent(QAction *action)
{
  qDebug()<<"GuiMainWindow::openRecent called";
  QString file_name = action->text().right(action->text().length()-23);
  this->open(file_name);
}

void GuiMainWindow::readRecentFiles()
{
  m_RecentFiles.clear();
  this->recentFileMenu()->clear();
  QStringList file_names = m_qset.value("FileNames").toStringList();
  QStringList file_dates = m_qset.value("FileDates").toStringList();
  int N = min(10,m_qset.value("NumberOfFiles").toInt());
//   cout << "NumberOfFiles=" << N << endl;
  for (int i = 0; i < N; ++i) {
    QString new_file = file_names.at(i);
    QString date_text = file_dates.at(i);
    QDateTime date = QDateTime::fromString(date_text,"dd.MM.yyyy_hh:mm:ss");
    addRecentFile(new_file,date);
  }
}

void GuiMainWindow::writeRecentFiles()
{
  m_qset.setValue("NumberOfFiles",m_RecentFiles.size());
  QStringList file_names;
  QStringList file_dates;
  for (QMap<QString,QDateTime>::iterator i = m_RecentFiles.begin(); i != m_RecentFiles.end(); ++i) {
    QString file_name = i.key();
    QString date_text = i.value().toString("dd.MM.yyyy_hh:mm:ss");
    file_names.append(file_name);
    file_dates.append(date_text);
  }
  m_qset.setValue("FileNames",file_names);
  m_qset.setValue("FileDates",file_dates);
}

void GuiMainWindow::addRecentFile(QString file_name, QDateTime date)
{
  m_RecentFiles[file_name] = date;
  while (m_RecentFiles.size() > 10) {
    QMap<QString,QDateTime>::iterator i,j;
    QDateTime old = QDateTime::currentDateTime();
    for (i = m_RecentFiles.begin(); i != m_RecentFiles.end(); ++i) {
      if (i.value() <= old) {
        old = i.value();
        j = i;
      }
    }
    m_RecentFiles.erase(j);
  }
  this->recentFileMenu()->clear();
  QMap<int,QString> menu_map;
  QDateTime now = QDateTime::currentDateTime();
  for (QMap<QString,QDateTime>::iterator i = m_RecentFiles.begin(); i != m_RecentFiles.end(); ++i) {
    QString action_text = i.value().toString("dd.MM.yyyy hh:mm:ss");
    action_text += " -> ";
    action_text += i.key();
    menu_map[i.value().secsTo(now)] = action_text;
  }
  {
    for (QMap<int,QString>::iterator i = menu_map.begin(); i != menu_map.end(); ++i) {
      QAction *action = new QAction(i.value(),this);
      this->recentFileMenu()->addAction(action);
    }
  }
}

void GuiMainWindow::callInsertNewCell()
{
  bool ok1,ok2,ok3,ok4;
  vtkIdType pts[3];
  pts[0] = QInputDialog::getInt(this, tr("id_node1"),tr("id_node1:"), 0, 0, m_Grid->GetNumberOfPoints(), 1, &ok1);
  pts[1] = QInputDialog::getInt(this, tr("id_node2"),tr("id_node2:"), 0, 0, m_Grid->GetNumberOfPoints(), 1, &ok2);
  pts[2] = QInputDialog::getInt(this, tr("id_node3"),tr("id_node3:"), 0, 0, m_Grid->GetNumberOfPoints(), 1, &ok3);
  vtkIdType id_cell = QInputDialog::getInt(this, tr("copy cell data from id_cell"),tr("copy cell data from id_cell:"), 0, 0, m_Grid->GetNumberOfCells(), 1, &ok4);
  if (ok1 && ok2 && ok3 && ok4) {
    EG_VTKSP( vtkUnstructuredGrid, new_grid );
    allocateGrid( new_grid, m_Grid->GetNumberOfCells() + 1, m_Grid->GetNumberOfPoints() );
    makeCopyNoAlloc(m_Grid, new_grid);
    vtkIdType id_new_cell = new_grid->InsertNextCell(VTK_TRIANGLE, 3, pts);
    copyCellData(m_Grid, id_cell, new_grid, id_new_cell);
    makeCopy(new_grid, m_Grid);
    m_Grid->Modified();
    QMessageBox::information(NULL, "new cell", tr("The new cell has ID = %1").arg(id_new_cell));
    qDebug()<<tr("The new cell has ID = %1").arg(id_new_cell);
  }
}

void GuiMainWindow::callMergeNodes()
{
  bool ok1,ok2;
  vtkIdType id_node1 = QInputDialog::getInt(this, tr("id_node1"),tr("id_node1:"), 0, 0, m_Grid->GetNumberOfPoints(), 1, &ok1);
  vtkIdType id_node2 = QInputDialog::getInt(this, tr("id_node2"),tr("id_node2:"), 0, 0, m_Grid->GetNumberOfPoints(), 1, &ok2);
  if (ok1 && ok2) {
    EG_VTKSP( vtkUnstructuredGrid, new_grid );
    allocateGrid( new_grid, m_Grid->GetNumberOfCells(), m_Grid->GetNumberOfPoints() - 1 );
    
    QVector<vtkIdType> old2new_nodes(m_Grid->GetNumberOfPoints(), -1);
    QVector<vtkIdType> old2new_cells(m_Grid->GetNumberOfCells(), -1);
    
    vtkIdType id_new_node = 0;
    for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
      if(id_node!=id_node1 && id_node!=id_node2) {
        vec3_t x;
        m_Grid->GetPoints()->GetPoint(id_node, x.data());
        new_grid->GetPoints()->SetPoint(id_new_node, x.data());
        copyNodeData(m_Grid, id_node, new_grid, id_new_node);
        old2new_nodes[id_node] = id_new_node;
        id_new_node++;
      }
      else if(id_node==id_node1) {
        vec3_t x1;
        m_Grid->GetPoints()->GetPoint(id_node1, x1.data());
        vec3_t x2;
        m_Grid->GetPoints()->GetPoint(id_node2, x2.data());
        vec3_t x = 0.5*(x1+x2);
        new_grid->GetPoints()->SetPoint(id_new_node, x.data());
        copyNodeData(m_Grid, id_node, new_grid, id_new_node);
        old2new_nodes[id_node1] = id_new_node;
        old2new_nodes[id_node2] = id_new_node;
        id_new_node++;
      }
      else {
      }
    }
    
    for (vtkIdType id_cell = 0; id_cell < m_Grid->GetNumberOfCells(); ++id_cell) {
      vtkIdType N_pts, *pts;
      vtkIdType type_cell = m_Grid->GetCellType(id_cell);
      m_Grid->GetCellPoints(id_cell, N_pts, pts);
      QVector<vtkIdType> new_pts(N_pts);
      for (int i = 0; i < N_pts; ++i) {
        new_pts[i] = old2new_nodes[pts[i]];
      }
      vtkIdType id_new_cell = new_grid->InsertNextCell(type_cell, N_pts, new_pts.data());
      copyCellData(m_Grid, id_cell, new_grid, id_new_cell);
    }
    
    makeCopy(new_grid, m_Grid);
    m_Grid->Modified();
    qDebug()<<"The fusion is complete.";
  }
  
}
