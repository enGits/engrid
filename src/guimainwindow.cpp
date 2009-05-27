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

// GuiOutputWindow::GuiOutputWindow()
// {
//   ui.setupUi(this);
// };

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
  
  QPoint pos = qset.value("pos", QPoint(200, 200)).toPoint();
  QSize size = qset.value("size", QSize(400, 400)).toSize();
  resize(size);
  move(pos);
  
//   dock_widget = new QDockWidget(tr("output"), this);
//   dock_widget->setFeatures(QDockWidget::AllDockWidgetFeatures);
//   output_window = new GuiOutputWindow();
//   dock_widget->setWidget(output_window);
//   addDockWidget(Qt::LeftDockWidgetArea, dock_widget);
//   ui.menuView->addAction(dock_widget->toggleViewAction());
  
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
  connect(ui.actionZoomOnPickedObject,     SIGNAL(activated()),       this, SLOT(ZoomOnPickedObject()));
  connect(ui.actionPrintGrid,              SIGNAL(activated()),       this, SLOT(PrintGrid()));
  connect(ui.actionShowInfo,               SIGNAL(activated()),       this, SLOT(Info()));
  connect(ui.actionDeselectAll,            SIGNAL(activated()),       this, SLOT(DeselectAll()));
  connect(ui.actionOpen,                   SIGNAL(activated()),       this, SLOT(open()));
  connect(ui.actionSave,                   SIGNAL(activated()),       this, SLOT(save()));
  connect(ui.actionSaveAs,                 SIGNAL(activated()),       this, SLOT(saveAs()));
  connect(ui.actionBoundaryCodes,          SIGNAL(activated()),       this, SLOT(selectBoundaryCodes()));
  connect(ui.actionNormalExtrusion,        SIGNAL(activated()),       this, SLOT(normalExtrusion()));
  connect(ui.actionViewAxes,               SIGNAL(changed()),         this, SLOT(setAxesVisibility()));
  connect(ui.actionViewOrthogonal,         SIGNAL(changed()),         this, SLOT(setViewingMode()));
  connect(ui.actionViewNodeIDs,            SIGNAL(changed()),         this, SLOT(ViewNodeIDs()));
  connect(ui.actionViewCellIDs,            SIGNAL(changed()),         this, SLOT(ViewCellIDs()));
  connect(ui.actionChangeOrientation,      SIGNAL(activated()),       this, SLOT(changeSurfaceOrientation()));
  connect(ui.actionCheckOrientation,       SIGNAL(activated()),       this, SLOT(checkSurfaceOrientation()));
  connect(ui.actionImproveAspectRatio,     SIGNAL(activated()),       this, SLOT(improveAspectRatio()));
  connect(ui.actionRedraw,                 SIGNAL(activated()),       this, SLOT(updateActors()));
  connect(ui.actionScaleToData,            SIGNAL(activated()),       this, SLOT(ScaleToData()));
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
  
  
#include "std_connections.h"
  
  if (qset.contains("working_directory")) {
    cwd = qset.value("working_directory").toString();
  };
  grid = vtkUnstructuredGrid::New();
  renderer = vtkRenderer::New();
  getRenderWindow()->AddRenderer(renderer);
  surface_actor = vtkActor::New();
  surface_wire_actor = vtkActor::New();
  
  tetra_mapper        = vtkPolyDataMapper::New();
  pyramid_mapper      = vtkPolyDataMapper::New();
  wedge_mapper        = vtkPolyDataMapper::New();
  hexa_mapper         = vtkPolyDataMapper::New();
  volume_wire_mapper  = vtkPolyDataMapper::New();
  surface_mapper      = vtkPolyDataMapper::New();
  surface_wire_mapper = vtkPolyDataMapper::New();
  
  backface_property = vtkProperty::New();
  
  tetra_actor       = NULL;
  pyramid_actor     = NULL;
  wedge_actor       = NULL;
  hexa_actor        = NULL;
  volume_wire_actor = NULL;
  iamlegend_actor = NULL;
  lut = NULL;
  field_mapper = NULL;
  
  surface_filter = vtkGeometryFilter::New();
  bcodes_filter = vtkEgBoundaryCodesFilter::New();
  renderer->AddActor(surface_actor);
  renderer->AddActor(surface_wire_actor);
  pick_sphere = vtkSphereSource::New();
  pick_mapper = vtkPolyDataMapper::New();
  pick_actor = NULL;
  
  extr_vol        = vtkEgExtractVolumeCells::New();
  extr_tetras     = vtkEgExtractVolumeCells::New();
  extr_pyramids   = vtkEgExtractVolumeCells::New();
  extr_wedges     = vtkEgExtractVolumeCells::New();
  extr_hexas      = vtkEgExtractVolumeCells::New();
  
  volume_geometry  = vtkGeometryFilter::New();
  tetra_geometry   = vtkGeometryFilter::New();
  pyramid_geometry = vtkGeometryFilter::New();
  wedge_geometry   = vtkGeometryFilter::New();
  hexa_geometry   = vtkGeometryFilter::New();
  
  extr_tetras->SetAllOff();
  extr_tetras->SetTetrasOn();
  extr_pyramids->SetAllOff();
  extr_pyramids->SetPyramidsOn();
  extr_wedges->SetAllOff();
  extr_wedges->SetWedgesOn();
  extr_hexas->SetAllOff();
  extr_hexas->SetHexasOn();
  
  boundary_pd = vtkPolyData::New();
  tetras_pd   = vtkPolyData::New();
  wedges_pd   = vtkPolyData::New();
  pyras_pd    = vtkPolyData::New();
  hexas_pd    = vtkPolyData::New();
  volume_pd   = vtkPolyData::New();
  
  current_filename = "untitled.vtu";
  ResetOperationCounter();//clears undo/redo list and disables undo/redo
  setWindowTitle(current_filename + " - enGrid - " + QString("%1").arg(current_operation) );
  
  status_bar = new QStatusBar(this);
  setStatusBar(status_bar);
//   status_label = new QLabel(this);
//   status_label->setWordWrap(true);
//   status_label->setSizePolicy(QSizePolicy::Minimum);
// //   status_label->setVerticalPolicy(QSizePolicy::Minimum);
//   status_bar->addWidget(status_label);

//   QVBoxLayout *status_layout = new QVBoxLayout;
//   status_layout->addWidget(status_label);
//   status_bar->setLayout(status_layout);

  QString txt = "0 volume cells (0 tetras, 0 hexas, 0 pyramids, 0 prisms), ";
  txt += "0 surface cells (0 triangles, 0 quads), 0 nodes";
//   status_label->setText(txt);
  status_bar->showMessage(txt);
  status_bar->setToolTip(txt);
  ui.label_node_cell_info->setText(txt);

  axes = vtkCubeAxesActor2D::New();
  axes->SetCamera(getRenderer()->GetActiveCamera());
  getRenderer()->AddActor(axes);
  setAxesVisibility();
  
  CellPicker = vtkCellPicker::New();
//   getInteractor()->SetPicker(CellPicker);
  PointPicker = vtkPointPicker::New();
//   getInteractor()->SetPicker(PointPicker);
  
  pick_sphere->SetRadius(0.25);//in case the user starts picking points instead of cells

  vtkCallbackCommand *cbc = vtkCallbackCommand::New();
  cbc->SetCallback(pickCellCallBack);
  CellPicker->AddObserver(vtkCommand::EndPickEvent, cbc);
  cbc->SetCallback(pickPointCallBack);
  PointPicker->AddObserver(vtkCommand::EndPickEvent, cbc);
  cbc->Delete();
  
  QString user = QString(getenv("USER"));
  QString basename="enGrid_output.txt";
  
// define temporary path
  QDir dir("/");
  if (qset.contains("tmp_directory")) {
    m_LogDir=qset.value("tmp_directory").toString();
  } else {
    m_LogDir=dir.tempPath();
  };
  m_LogDir = m_LogDir + "/" + "enGrid_"+user + "/";
  dir.mkpath(m_LogDir);
  
  log_file_name = m_LogDir + basename;
  cout<<"log_file_name="<<log_file_name.toLatin1().data()<<endl;

  system_stdout = stdout;
  freopen (log_file_name.toAscii().data(), "w", stdout);
  
  busy = false;
  
//   ui.checkBox_UseVTKInteractor->setCheckState(Qt::Checked);
//   ui.radioButton_CellPicker->setChecked(true);
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
  
  ReferenceSize=0.2;
  
  ui.doubleSpinBox_HueMin->setValue(0.667);
  ui.doubleSpinBox_HueMax->setValue(0);
  
  egvtkInteractorStyle *style = egvtkInteractorStyle::New();
  getInteractor()->SetInteractorStyle(style);
  style->Delete();
};
//end of GuiMainWindow::GuiMainWindow() : QMainWindow(NULL)

GuiMainWindow::~GuiMainWindow()
{
  qset.setValue("pos", pos());
  qset.setValue("size", size());
  
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
  
};

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
    };
    ui.textEditOutput->append(txt);
  };
};

void GuiMainWindow::exit()
{
  QCoreApplication::exit();
};

vtkRenderWindow* GuiMainWindow::getRenderWindow() 
{
  return ui.qvtkWidget->GetRenderWindow();
};

vtkRenderer* GuiMainWindow::getRenderer()
{
  return renderer;
};

QVTKInteractor* GuiMainWindow::getInteractor()
{
  return ui.qvtkWidget->GetInteractor();
};

QString GuiMainWindow::getCwd()
{
  return cwd;
};

void GuiMainWindow::setCwd(QString dir)
{
  cwd = dir;
  qset.setValue("working_directory",dir);
};

void GuiMainWindow::ScaleToData()
{
  int current_field=ui.comboBox_Field->currentIndex();
  if(current_field>0)
  {
    double range[2];
    boundary_pd->GetPointData()->GetArray(current_field-1)->GetRange(range);
    cout<<"current_field="<<current_field<<endl;
    cout<<"range[0]="<<range[0]<<endl;
    cout<<"range[1]="<<range[1]<<endl;
    ui.doubleSpinBox_FieldMin->setRange(range[0],range[1]);
    ui.doubleSpinBox_FieldMax->setRange(range[0],range[1]);
    ui.doubleSpinBox_FieldMin->setValue(range[0]);
    ui.doubleSpinBox_FieldMax->setValue(range[1]);
  }
}

void GuiMainWindow::updateActors()
{
  if (!tryLock()) return;
  try {
    {
      {
        if (!grid->GetCellData()->GetScalars("cell_index")) {
          EG_VTKSP(vtkLongArray_t, cell_idx);
          cell_idx->SetName("cell_index");
          cell_idx->SetNumberOfValues(grid->GetNumberOfCells());
          grid->GetCellData()->AddArray(cell_idx);
        };
      };
      EG_VTKDCC(vtkLongArray_t, cell_index, grid, "cell_index");
      for (vtkIdType cellId = 0; cellId < grid->GetNumberOfCells(); ++cellId) {
        cell_index->SetValue(cellId, cellId);
      };
    };
    
    axes->SetInput(grid);
    
    double xmin =  1e99;
    double xmax = -1e99;
    double ymin =  1e99;
    double ymax = -1e99;
    double zmin =  1e99;
    double zmax = -1e99;
    for (vtkIdType nodeId = 0; nodeId < grid->GetNumberOfPoints(); ++nodeId) {
      vec3_t x;
      grid->GetPoints()->GetPoint(nodeId, x.data());
      xmin = min(x[0], xmin);
      xmax = max(x[0], xmax);
      ymin = min(x[1], ymin);
      ymax = max(x[1], ymax);
      zmin = min(x[2], zmin);
      zmax = max(x[2], zmax);
    };
    
    if (surface_actor) {
      getRenderer()->RemoveActor(surface_actor);
      surface_actor->Delete();
      surface_actor = NULL;
    };
    if (surface_wire_actor) {
      getRenderer()->RemoveActor(surface_wire_actor);
      surface_wire_actor->Delete();
      surface_wire_actor = NULL;
    };
    if (tetra_actor) {
      getRenderer()->RemoveActor(tetra_actor);
      tetra_actor->Delete();
      tetra_actor = NULL;
    };
    if (pyramid_actor) {
      getRenderer()->RemoveActor(pyramid_actor);
      pyramid_actor->Delete();
      pyramid_actor = NULL;
    };
    if (wedge_actor) {
      getRenderer()->RemoveActor(wedge_actor);
      wedge_actor->Delete();
      wedge_actor = NULL;
    };
    if (hexa_actor) {
      getRenderer()->RemoveActor(hexa_actor);
      hexa_actor->Delete();
      hexa_actor = NULL;
    };
    if (volume_wire_actor) {
      getRenderer()->RemoveActor(volume_wire_actor);
      volume_wire_actor->Delete();
      volume_wire_actor = NULL;
    };
    if (pick_actor) {
      getRenderer()->RemoveActor(pick_actor);
      pick_actor->Delete();cout<<"Deleting pick_actor="<<pick_actor<<endl;
      pick_actor = NULL;
    };
    if (iamlegend_actor) {
      getRenderer()->RemoveActor(iamlegend_actor);
      iamlegend_actor->Delete();
      iamlegend_actor = NULL;
    };
    if (lut) {
      lut->Delete();
      lut = NULL;
    };
    if (field_mapper) {
      field_mapper->Delete();
      field_mapper = NULL;
    };
    
    if (ui.checkBoxSurface->isChecked()) {
      bcodes_filter->SetBoundaryCodes(&display_boundary_codes);
      
      bcodes_filter->SetInput(grid);
      
      surface_filter->SetInput(bcodes_filter->GetOutput());
      
      surface_filter->Update();
      
      boundary_pd->DeepCopy(surface_filter->GetOutput());
      
      surface_mapper->SetInput(boundary_pd);
      
      surface_wire_mapper->SetInput(boundary_pd);
      surface_actor = vtkActor::New();
      surface_actor->GetProperty()->SetRepresentationToSurface();
      
      //Fill node field combobox
      int current_field=ui.comboBox_Field->currentIndex();
      ui.comboBox_Field->clear();
      ui.comboBox_Field->addItem("None");
      int N_Arrays=boundary_pd->GetPointData()->GetNumberOfArrays();
//       cout<<"N_Arrays="<<N_Arrays<<endl;
      for (int i=0; i<N_Arrays; i++)
      {
        ui.comboBox_Field->addItem(boundary_pd->GetPointData()->GetArrayName(i));
      }
      if(current_field==-1) ui.comboBox_Field->setCurrentIndex(0);
      else ui.comboBox_Field->setCurrentIndex(current_field);
      
      //Fill cell field combobox
      int current_cell_field=ui.comboBox_CellTextField->currentIndex();
      ui.comboBox_CellTextField->clear();
      ui.comboBox_CellTextField->addItem("Cell ID");
      int N_CellArrays=boundary_pd->GetCellData()->GetNumberOfArrays();
//       cout<<"N_CellArrays="<<N_CellArrays<<endl;
      for (int i=0; i<N_CellArrays; i++)
      {
        ui.comboBox_CellTextField->addItem(grid->GetCellData()->GetArrayName(i));
      }
      if(current_cell_field==-1) ui.comboBox_CellTextField->setCurrentIndex(0);
      else ui.comboBox_CellTextField->setCurrentIndex(current_cell_field);
      
//       cout<<"index="<<ui.comboBox_Field->currentIndex()<<endl;
//       cout<<"name="<<ui.comboBox_Field->currentText().toLatin1().data()<<endl;
      
      current_field=ui.comboBox_Field->currentIndex();
      if(current_field>0)
      {
        double range[2];
        boundary_pd->GetPointData()->GetArray(current_field-1)->GetRange(range);
        cout<<"current_field="<<current_field<<endl;
        cout<<"range[0]="<<range[0]<<endl;
        cout<<"range[1]="<<range[1]<<endl;
        ui.doubleSpinBox_FieldMin->setRange(range[0],range[1]);
        ui.doubleSpinBox_FieldMax->setRange(range[0],range[1]);
      }
      
      if(ui.comboBox_Field->currentIndex()<=0) {
        surface_actor->SetBackfaceProperty(backface_property);
        surface_actor->GetProperty()->SetColor(0.5,1,0.5);
        surface_actor->GetBackfaceProperty()->SetColor(1,1,0.5);
        surface_actor->SetMapper(surface_mapper);
      }
      else {
        lut=vtkLookupTable::New();
        lut->SetNumberOfColors(ui.spinBox_Color->value());
        lut->SetHueRange(ui.doubleSpinBox_HueMin->value(),ui.doubleSpinBox_HueMax->value());
        lut->Build();
        field_mapper=vtkPolyDataMapper::New();
        field_mapper->SetLookupTable(lut);
        field_mapper->SetInput(boundary_pd);
        field_mapper->SetColorModeToMapScalars();
        field_mapper->SetScalarModeToUsePointFieldData();
        field_mapper->ColorByArrayComponent(ui.comboBox_Field->currentText().toLatin1().data(),0);
        field_mapper->SetScalarRange(ui.doubleSpinBox_FieldMin->value(),ui.doubleSpinBox_FieldMax->value());
        surface_actor->SetMapper(field_mapper);
        
        if(ui.checkBox_Legend->checkState()) {
          iamlegend_actor = vtkScalarBarActor::New();
          iamlegend_actor->SetLookupTable (lut);
          getRenderer()->AddActor(iamlegend_actor);
        }
      }
      
/*      vtkDoubleArray *newScalars = vtkDoubleArray::New();
      int index;
      newScalars=(vtkDoubleArray *)boundary_pd->GetPointData()->GetArray("node_meshdensity_current",index);
      cout<<"index="<<index<<endl;*/
      
/*      cout<<"=========="<<endl;
      boundary_pd->GetPointData()->GetArray("node_status",index);
      cout<<"index="<<index<<endl;
      boundary_pd->GetPointData()->GetArray("node_layer",index);
      cout<<"index="<<index<<endl;
      boundary_pd->GetPointData()->GetArray("node_index",index);
      cout<<"index="<<index<<endl;
      boundary_pd->GetPointData()->GetArray("node_meshdensity_desired",index);
      cout<<"index="<<index<<endl;
      boundary_pd->GetPointData()->GetArray("node_meshdensity_current",index);
      cout<<"index="<<index<<endl;
      boundary_pd->GetPointData()->GetArray("node_type",index);
      cout<<"index="<<index<<endl;
      cout<<"=========="<<endl;*/
      
/*      int N2=newScalars->GetNumberOfComponents();
      int N3=newScalars->GetNumberOfTuples();
      cout<<"Number of components=N2="<<N2<<endl;
      cout<<"Number of tuples=N3="<<N3<<endl;*/
      
      
/*      for (int i=0; i<N3; i++)
      {
        double D=newScalars->GetComponent(i,0);//strange, but works. O.o
        cout<<"D["<<i<<"]="<<D<<endl;
      }*/
      
      surface_wire_actor = vtkActor::New();
      surface_wire_actor->GetProperty()->SetRepresentationToWireframe();
      surface_wire_actor->GetProperty()->SetColor(0,0,1);
      
      surface_wire_actor->SetMapper(surface_wire_mapper);
      getRenderer()->AddActor(surface_actor);
      getRenderer()->AddActor(surface_wire_actor);
      bcodes_filter->Update();

      if(ui.checkBox_ShowPickSphere->checkState())
      {
        if(m_UseVTKInteractor)
        {
          
          if(ui.radioButton_CellPicker->isChecked())
          {
            getInteractor()->SetPicker(CellPicker);
            vtkIdType cellId = getPickedCell();
            pickCell(cellId);
          }
          else
          {
            getInteractor()->SetPicker(PointPicker);
            vtkIdType nodeId = getPickedPoint();
            pickPoint(nodeId);
          }
        }
        else
        {
          if(ui.radioButton_CellPicker->isChecked()) pickCell(PickedCell);
          else pickPoint(PickedPoint);
        }
      }
    };
    
    vec3_t x, n;
    x[0] = ui.lineEditClipX->text().toDouble();
    x[1] = ui.lineEditClipY->text().toDouble();
    x[2] = ui.lineEditClipZ->text().toDouble();
    n[0] = ui.lineEditClipNX->text().toDouble();
    n[1] = ui.lineEditClipNY->text().toDouble();
    n[2] = ui.lineEditClipNZ->text().toDouble();
    n.normalise();
    x = x + ui.lineEditOffset->text().toDouble()*n;
    extr_vol->SetAllOff();
    if (ui.checkBoxTetra->isChecked()) {
      extr_vol->SetTetrasOn();
      extr_tetras->SetInput(grid);
      if (ui.checkBoxClip->isChecked()) {
        extr_tetras->SetClippingOn();
        extr_tetras->SetX(x);
        extr_tetras->SetN(n);
      } else {
        extr_tetras->SetClippingOff();
      };
      tetra_actor = vtkActor::New();
      tetra_geometry->SetInput(extr_tetras->GetOutput());
      tetra_geometry->Update();
      tetras_pd->DeepCopy(tetra_geometry->GetOutput());
      tetra_mapper->SetInput(tetras_pd);
      tetra_actor = vtkActor::New();
      tetra_actor->SetMapper(tetra_mapper);
      tetra_actor->GetProperty()->SetColor(1,0,0);
      getRenderer()->AddActor(tetra_actor);
    };
    if (ui.checkBoxPyramid->isChecked()) {
      extr_vol->SetPyramidsOn();
      extr_pyramids->SetInput(grid);
      if (ui.checkBoxClip->isChecked()) {
        extr_pyramids->SetClippingOn();
        extr_pyramids->SetX(x);
        extr_pyramids->SetN(n);
      } else {
        extr_pyramids->SetClippingOff();
      };
      pyramid_actor = vtkActor::New();
      pyramid_geometry->SetInput(extr_pyramids->GetOutput());
      pyramid_geometry->Update();
      pyras_pd->DeepCopy(pyramid_geometry->GetOutput());
      pyramid_mapper->SetInput(pyras_pd);
      pyramid_actor = vtkActor::New();
      pyramid_actor->SetMapper(pyramid_mapper);
      pyramid_actor->GetProperty()->SetColor(1,1,0);
      getRenderer()->AddActor(pyramid_actor);
    };
    if (ui.checkBoxWedge->isChecked()) {
      extr_vol->SetWedgesOn();
      extr_wedges->SetInput(grid);
      if (ui.checkBoxClip->isChecked()) {
        extr_wedges->SetClippingOn();
        extr_wedges->SetX(x);
        extr_wedges->SetN(n);
      } else {
        extr_wedges->SetClippingOff();
      };
      wedge_actor = vtkActor::New();
      wedge_geometry->SetInput(extr_wedges->GetOutput());
      wedge_geometry->Update();
      wedges_pd->DeepCopy(wedge_geometry->GetOutput());
      wedge_mapper->SetInput(wedges_pd);
      wedge_actor = vtkActor::New();
      wedge_actor->SetMapper(wedge_mapper);
      wedge_actor->GetProperty()->SetColor(0,1,0);
      getRenderer()->AddActor(wedge_actor);
    };
    if (ui.checkBoxHexa->isChecked()) {
      extr_vol->SetHexasOn();
      extr_hexas->SetInput(grid);
      if (ui.checkBoxClip->isChecked()) {
        extr_hexas->SetClippingOn();
        extr_hexas->SetX(x);
        extr_hexas->SetN(n);
      } else {
        extr_hexas->SetClippingOff();
      };
      hexa_actor = vtkActor::New();
      hexa_geometry->SetInput(extr_hexas->GetOutput());
      hexa_geometry->Update();
      hexas_pd->DeepCopy(hexa_geometry->GetOutput());
      hexa_mapper->SetInput(hexas_pd);
      hexa_actor = vtkActor::New();
      hexa_actor->SetMapper(hexa_mapper);
      hexa_actor->GetProperty()->SetColor(0,0.7,1);
      getRenderer()->AddActor(hexa_actor);
    };
    
    // wireframe
    extr_vol->SetInput(grid);
    if (ui.checkBoxClip->isChecked()) {
      extr_vol->SetClippingOn();
      extr_vol->SetX(x);
      extr_vol->SetN(n);
    } else {
      extr_vol->SetClippingOff();
    };
    volume_wire_actor = vtkActor::New();
    volume_geometry->SetInput(extr_vol->GetOutput());
    volume_geometry->Update();
    volume_pd->DeepCopy(volume_geometry->GetOutput());
    volume_wire_mapper->SetInput(volume_pd);
    volume_wire_actor->SetMapper(volume_wire_mapper);
    volume_wire_actor->GetProperty()->SetRepresentationToWireframe();
    volume_wire_actor->GetProperty()->SetColor(0,0,1);
    getRenderer()->AddActor(volume_wire_actor);
    
    
    updateStatusBar();
    getRenderWindow()->Render();
  } catch (Error err) {
    err.display();
  };
  unlock();
};

void GuiMainWindow::setPickMode(bool a_UseVTKInteractor,bool a_CellPickerMode)
{
  m_UseVTKInteractor=a_UseVTKInteractor;
  if(a_UseVTKInteractor) ui.checkBox_UseVTKInteractor->setCheckState(Qt::Checked);
  else ui.checkBox_UseVTKInteractor->setCheckState(Qt::Unchecked);
  if(a_CellPickerMode) ui.radioButton_CellPicker->toggle();
  else ui.radioButton_PointPicker->toggle();
//   cout<<"m_UseVTKInteractor="<<m_UseVTKInteractor<<endl;
}

void GuiMainWindow::setUseVTKInteractor(int a_UseVTKInteractor)
{
  m_UseVTKInteractor=a_UseVTKInteractor;
//   cout<<"m_UseVTKInteractor="<<m_UseVTKInteractor<<endl;
}

bool GuiMainWindow::pickPoint(vtkIdType nodeId)
{
  cout<<"pickPoint called with =nodeId"<<nodeId<<endl;
  if (nodeId >= 0 && nodeId<grid->GetNumberOfPoints()) {
    vec3_t x(0,0,0);
    grid->GetPoints()->GetPoint(nodeId, x.data());
    pick_sphere->SetCenter(x.data());
    pick_mapper->SetInput(pick_sphere->GetOutput());
    
    if (pick_actor) {
      getRenderer()->RemoveActor(pick_actor);
      pick_actor->Delete();cout<<"Deleting pick_actor="<<pick_actor<<endl;
      pick_actor = NULL;
    };
    pick_actor = vtkActor::New();cout<<"New pick_actor="<<pick_actor<<endl;
//     pick_actor->VisibilityOn();
    
    pick_actor->SetMapper(pick_mapper);
    pick_actor->GetProperty()->SetRepresentationToSurface();
    pick_actor->GetProperty()->SetColor(0,0,1);
    getRenderer()->AddActor(pick_actor);
    PickedPoint=nodeId;
    return(true);
  }
  else return(false);
}

bool GuiMainWindow::pickCell(vtkIdType cellId)
{
  cout<<"pickCell called with cellId="<<cellId<<endl;
  if (cellId >= 0 && cellId<grid->GetNumberOfCells()) {
    vtkIdType *pts, Npts;
    grid->GetCellPoints(cellId, Npts, pts);
    vec3_t x(0,0,0);
    for (vtkIdType i = 0; i < Npts; ++i) {
      vec3_t xp;
      grid->GetPoints()->GetPoint(pts[i], xp.data());
      x += double(1)/Npts * xp;
    };
    pick_sphere->SetCenter(x.data());
    double R = 1e99;
//     for (vtkIdType i = 0; i < Npts; ++i) {
//       vec3_t xp;
//       grid->GetPoints()->GetPoint(pts[i], xp.data());
//       R = min(R, 0.25*(xp-x).abs());
//     };
    for (vtkIdType i = 0; i < Npts; ++i) {
      vec3_t xp;
      grid->GetPoints()->GetPoint(pts[i], xp.data());
      R = min(R, 0.25*(xp-x).abs());
    };
//     R=getShortestSide(cellId,grid);
    ReferenceSize=R;//Used for text annotations too!
    pick_sphere->SetRadius(R);
    pick_mapper->SetInput(pick_sphere->GetOutput());
    
    if (pick_actor) {
      getRenderer()->RemoveActor(pick_actor);
      pick_actor->Delete();cout<<"Deleting pick_actor="<<pick_actor<<endl;
      pick_actor = NULL;
    };
    pick_actor = vtkActor::New();cout<<"New pick_actor="<<pick_actor<<endl;
    
    pick_actor->SetMapper(pick_mapper);
    pick_actor->GetProperty()->SetRepresentationToSurface();
    pick_actor->GetProperty()->SetColor(1,0,0);
    getRenderer()->AddActor(pick_actor);
    PickedCell=cellId;
    return(true);
  }
  else return(false);
}

void GuiMainWindow::importSTL()
{
  StlReader stl;
  stl();
  updateBoundaryCodes(true);
  updateActors();
  updateStatusBar();
  zoomAll();
};

void GuiMainWindow::importGmsh1Ascii()
{
  GmshReader gmsh;
  gmsh.setV1Ascii();
  gmsh();
  updateBoundaryCodes(true);
  updateActors();
  updateStatusBar();
  zoomAll();
};

void GuiMainWindow::exportGmsh1Ascii()
{
  GmshWriter gmsh;
  gmsh.setV1Ascii();
  gmsh();
};

void GuiMainWindow::importGmsh2Ascii()
{
  GmshReader gmsh;
  gmsh.setV2Ascii();
  gmsh();
  updateBoundaryCodes(true);
  updateActors();
  updateStatusBar();
  zoomAll();
};

void GuiMainWindow::exportGmsh2Ascii()
{
  GmshWriter gmsh;
  gmsh.setV2Ascii();
  gmsh();
};

void GuiMainWindow::exportNeutral()
{
  NeutralWriter neutral;
  neutral();
};

void GuiMainWindow::zoomAll()
{
  getRenderer()->ResetCamera();
  getRenderWindow()->Render();
}

void GuiMainWindow::ZoomOnPickedObject()
{
  if(pick_actor!=NULL)
  {
    getRenderer()->ResetCamera(pick_actor->GetBounds());
    getRenderWindow()->Render();
  }
}

void GuiMainWindow::DeselectAll()
{
  cout<<"void GuiMainWindow::DeselectAll()"<<endl;
  pick_actor->VisibilityOff();
//   goo;
/*  if (pick_actor) {
    cout<<"Terminating pick_actor"<<endl;
    getRenderer()->RemoveActor(pick_actor);
    pick_actor->Delete();
    pick_actor = NULL;
  };
  setPickMode(false,true);
  PickedCell=-1;
  PickedPoint=-1;*/
  updateActors();
}

///@@@  TODO: Should display a window
void GuiMainWindow::Info()
{
  ShowInfo info(ui.radioButton_CellPicker->isChecked(),PickedPoint,PickedCell);
  info();
}

int GuiMainWindow::QuickSave()
{
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
}

void GuiMainWindow::QuickLoad(int a_operation)
{
  QFileInfo fileinfo(current_filename);
  QString l_filename = m_LogDir + fileinfo.completeBaseName() + "_" + QString("%1").arg(a_operation) + ".vtu";
//   cout<<"Loading l_filename="<<l_filename.toLatin1().data()<<endl;
  QuickLoad(l_filename);
  setWindowTitle(current_filename + " - enGrid - " + QString("%1").arg(current_operation) );
}

void GuiMainWindow::Undo()
{
  cout << "Undoing operation " << current_operation << endl;
  current_operation--;
  QuickLoad(current_operation);
  ui.actionRedo->setEnabled(true);
  if(current_operation<=0) ui.actionUndo->setEnabled(false);
};

void GuiMainWindow::Redo()
{
  current_operation++;
  cout << "Redoing operation " << current_operation << endl;
  QuickLoad(current_operation);
  ui.actionUndo->setEnabled(true);
  if(current_operation>=last_operation) ui.actionRedo->setEnabled(false);
};

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
};

void GuiMainWindow::saveBC()
{
  saveBC(current_filename);
};

void GuiMainWindow::openBC(QString a_file)
{
  QString bc_file = a_file + ".bcs";
  QFile file(bc_file);
  bcmap.clear();
  if (file.exists()) {
    file.open(QIODevice::ReadOnly | QIODevice::Text);
    QTextStream f(&file);
    while (!f.atEnd()) {
      QString name, type;
      int i;
      f >> i >> name >> type;
      bcmap[i] = BoundaryCondition(name,type);
    };
  };
}

void GuiMainWindow::saveBC(QString a_file)
{
  QString bc_file = a_file + ".bcs";
  QFile file(bc_file);
  file.open(QIODevice::WriteOnly | QIODevice::Text);
  QTextStream f(&file);
  foreach(int i, all_boundary_codes) {
    BoundaryCondition bc = bcmap[i];
    f << i << " " << bc.getName() << " " << bc.getType() << "\n";
  };
}

///@@@  TODO: I think this should also be a done by a subclass of IOOperation just like for import operations
void GuiMainWindow::open()
{
  current_filename = QFileDialog::getOpenFileName
    (
      NULL,
      "open grid from file",
      getCwd(),
      "VTK unstructured grid files (*.vtu *.VTU)"
    );
  if (!current_filename.isNull()) {
    GuiMainWindow::setCwd(QFileInfo(current_filename).absolutePath());
    EG_VTKSP(vtkXMLUnstructuredGridReader,vtu);
    vtu->SetFileName(current_filename.toAscii().data());
    vtu->Update();
    grid->DeepCopy(vtu->GetOutput());
    createBasicFields(grid, grid->GetNumberOfCells(), grid->GetNumberOfPoints(), false);
    setWindowTitle(current_filename + " - enGrid - " + QString("%1").arg(current_operation) );
    openBC();
    updateBoundaryCodes(true);
    updateActors();
    updateStatusBar();
    zoomAll();
    ResetOperationCounter();
    QuickSave();
  };
};

void GuiMainWindow::save()
{
  cout << current_filename.toAscii().data() << endl;
  if (current_filename == "untitled.vtu") {
    saveAs();
  } else {
    EG_VTKDCC(vtkDoubleArray, cell_VA, grid, "cell_VA");
    for (vtkIdType cellId = 0; cellId < grid->GetNumberOfCells(); ++cellId) {
      cell_VA->SetValue(cellId, GeometryTools::cellVA(grid, cellId, true));
    };
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
  };
};

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
      };
    };
    EG_VTKDCC(vtkDoubleArray, cell_VA, grid, "cell_VA");
    for (vtkIdType cellId = 0; cellId < grid->GetNumberOfCells(); ++cellId) {
      cell_VA->SetValue(cellId, GeometryTools::cellVA(grid, cellId, true));
    };
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
      QuickSave();
    }
  };
};

void GuiMainWindow::QuickSave(QString a_filename)
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
    };
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
};

void GuiMainWindow::QuickLoad(QString a_filename)
{
    cout << "Loading " << a_filename.toAscii().data() << endl;

    if (!a_filename.isNull()) {
//       GuiMainWindow::setCwd(QFileInfo(a_filename).absolutePath());
      EG_VTKSP(vtkXMLUnstructuredGridReader,vtu);
      vtu->SetFileName(a_filename.toAscii().data());
      vtu->Update();
      grid->DeepCopy(vtu->GetOutput());
      createBasicFields(grid, grid->GetNumberOfCells(), grid->GetNumberOfPoints(), false);
      setWindowTitle(a_filename + " - enGrid - " + QString("%1").arg(current_operation) );
      openBC(a_filename);
      updateBoundaryCodes(true);
      updateActors();
      updateStatusBar();
      zoomAll();
    };
  
};

void GuiMainWindow::updateStatusBar()
{
  QString num, txt = "enGrid is currently busy with an operation ...";
  if (!busy) {
    txt = "";
  };
  if (!tryLock()) {
//     status_label->setText(txt);
    status_bar->showMessage(txt);
    status_bar->setToolTip(txt);
    ui.label_node_cell_info->setText(txt);
    return;
  };
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
  };
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
      if      (type_cell == VTK_TRIANGLE)   pick_txt += "triangle";
      else if (type_cell == VTK_QUAD)       pick_txt += "quad";
      else if (type_cell == VTK_TETRA)      pick_txt += "tetrahedron";
      else if (type_cell == VTK_PYRAMID)    pick_txt += "pyramid";
      else if (type_cell == VTK_WEDGE)      pick_txt += "prism";
      else if (type_cell == VTK_HEXAHEDRON) pick_txt += "hexahedron";
      vtkIdType N_pts, *pts;
      grid->GetCellPoints(id_cell, N_pts, pts);
      pick_txt += " [";
      for (int i_pts = 0; i_pts < N_pts; ++i_pts) {
        QString num;
        num.setNum(pts[i_pts]);
        pick_txt += num;
        if (i_pts < N_pts-1) {
          pick_txt += ",";
        };
      };
      pick_txt += "]";
      QString tmp;
      EG_VTKDCC(vtkIntArray, cell_code, grid, "cell_code");
      tmp.setNum(cell_code->GetValue(id_cell));
      pick_txt += " cell_code=" + tmp;
      tmp.setNum(id_cell);
      pick_txt += " id_cell=" + tmp;
    };
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
        };
      };
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
    };
    
    txt += pick_txt;
  }
  
//   status_label->setText(txt);
  status_bar->showMessage(txt);
  status_bar->setToolTip(txt);
  ui.label_node_cell_info->setText(txt);
  unlock();
};

void GuiMainWindow::selectBoundaryCodes()
{
  GuiSelectBoundaryCodes bcodes;
  bcodes.setDisplayBoundaryCodes(display_boundary_codes);
  cout<<"void GuiMainWindow::selectBoundaryCodes(): all_boundary_codes="<<all_boundary_codes<<endl;

  bcodes.setBoundaryCodes(all_boundary_codes);
  bcodes();
  bcodes.getThread().wait();
  bcodes.getSelectedBoundaryCodes(display_boundary_codes);
  updateActors();
};

void GuiMainWindow::updateBoundaryCodes(bool all_on)
{
//   cout<<"void GuiMainWindow::updateBoundaryCodes(bool all_on)"<<endl;
  try {
    all_boundary_codes.clear();
    EG_VTKDCC(vtkIntArray, cell_code, grid, "cell_code");
    for (vtkIdType i = 0; i < grid->GetNumberOfCells(); ++i) {
      int ct = grid->GetCellType(i);
      if ((ct == VTK_TRIANGLE) || (ct == VTK_QUAD)) {
        all_boundary_codes.insert(cell_code->GetValue(i));
      };
    };
    if (all_on) {
      display_boundary_codes.clear();
      foreach (int bc, all_boundary_codes) {
        display_boundary_codes.insert(bc);
      };
    } else {
      QSet<int> dbcs;
      foreach (int bc, display_boundary_codes) {
        if (all_boundary_codes.contains(bc)) {
          dbcs.insert(bc);
        };
      };
      display_boundary_codes.clear();
      foreach (int bc, all_boundary_codes) {
        if (dbcs.contains(bc)) {
          display_boundary_codes.insert(bc);
        };
      };
    };
  } catch (Error err) {
    err.display();
  };
//   cout<<"void GuiMainWindow::updateBoundaryCodes(bool all_on): all_boundary_codes="<<all_boundary_codes<<endl;
};

void GuiMainWindow::normalExtrusion()
{
  GuiNormalExtrusion extr;
  extr();
  updateBoundaryCodes(false);
  updateActors();
};

void GuiMainWindow::setAxesVisibility()
{
  if (ui.actionViewAxes->isChecked()) axes->SetVisibility(1);
  else                                axes->SetVisibility(0);
  getRenderWindow()->Render();
};

void GuiMainWindow::setViewingMode()
{
  if (ui.actionViewOrthogonal->isChecked()) getRenderer()->GetActiveCamera()->ParallelProjectionOn();
  else getRenderer()->GetActiveCamera()->ParallelProjectionOff();
  getRenderWindow()->Render();
};

void GuiMainWindow::ViewNodeIDs()
{
  int N=grid->GetNumberOfPoints();
  cout<<"N="<<N<<endl;
  if (ui.actionViewNodeIDs->isChecked()) {
    cout<<"Activating node ID view"<<endl;
    NodeText_VectorText.resize(N);
    NodeText_PolyDataMapper.resize(N);
    NodeText_Follower.resize(N);
    for(int i=0;i<N;i++){
      NodeText_VectorText[i]=vtkVectorText::New();
      QString tmp;
      tmp.setNum(i);
      NodeText_VectorText[i]->SetText(tmp.toLatin1().data());
      NodeText_PolyDataMapper[i]=vtkPolyDataMapper::New();
      NodeText_PolyDataMapper[i]->SetInputConnection(NodeText_VectorText[i]->GetOutputPort());
      NodeText_Follower[i]=vtkFollower::New();
      NodeText_Follower[i]->SetMapper(NodeText_PolyDataMapper[i]);
      NodeText_Follower[i]->SetScale(ReferenceSize,ReferenceSize,ReferenceSize);
      vec3_t M;
      grid->GetPoint(i,M.data());
      vec3_t tmp_M=M;
      vec3_t OffSet=ReferenceSize*tmp_M.normalise();
      M=M+OffSet;
      NodeText_Follower[i]->AddPosition(M[0],M[1],M[2]);
      NodeText_Follower[i]->SetCamera(getRenderer()->GetActiveCamera());
      NodeText_Follower[i]->GetProperty()->SetColor(0,0,1);
      getRenderer()->AddActor(NodeText_Follower[i]);
    }
  }
  else {
    cout<<"Deactivating node ID view"<<endl;
    for(unsigned int i=0;i<NodeText_Follower.size();i++){
      getRenderer()->RemoveActor(NodeText_Follower[i]);
      NodeText_Follower[i]->Delete();
      NodeText_PolyDataMapper[i]->Delete();
      NodeText_VectorText[i]->Delete();
    }
    NodeText_Follower.clear();
    NodeText_PolyDataMapper.clear();
    NodeText_VectorText.clear();
  }
  
  getRenderWindow()->Render();
};

void GuiMainWindow::ViewCellIDs()
{
  vtkIdType N=grid->GetNumberOfCells();
  cout<<"N="<<N<<endl;
  if (ui.actionViewCellIDs->isChecked()) {
    cout<<"Activating cell ID view"<<endl;
    CellText_VectorText.resize(N);
    CellText_PolyDataMapper.resize(N);
    CellText_Follower.resize(N);
    for(vtkIdType id_cell=0;id_cell<N;id_cell++){
      CellText_VectorText[id_cell]=vtkVectorText::New();
      
      QString tmp;
      
      if(ui.comboBox_CellTextField->currentIndex()==0) {
        tmp.setNum(id_cell);
      }
      else if(ui.comboBox_CellTextField->currentIndex()>0) {
        EG_VTKDCC(vtkIntArray, current_cell_field, grid, ui.comboBox_CellTextField->currentText().toLatin1().data());
        tmp.setNum(current_cell_field->GetValue(id_cell));
      }
      else EG_BUG;
      
      CellText_VectorText[id_cell]->SetText(tmp.toLatin1().data());
      CellText_PolyDataMapper[id_cell]=vtkPolyDataMapper::New();
      CellText_PolyDataMapper[id_cell]->SetInputConnection(CellText_VectorText[id_cell]->GetOutputPort());
      CellText_Follower[id_cell]=vtkFollower::New();
      CellText_Follower[id_cell]->SetMapper(CellText_PolyDataMapper[id_cell]);
      CellText_Follower[id_cell]->SetScale(ReferenceSize,ReferenceSize,ReferenceSize);
      vtkIdType N_pts,*pts;
      grid->GetCellPoints(id_cell,N_pts,pts);
      vec3_t Center(0,0,0);
      for(int p=0;p<N_pts;p++)
      {
        vec3_t M;
        grid->GetPoint(pts[p],M.data());
        Center+=M.data();
      }
      vec3_t OffSet=ReferenceSize*triNormal(grid,pts[0],pts[1],pts[2]).normalise();
      Center=1.0/(double)N_pts*Center+OffSet;
      CellText_Follower[id_cell]->AddPosition(Center[0],Center[1],Center[2]);
      CellText_Follower[id_cell]->SetCamera(getRenderer()->GetActiveCamera());
      CellText_Follower[id_cell]->GetProperty()->SetColor(1,0,0);
      getRenderer()->AddActor(CellText_Follower[id_cell]);
    }
  }
  else {
    cout<<"Deactivating cell ID view"<<endl;
    for(vtkIdType id_cell=0;id_cell<(vtkIdType)CellText_Follower.size();id_cell++){
      getRenderer()->RemoveActor(CellText_Follower[id_cell]);
      CellText_Follower[id_cell]->Delete();
      CellText_PolyDataMapper[id_cell]->Delete();
      CellText_VectorText[id_cell]->Delete();
    }
    CellText_Follower.clear();
    CellText_PolyDataMapper.clear();
    CellText_VectorText.clear();
  }
  
  getRenderWindow()->Render();
};

void GuiMainWindow::addVtkTypeInfo()
{
  EG_VTKSP(vtkIntArray, vtk_type);
  vtk_type->SetName("vtk_type");
  vtk_type->SetNumberOfValues(grid->GetNumberOfCells());
  for (vtkIdType cellId = 0; cellId < grid->GetNumberOfCells(); ++cellId) {
    vtk_type->SetValue(cellId, grid->GetCellType(cellId));
  };
  grid->GetCellData()->AddArray(vtk_type);
};

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
};

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
};

vtkIdType GuiMainWindow::getPickedCell()
{
  vtkIdType picked_cell = -1;
  if (THIS->grid->GetNumberOfCells() > 0) {
    THIS->bcodes_filter->Update();
    EG_VTKDCC(vtkLongArray_t, cell_index, THIS->bcodes_filter->GetOutput(), "cell_index");
    
    vtkIdType cellId;
    if(m_UseVTKInteractor) cellId = THIS->CellPicker->GetCellId();
    else cellId = PickedCell;
    
    if (cellId >= 0) {
      if (cellId < THIS->bcodes_filter->GetOutput()->GetNumberOfCells()) {
        picked_cell = cell_index->GetValue(cellId);
      };
    };
  };
  return picked_cell;
};

vtkIdType GuiMainWindow::getPickedPoint()
{
  vtkIdType picked_point = -1;
  if (THIS->grid->GetNumberOfCells() > 0) {
    THIS->bcodes_filter->Update();
    
    vtkIdType pointId;
    if(m_UseVTKInteractor) pointId = THIS->PointPicker->GetPointId();
    else pointId = PickedPoint;
    
    if (pointId >= 0) {
      picked_point = pointId;
    }
  };
  return picked_point;
};

void GuiMainWindow::changeSurfaceOrientation()
{
  for (vtkIdType cellId = 0; cellId < grid->GetNumberOfCells(); ++cellId) {
    vtkIdType Npts, *pts;
    grid->GetCellPoints(cellId, Npts, pts);
    QVector<vtkIdType> nodes(Npts);
    for (vtkIdType j = 0; j < Npts; ++j) nodes[j]          = pts[j];
    for (vtkIdType j = 0; j < Npts; ++j) pts[Npts - j - 1] = nodes[j];
  };
  updateActors();
};

void GuiMainWindow::checkSurfaceOrientation()
{
  CorrectSurfaceOrientation corr_surf;
  vtkIdType picked_cell = getPickedCell();
  if (picked_cell >= 0) {
    corr_surf.setStart(picked_cell);
  };
  corr_surf();
  updateActors();
};

void GuiMainWindow::improveAspectRatio()
{
  GuiImproveAspectRatio impr_ar;
  impr_ar();
  updateActors();
};

void GuiMainWindow::exportAsciiStl()
{
  StlWriter stl;
  stl();
};

void GuiMainWindow::exportBinaryStl()
{
};

void GuiMainWindow::periodicUpdate()
{
  Operation::collectGarbage(); 
  updateStatusBar();
};

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
};

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
};

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
};

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
};

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
};

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
};

void GuiMainWindow::callFixSTL()
{
  FixSTL *fix;
  fix = new FixSTL();
  fix->setGui();
  (*fix)();
  updateBoundaryCodes(false);
  updateActors();
};

void GuiMainWindow::editBoundaryConditions()
{
  GuiEditBoundaryConditions editbcs;
  editbcs.setBoundaryCodes(all_boundary_codes);
  editbcs.setMap(&bcmap);
  editbcs();
};

void GuiMainWindow::configure()
{
  {
    // Just to create initial entries in the settings file 
    // so that the options menu isn't empty at first start.
    GridSmoother tmp01;
    GuiCreateBoundaryLayer tmp02;
  };
  GuiSettingsViewer settings(&qset);
  settings.exec();
};

void GuiMainWindow::about()
{
  QMessageBox box(this);
  
  QString title="ENGRID";
  QString version = QString("version ") + ENGRID_VERSION;
  
/*  if (version == "version CVS") {
  };*/
  
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
  
};

void GuiMainWindow::getAllBoundaryCodes(QSet<int> &bcs)
{
  bcs.clear();
  foreach (int bc, all_boundary_codes) {
    bcs.insert(bc);
  };
};

void GuiMainWindow::getDisplayBoundaryCodes(QSet<int> &bcs)
{
  bcs.clear();
  foreach (int bc, display_boundary_codes) {
    bcs.insert(bc);
  };
};
