// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008-2014 enGits GmbH                                      +
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
#include "guimainwindow.h"
#include "engrid.h"
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
#include "triangularcadinterface.h"
#include "cgaltricadinterface.h"

#include <qmessagebox.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkProperty.h>
#include <vtkCallbackCommand.h>
#include <vtkCamera.h>
#include <vtkSignedCharArray.h>
#include <vtkTextActor.h>
#include <vtkVectorText.h>
#include <vtkFollower.h>
#include <vtkLookupTable.h>
#include <vtkScalarBarActor.h>
#include <vtkFileOutputWindow.h>
#include <QVTKInteractor.h>

#include <QFileDialog>
#include <QFileSystemWatcher>
#include <QFileInfo>
#include <QAction>
#include <QStatusBar>
#include <QInputDialog>

#include <stdlib.h>
#include <stdio.h>

#if !defined( _WIN32 ) //required for "dup" and "dup2" on POSIX systems
#include <unistd.h>
#endif

#include "geometrytools.h"
#include "engrid_version.h"

using namespace GeometryTools;

#include "guisettingsviewer.h"
#include "guitransform.h"
#include "egvtkinteractorstyle.h"
#include "showinfo.h"
#include "engrid_version.h"

#include <csignal>

QString GuiMainWindow::m_cwd = ".";
QSettings GuiMainWindow::m_qset("enGits", QString("enGrid-") + ENGRID_VERSION_STRING);
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
  //QMenuBar *menubar = new QMenuBar(0);
  ui.setupUi(this);
#ifdef __APPLE__
  QMenuBar *menubar = ui.menubar;
  menubar->setNativeMenuBar(true);
#endif
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

  m_StatusInfoLabel = new QLabel(this);
  statusBar()->addWidget(m_StatusInfoLabel);
  m_StatusInfoLabel->setText("");

  m_StatusProgressBar = new QProgressBar(this);
  statusBar()->addWidget(m_StatusProgressBar);

  m_StatusLabel = new QLabel(this);
  statusBar()->addWidget(m_StatusLabel);

  QString txt = "0 volume cells (0 tetras, 0 hexas, 0 pyramids, 0 prisms, 0 polys), ";
  txt += "0 surface cells (0 triangles, 0 quads, 0 polys), 0 nodes";
  m_StatusLabel->setText(txt);
  ui.label_node_cell_info->setText(txt);

  m_OriginalCoutBuffer = cout.rdbuf();
  cout.rdbuf(&m_CoutBuffer);
  m_Log = true;

  m_Busy = false;

  setPickMode(true,true);
  m_PickedPoint = -1;
  m_PickedCell = -1;

  updateStatusBar();

  connect(&m_GarbageTimer, SIGNAL(timeout()), this, SLOT(periodicUpdate()));
  m_GarbageTimer.start(1000);

  connect(&m_LogTimer, SIGNAL(timeout()), this, SLOT(updateOutput()));
  m_LogTimer.start(1000);

  bool exp_features=false;
  getSet("General","enable experimental features",false,exp_features);
  getSet("General","enable undo+redo",false,m_undo_redo_enabled);
  bool undo_redo_mode;
  getSet("General","use RAM for undo+redo operations",false,undo_redo_mode);
  getSet("General", "open last used file on startup", false, m_open_last);

  ui.actionMirrorMesh->setEnabled(exp_features);
  ui.actionCreateHexCore->setEnabled(exp_features);
  ui.actionImportFluentCase->setEnabled(exp_features);

  m_ReferenceSize=0.2;

  ui.doubleSpinBox_HueMin->setValue(0.667);
  ui.doubleSpinBox_HueMax->setValue(0);

  // egvtkInteractorStyle *style = egvtkInteractorStyle::New();
  // getInteractor()->SetInteractorStyle(style);
  // style->Delete();

  // initialise XML document
  m_XmlHandler = new XmlHandler("engridcase");
//   this->resetXmlDoc();

  m_SolverIndex = 0;

  readRecentFiles();

  // load plugins
  QString plugin_path;
  getSet("General", "plugin path", "/usr/lib/engrid", plugin_path);
  QDir plugin_dir(plugin_path);
  m_PluginOperations.clear();
  foreach (QString file_name, plugin_dir.entryList(QDir::Files)) {
    if (file_name.right(3) == ".so") {
      cout << qPrintable(plugin_dir.absoluteFilePath(file_name)) << endl;
      QPluginLoader loader(plugin_dir.absoluteFilePath(file_name));
      QObject *qobject = loader.instance();
      if (!qobject) {
        cout << "an error occurred while loading the plugins:\n";
        cout << qPrintable(loader.errorString()) << "\n" << endl;
      }
      if (Operation *operation = qobject_cast<Operation*>(qobject)) {
        //operation->setLockGui();
        QAction *action = new QAction(operation->getMenuText(), this);
        connect(action, SIGNAL(triggered()), this, SLOT(pluginCalled()));
        m_PluginOperations[action] = operation;
        ui.menuPlugins->addAction(action);
      }
    }
  }

  m_EscAction = new QAction("escape", this);
  addAction(m_EscAction);
  m_EscAction->setShortcut(QKeySequence(Qt::Key_Escape));
  connect(m_EscAction, SIGNAL(triggered()), this, SLOT(onEsc()));

  m_UniCadInterface = NULL;
}
//end of GuiMainWindow::GuiMainWindow() : QMainWindow(NULL)

void GuiMainWindow::pluginCalled()
{
  QAction *action = qobject_cast<QAction*>(QObject::sender());
  if (action) {
    Operation *operation = m_PluginOperations[action];
    operation->operator()();
  }
}

void GuiMainWindow::resetXmlDoc()
{
  m_XmlHandler->resetXmlDoc();
}

GuiMainWindow::~GuiMainWindow()
{
  writeRecentFiles();

  m_qset.setValue("GuiMainWindow", this->geometry());
  m_qset.setValue("dockWidget_states", this->saveState());

  cout.rdbuf(m_OriginalCoutBuffer);

  delete m_XmlHandler;
}

void GuiMainWindow::setupVtk()
{

// avoid VTK pop-up window on Windows
#ifdef WIN32
  vtkFileOutputWindow *w = vtkFileOutputWindow::New();
  QString vtk_log_file = m_qset.value("tmp_directory").toString() + "/enGrid-vtk-errors.txt";
  w->SetFileName(qPrintable(vtk_log_file));
  vtkOutputWindow::SetInstance(w);
  w->Delete();
#endif

  // colour settings
  getSet("Colours", "'A' faces (1-red)",   0.5, m_ColAR);
  getSet("Colours", "'A' faces (2-green)", 1.0, m_ColAG);
  getSet("Colours", "'A' faces (3-blue)",  0.5, m_ColAB);

  getSet("Colours", "'B' faces (1-red)",   1.0, m_ColBR);
  getSet("Colours", "'B' faces (2-green)", 1.0, m_ColBG);
  getSet("Colours", "'B' faces (3-blue)",  0.5, m_ColBB);

  getSet("Colours", " tetras (1-red)",     1.0, m_ColTetraR);
  getSet("Colours", " tetras (2-green)",   0.0, m_ColTetraG);
  getSet("Colours", " tetras (3-blue)",    0.0, m_ColTetraB);

  getSet("Colours", " prisms (1-red)",     0.0, m_ColPrismR);
  getSet("Colours", " prisms (2-green)",   1.0, m_ColPrismG);
  getSet("Colours", " prisms (3-blue)",    0.0, m_ColPrismB);

  getSet("Colours", " pyramids (1-red)",   1.0, m_ColPyraR);
  getSet("Colours", " pyramids (2-green)", 1.0, m_ColPyraG);
  getSet("Colours", " pyramids (3-blue)",  0.0, m_ColPyraB);

  getSet("Colours", " hexes (1-red)",      0.0, m_ColHexR);
  getSet("Colours", " hexes (2-green)",    0.7, m_ColHexG);
  getSet("Colours", " hexes (3-blue)",     1.0, m_ColHexB);

  getSet("Colours", " polys (1-red)",      0.7, m_ColPolyR);
  getSet("Colours", " polys (2-green)",    1.0, m_ColPolyG);
  getSet("Colours", " polys (3-blue)",     1.0, m_ColPolyB);

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
  m_SurfaceFilter     = vtkDataSetSurfaceFilter::New();
  m_SurfaceMapper     = vtkPolyDataMapper::New();
  m_BCodesFilter      = vtkEgBoundaryCodesFilter::New();
  m_LookupTable       = vtkLookupTable::New();
  m_SurfaceActor      = vtkActor::New();
  m_LegendActor       = vtkScalarBarActor::New();
  //
  m_BCodesFilter->SetBoundaryCodes(m_DisplayBoundaryCodes);
  m_BCodesFilter->SetInputData(m_Grid);
  m_SurfaceFilter->SetInputConnection(m_BCodesFilter->GetOutputPort());
  m_SurfaceMapper->SetInputConnection(m_SurfaceFilter->GetOutputPort());
  m_SurfaceMapper->SetLookupTable(m_LookupTable);
  m_SurfaceActor->GetProperty()->SetRepresentationToSurface();
  m_SurfaceActor->GetProperty()->SetColor(m_ColAR, m_ColAG, m_ColAB);
  m_SurfaceActor->SetBackfaceProperty(m_BackfaceProperty);
  m_SurfaceActor->GetBackfaceProperty()->SetColor(m_ColBR, m_ColBG, m_ColBB);
  m_SurfaceActor->GetProperty()->EdgeVisibilityOn();
  m_SurfaceActor->GetProperty()->SetEdgeColor(0,0,1);
  m_SurfaceActor->SetMapper(m_SurfaceMapper);
  getRenderer()->AddActor(m_SurfaceActor);
  m_SurfaceActor->SetVisibility(1);
  m_LegendActor->SetLookupTable(m_LookupTable);
  getRenderer()->AddActor(m_LegendActor);
  m_LegendActor->SetVisibility(0);

  // tetra pipline
  m_ExtrTetras   = vtkEgExtractVolumeCells::New();
  m_TetraActor   = vtkActor::New();
  m_TetraGeometry = vtkDataSetSurfaceFilter::New();
  m_TetraMapper  = vtkPolyDataMapper::New();
  //
  m_ExtrTetras->SetInputData(m_Grid);
  m_ExtrTetras->SetAllOff();
  m_ExtrTetras->SetTetrasOn();;
  m_TetraGeometry->SetInputConnection(m_ExtrTetras->GetOutputPort());
  m_TetraMapper->SetInputConnection(m_TetraGeometry->GetOutputPort());
  m_TetraActor->SetMapper(m_TetraMapper);
  m_TetraActor->GetProperty()->SetColor(m_ColTetraR, m_ColTetraG, m_ColTetraB);
  m_TetraActor->GetProperty()->EdgeVisibilityOn();
  m_TetraActor->GetProperty()->SetEdgeColor(0,0,1);
  getRenderer()->AddActor(m_TetraActor);
  m_TetraActor->SetVisibility(0);

  // pyramid pipeline
  m_PyramidActor   = vtkActor::New();
  m_ExtrPyramids   = vtkEgExtractVolumeCells::New();
  m_PyramidGeometry = vtkDataSetSurfaceFilter::New();
  m_PyramidMapper  = vtkPolyDataMapper::New();
  //
  m_ExtrPyramids->SetInputData(m_Grid);
  m_ExtrPyramids->SetAllOff();
  m_ExtrPyramids->SetPyramidsOn();
  m_PyramidGeometry->SetInputConnection(m_ExtrPyramids->GetOutputPort());
  m_PyramidMapper->SetInputConnection(m_PyramidGeometry->GetOutputPort());
  m_PyramidActor->SetMapper(m_PyramidMapper);
  m_PyramidActor->GetProperty()->SetColor(m_ColPyraR, m_ColPyraG, m_ColPyraB);
  m_PyramidActor->GetProperty()->EdgeVisibilityOn();
  m_PyramidActor->GetProperty()->SetEdgeColor(0,0,1);
  getRenderer()->AddActor(m_PyramidActor);
  m_PyramidActor->SetVisibility(0);

  // wedge pipeline
  m_WedgeActor   = vtkActor::New();
  m_ExtrWedges   = vtkEgExtractVolumeCells::New();
  m_WedgeGeometry = vtkDataSetSurfaceFilter::New();
  m_WedgeMapper  = vtkPolyDataMapper::New();
  //
  m_ExtrWedges->SetInputData(m_Grid);
  m_ExtrWedges->SetAllOff();
  m_ExtrWedges->SetWedgesOn();
  m_WedgeGeometry->SetInputConnection(m_ExtrWedges->GetOutputPort());
  m_WedgeMapper->SetInputConnection(m_WedgeGeometry->GetOutputPort());
  m_WedgeActor->SetMapper(m_WedgeMapper);
  m_WedgeActor->GetProperty()->SetColor(m_ColPrismR, m_ColPrismG, m_ColPrismB);
  m_WedgeActor->GetProperty()->EdgeVisibilityOn();
  m_WedgeActor->GetProperty()->SetEdgeColor(0,0,1);
  getRenderer()->AddActor(m_WedgeActor);
  m_WedgeActor->SetVisibility(0);

  // hexa pipeline
  m_HexaActor   = vtkActor::New();
  m_ExtrHexes   = vtkEgExtractVolumeCells::New();
  m_HexaGeometry = vtkDataSetSurfaceFilter::New();
  m_HexaMapper  = vtkPolyDataMapper::New();
  //
  m_ExtrHexes->SetInputData(m_Grid);
  m_ExtrHexes->SetAllOff();
  m_ExtrHexes->SetHexesOn();
  m_HexaGeometry->SetInputConnection(m_ExtrHexes->GetOutputPort());
  m_HexaMapper->SetInputConnection(m_HexaGeometry->GetOutputPort());
  m_HexaActor->SetMapper(m_HexaMapper);
  m_HexaActor->GetProperty()->SetColor(m_ColHexR, m_ColHexG, m_ColHexB);
  m_HexaActor->GetProperty()->EdgeVisibilityOn();
  m_HexaActor->GetProperty()->SetEdgeColor(0,0,1);
  getRenderer()->AddActor(m_HexaActor);
  m_HexaActor->SetVisibility(0);

  // polyhedra pipeline
  m_PolyhedraActor    = vtkActor::New();
  m_ExtrPolyhedra     = vtkEgExtractVolumeCells::New();
  m_PolyhedraGeometry = vtkDataSetSurfaceFilter::New();
  m_PolyhedraMapper   = vtkPolyDataMapper::New();
  //
  m_ExtrPolyhedra->SetInputData(m_Grid);
  m_ExtrPolyhedra->SetAllOff();
  m_ExtrPolyhedra->SetPolysOn();
  m_PolyhedraGeometry->SetInputConnection(m_ExtrPolyhedra->GetOutputPort());
  m_PolyhedraMapper->SetInputConnection(m_PolyhedraGeometry->GetOutputPort());
  m_PolyhedraActor->SetMapper(m_PolyhedraMapper);
  m_PolyhedraActor->GetProperty()->SetColor(m_ColPolyR, m_ColPolyG, m_ColPolyB);
  m_PolyhedraActor->GetProperty()->EdgeVisibilityOn();
  m_PolyhedraActor->GetProperty()->SetEdgeColor(0,0,1);
  getRenderer()->AddActor(m_PolyhedraActor);
  m_PolyhedraActor->SetVisibility(0);

  // picker stuff
  m_PickSphere  = vtkSphereSource::New();
  m_PickMapper  = vtkPolyDataMapper::New();
  m_PickActor   = vtkActor::New();
  m_CellPicker  = vtkCellPicker::New();
  m_PointPicker = vtkPointPicker::New();

  m_PickSphere->SetRadius(0.25); //in case the user starts picking points instead of cells
  m_PickMapper->SetInputConnection(m_PickSphere->GetOutputPort());
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

  viewFront();

  // Create a new QVTKInteractor instance
  //vtkSmartPointer<QVTKInteractor> interactor = vtkSmartPointer<QVTKInteractor>::New();

  // Set the interactor for your QVTKOpenGLNativeWidget or QVTKOpenGLStereoWidget
  //ui.centralwidget->GetRenderWindow()->SetInteractor(interactor);

  // Initialize the interactor
  //interactor->Initialize();

  // Render the scene
  //ui.myVtkWidget->GetRenderWindow()->Render();


  // ui.centralwidget->setFocusPolicy(Qt::StrongFocus);
  // vtkSmartPointer<vtkRenderWindowInteractor> interactor = ui.centralwidget->GetInteractor();
  // interactor->Initialize();
  //interactor->Start();


}

void GuiMainWindow::updateOutput()
{
  if (m_Log) {
    QString txt = m_CoutBuffer.str().c_str();
    m_CoutBuffer.str("");
    if (txt.right(1) == "\n") {
      txt = txt.left(txt.size()-1);
    }
    if (txt.size() > 0) {
      ui.textEditOutput->append(txt);
    }
  }
}

void GuiMainWindow::exit()
{
  QCoreApplication::exit();
}

vtkRenderWindow* GuiMainWindow::getRenderWindow()
{
  return ui.centralwidget->GetRenderWindow();
}

vtkRenderer* GuiMainWindow::getRenderer()
{
  return m_Renderer;
}

QVTKInteractor* GuiMainWindow::getInteractor()
{
  return ui.centralwidget->GetInteractor();
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
  m_ExtrTetras->Setx(txt.toDouble());
  m_ExtrPyramids->Setx(txt.toDouble());
  m_ExtrWedges->Setx(txt.toDouble());
  m_ExtrHexes->Setx(txt.toDouble());
  m_ExtrPolyhedra->Setx(txt.toDouble());
}

void GuiMainWindow::setClipY(const QString &txt)
{
  m_ExtrTetras->Sety(txt.toDouble());
  m_ExtrPyramids->Sety(txt.toDouble());
  m_ExtrWedges->Sety(txt.toDouble());
  m_ExtrHexes->Sety(txt.toDouble());
  m_ExtrPolyhedra->Sety(txt.toDouble());
}

void GuiMainWindow::setClipZ(const QString &txt)
{
  m_ExtrTetras->Setz(txt.toDouble());
  m_ExtrPyramids->Setz(txt.toDouble());
  m_ExtrWedges->Setz(txt.toDouble());
  m_ExtrHexes->Setz(txt.toDouble());
  m_ExtrPolyhedra->Setz(txt.toDouble());
}

void GuiMainWindow::setClipNX(const QString &txt)
{
  m_ExtrTetras->Setnx(txt.toDouble());
  m_ExtrPyramids->Setnx(txt.toDouble());
  m_ExtrWedges->Setnx(txt.toDouble());
  m_ExtrHexes->Setnx(txt.toDouble());
  m_ExtrPolyhedra->Setnx(txt.toDouble());
}

void GuiMainWindow::setClipNY(const QString &txt)
{
  m_ExtrTetras->Setny(txt.toDouble());
  m_ExtrPyramids->Setny(txt.toDouble());
  m_ExtrWedges->Setny(txt.toDouble());
  m_ExtrHexes->Setny(txt.toDouble());
  m_ExtrPolyhedra->Setny(txt.toDouble());
}

void GuiMainWindow::setClipNZ(const QString &txt)
{
  m_ExtrTetras->Setnz(txt.toDouble());
  m_ExtrPyramids->Setnz(txt.toDouble());
  m_ExtrWedges->Setnz(txt.toDouble());
  m_ExtrHexes->Setnz(txt.toDouble());
  m_ExtrPolyhedra->Setnz(txt.toDouble());
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
  } else {
    m_SurfaceActor->SetVisibility(0);
  }
}

void GuiMainWindow::updateVolumeActors(bool forced)
{
  if (ui.checkBoxVolume->isChecked()) {
    if (ui.checkBoxTetra->isChecked()) {
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
      m_TetraActor->SetVisibility(0);
    }
    if (ui.checkBoxPyramid->isChecked()) {
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
      m_PyramidActor->SetVisibility(0);
    }
    if (ui.checkBoxWedge->isChecked()) {
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
      m_WedgeActor->SetVisibility(0);
    }
    if (ui.checkBoxHexa->isChecked()) {
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
      m_HexaActor->SetVisibility(0);
    }
    if (ui.checkBoxPoly->isChecked()) {
      if (ui.checkBoxClip->isChecked()) {
        m_ExtrPolyhedra->SetClippingOn();
      } else {
        m_ExtrPolyhedra->SetClippingOff();
      }
      if (forced) {
        m_PolyhedraGeometry->Update();
      }
      m_PolyhedraActor->SetVisibility(1);
    } else {
      m_PolyhedraActor->SetVisibility(0);
    }

  } else {
    m_TetraActor->VisibilityOff();
    m_PyramidActor->VisibilityOff();
    m_WedgeActor->VisibilityOff();
    m_HexaActor->VisibilityOff();
    m_PolyhedraActor->VisibilityOff();
  }
}

void GuiMainWindow::updateActors(bool forced)
{
//   qDebug()<<"QApplication::setOverrideCursor(QCursor(Qt::WaitCursor)); called()";
  QApplication::setOverrideCursor(QCursor(Qt::WaitCursor));

  //if (!tryLock()) return;
  try {
    m_Axes->SetInputData(m_Grid);
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
    EG_GET_CELL(id_cell, m_Grid);
    vec3_t x(0,0,0);
    for (vtkIdType i = 0; i < num_pts; ++i) {
      vec3_t xp;
      m_Grid->GetPoints()->GetPoint(pts[i], xp.data());
      x += double(1)/num_pts * xp;
    }
    m_PickSphere->SetCenter(x.data());
    double R = 1e99;
    for (vtkIdType i = 0; i < num_pts; ++i) {
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
  //FIXME: emits an error if no file is imported, so check if there is a valid file
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
      QString l_filename = m_qset.value("tmp_directory").toString() + fileinfo.completeBaseName() + "_" + QString("%1").arg(m_CurrentOperation);
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
    QString l_filename = m_qset.value("tmp_directory").toString() + fileinfo.completeBaseName() + "_" + QString("%1").arg(a_operation) + ".egc";
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
}

void GuiMainWindow::setXmlSection(QString name, QString contents)
{
  m_XmlHandler->setXmlSection(name,contents);
}

void GuiMainWindow::openPhysicalBoundaryConditions()
{
  m_PhysicalBoundaryConditionsMap.clear();
  QString buffer = getXmlSection("engrid/physical").trimmed();
  QStringList lines = buffer.split("\n");
  foreach (QString line, lines) {
    line = line.trimmed();
    QStringList parts = line.split(";", QString::SkipEmptyParts);
    if (parts.size() > 0) {
      QStringList words = parts[0].split(" ", QString::SkipEmptyParts);
      int index = words[0].trimmed().toInt();
      QString name = words[1].trimmed();
      QString type = words[2].trimmed();
      if ((name != "") && (type != "")) {
        PhysicalBoundaryCondition PBC;
        PBC.setName(name);
        PBC.setType(type);
        if (PBC.getNumVars() == parts.size() - 1) {
          PBC.setIndex(index);
          for (int i = 0; i < PBC.getNumVars(); ++i) {
            QStringList words = parts[i+1].split("=");
            PBC.setValueFromString(i, words[1].trimmed());
          }
          m_PhysicalBoundaryConditionsMap[name] = PBC;
        }
      }
    }
  }
}

void GuiMainWindow::savePhysicalBoundaryConditions()
{
  QString buffer("");
  QTextStream f(&buffer, QIODevice::WriteOnly);
  foreach (PhysicalBoundaryCondition PBC, m_PhysicalBoundaryConditionsMap) {
    f << PBC.xmlText() << "\n";
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
      if (i > 0) {
        m_bcmap[i] = BoundaryCondition(name,type,i);
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
  if (m_Grid->GetPointData()->GetArray("node_meshdensity_current")) {
    m_Grid->GetPointData()->RemoveArray("node_meshdensity_current");
  }
  if (m_Grid->GetCellData()->GetArray("cell_VA")) {
    m_Grid->GetCellData()->RemoveArray("cell_VA");
  }
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
  QString file_name = QFileDialog::getOpenFileName(NULL, "open grid from file", getCwd(), "enGrid case files (*.egc *.EGC);; legacy grid files(*.vtu *.VTU)");
  if (!file_name.isNull()) {
    this->open(file_name);
  }
}

void GuiMainWindow::open(QString file_name, bool update_current_filename)
{
  cout << "Opening " << qPrintable(file_name) << endl;

  //QFileInfo file_info(file_name);
  bool no_case_file = false;
  QString file_extension = getExtension(file_name);
  QString grid_file_name = file_name;
  if (file_extension.toLower() == "vtu") {
    no_case_file = true;
    grid_file_name = stripFromExtension(file_name);
  }
  if (!no_case_file) {
    if(!m_XmlHandler->openXml(file_name)) {
      QMessageBox::critical(this, tr("Open failed"), tr("Error reading enGrid case file:\n%1").arg(file_name));
      return;
    }
  }
  if(update_current_filename) {
    GuiMainWindow::setCwd(QFileInfo(file_name).absolutePath());
  }
  resetCadInterfaces();
  {
    QFile geo_file(file_name + ".geo.vtu");
    if (geo_file.exists()) {
      openGrid(file_name + ".geo");
      storeCadInterfaces(true);
    }
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
    m_qset.setValue("LatestFile",file_name);
    resetOperationCounter();
    quickSave();
  }
}

QString GuiMainWindow::saveAs(QString file_name, bool update_current_filename)
{
  QString buffer = m_XmlHandler->getBuffer(0);
  if(update_current_filename) {
    QApplication::setOverrideCursor(QCursor(Qt::WaitCursor));
  }
  QFileInfo file_info(file_name);
  if (file_info.suffix().toLower() != "egc") {
    file_name += ".egc";
  }
  if(update_current_filename) {
    GuiMainWindow::setCwd(file_info.absolutePath());
    m_CurrentFilename = file_name;
  }
  if(!saveGrid(m_Grid, file_name)) {
    QMessageBox::critical(this, QObject::tr("Save failed"), QObject::tr("The grid could not be saved as:\n%1").arg(file_name));
  }
  saveBC();
  savePhysicalBoundaryConditions();
  m_XmlHandler->saveXml(file_name);
  setWindowTitle(m_CurrentFilename + " - enGrid - " + QString("%1").arg(m_CurrentOperation) );
  setUnsaved(false);
  if(update_current_filename) {
    QApplication::restoreOverrideCursor();
  }
  if(update_current_filename) {
    this->addRecentFile(file_name,QDateTime::currentDateTime());
    m_qset.setValue("LatestFile",file_name);
  }
  return(file_name);
}

void GuiMainWindow::save()
{
  if ( m_CurrentFilename == "untitled.egc" || m_UnSaved ) {

    //FIXME: This is more of a hack than a fix...
    if(GuiMainWindow::tryLock()) {
      GuiMainWindow::unlock(); //must unlock before continuing.
      saveAs();
    } else {
      cout <<endl
           << "WARNING: Please save the project before running the requested operation "
              "or after the current operation is complete."
           <<endl;
    }
  } else {
    saveAs(m_CurrentFilename);
  }
}

void GuiMainWindow::saveAs()
{
  QApplication::restoreOverrideCursor();
  //saveGrid(m_Grid, m_CurrentFilename + ".geo");
  bool geo_file_exists = false;
  QString old_geo_file = m_CurrentFilename + ".geo.vtu";
  {
    if (QFileInfo(old_geo_file).exists()) {
      geo_file_exists = true;
    }
  }
  QFileDialog dialog(NULL, "write case to file", getCwd(), "enGrid case files (*.egc)");
  QFileInfo file_info(m_CurrentFilename);
  dialog.selectFile(file_info.completeBaseName() + ".egc");
  dialog.setAcceptMode(QFileDialog::AcceptSave);
  dialog.setConfirmOverwrite(true);
  if (dialog.exec()) {
    QApplication::setOverrideCursor(QCursor(Qt::WaitCursor));
    QStringList selected_files = dialog.selectedFiles();
    QString file_name = selected_files[0];
    if (!file_name.isNull()) {
      QString new_geo_file = file_name + ".geo.vtu";
      {
        QFile file(new_geo_file);
        file.remove();
      }
      QFile geo_file(old_geo_file);
      geo_file.copy(new_geo_file);
      saveAs(file_name);
      //for the undo/redo operations
      resetOperationCounter();
      quickSave();
    }
  }
  QApplication::restoreOverrideCursor();
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
  vtkIdType Nplgs  = 0;
  vtkIdType Ntets  = 0;
  vtkIdType Npyras = 0;
  vtkIdType Nprism = 0;
  vtkIdType Nhexas = 0;
  vtkIdType Npolys = 0;
  for (vtkIdType i = 0; i < Ncells; ++i) {
    int ct = m_Grid->GetCellType(i);
    if      (ct == VTK_TRIANGLE)   ++Ntris;
    else if (ct == VTK_QUAD)       ++Nquads;
    else if (ct == VTK_POLYGON)    ++Nplgs;
    else if (ct == VTK_TETRA)      ++Ntets;
    else if (ct == VTK_WEDGE)      ++Nprism;
    else if (ct == VTK_PYRAMID)    ++Npyras;
    else if (ct == VTK_HEXAHEDRON) ++Nhexas;
    else if (ct == VTK_POLYHEDRON) ++Npolys;
  }
  num.setNum(Ntets + Npyras + Nprism + Nhexas + Npolys); txt += num + " volume cells(";
  num.setNum(Ntets);  txt += num + " tetras, ";
  num.setNum(Npyras); txt += num + " pyramids, ";
  num.setNum(Nprism); txt += num + " prisms, ";
  num.setNum(Nhexas); txt += num + " hexas, ";
  num.setNum(Npolys); txt += num + " polys), ";
  num.setNum(Ntris + Nquads + Nplgs); txt += num + " surface cells(";
  num.setNum(Ntris);  txt += num + " triangles, ";
  num.setNum(Nquads); txt += num + " quads, ";
  num.setNum(Nplgs);  txt += num + " polys), ";
  num.setNum(Nnodes); txt += num + " nodes";

  if(ui.radioButton_CellPicker->isChecked())
  {
    QString pick_txt = ", picked cell: ";
    vtkIdType id_cell = m_PickedCell;
    if (id_cell < 0 || id_cell>=m_Grid->GetNumberOfCells()) {
      pick_txt += "no cell picked";
    } else {
      EG_GET_CELL(id_cell, m_Grid);
      if      (type_cell == VTK_TRIANGLE)   pick_txt += "tri";
      else if (type_cell == VTK_QUAD)       pick_txt += "qua";
      else if (type_cell == VTK_POLYGON)    pick_txt += "plg";
      else if (type_cell == VTK_TETRA)      pick_txt += "tet";
      else if (type_cell == VTK_PYRAMID)    pick_txt += "pyr";
      else if (type_cell == VTK_WEDGE)      pick_txt += "pri";
      else if (type_cell == VTK_HEXAHEDRON) pick_txt += "hex";
      else if (type_cell == VTK_POLYHEDRON) pick_txt += "pol";
      pick_txt += " [";
      for (int i_pts = 0; i_pts < num_pts; ++i_pts) {
        QString num;
        num.setNum(pts[i_pts]);
        pick_txt += num;
        if (i_pts < num_pts-1) {
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
      QString tmp;
      EG_VTKDCN(vtkDoubleArray, characteristic_length_desired, m_Grid, "node_meshdensity_desired");
      tmp.setNum(characteristic_length_desired->GetValue(id_node));
      pick_txt += " wanted density=" + tmp;
      EG_VTKDCN(vtkIntArray, node_specified_density, m_Grid, "node_specified_density");
      tmp.setNum(node_specified_density->GetValue(id_node));
      pick_txt += " node_specified_density=" + tmp;
      EG_VTKDCN(vtkCharArray_t, node_type, m_Grid, "node_type");
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
      if ((ct == VTK_TRIANGLE) || (ct == VTK_QUAD) || (ct == VTK_POLYGON)) {
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
      EG_GET_CELL(id_cell, m_Grid);
      vec3_t Center(0,0,0);
      for (int p = 0; p < num_pts; ++p) {
        vec3_t M;
        m_Grid->GetPoint(pts[p],M.data());
        Center+=M.data();
      }
      vec3_t OffSet = m_ReferenceSize*triNormal(m_Grid, pts[0], pts[1], pts[2]).normalise();
      Center = 1.0/(double)num_pts*Center+OffSet;
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
        picked_cell = m_CellPicker->GetCellId();
        if (picked_cell >= 0) {
          picked_cell = cell_index->GetValue(picked_cell);
        }
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
        picked_point = m_PointPicker->GetPointId();
        if (picked_point >= 0) {
          picked_point = node_index->GetValue(picked_point);
        }
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
    m_Grid->GetCells()->ReverseCellAtId(cellId);
    // EG_GET_CELL(cellId, m_Grid);
    // QVector<vtkIdType> nodes(num_pts);
    // for (vtkIdType j = 0; j < num_pts; ++j) nodes[j]          = pts[j];
    // for (vtkIdType j = 0; j < num_pts; ++j) pts[num_pts - j - 1] = nodes[j];
  }
  updateActors();
  m_Grid->Modified();// to make sure VTK notices the changes and changes the cell colors
  //m_Renderer->GetRenderWindow()->Render();
}

void GuiMainWindow::checkSurfaceOrientation()
{
  CorrectSurfaceOrientation corr_surf;
  corr_surf();
  updateActors();
  m_Grid->Modified();// to make sure VTK notices the changes and changes the cell colors
  //m_Renderer->GetRenderWindow()->Render();
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

void GuiMainWindow::viewRight()
{
  bool use_blender;
  getSet("General","use Blender definition for front, top, etc.", true, use_blender);
  getRenderer()->ResetCamera();
  double x[3];
  getRenderer()->GetActiveCamera()->GetFocalPoint(x);
  if (use_blender) {
    x[0] += 1;
    getRenderer()->GetActiveCamera()->SetPosition(x);
    getRenderer()->GetActiveCamera()->ComputeViewPlaneNormal();
    getRenderer()->GetActiveCamera()->SetViewUp(0,0,1);
  } else {
    x[0] += 1;
    getRenderer()->GetActiveCamera()->SetPosition(x);
    getRenderer()->GetActiveCamera()->ComputeViewPlaneNormal();
    getRenderer()->GetActiveCamera()->SetViewUp(0,1,0);
  }
  getRenderer()->ResetCamera();
  getRenderWindow()->Render();
}

void GuiMainWindow::viewLeft()
{
  bool use_blender;
  getSet("General","use Blender definition for front, top, etc.", true, use_blender);
  getRenderer()->ResetCamera();
  double x[3];
  getRenderer()->GetActiveCamera()->GetFocalPoint(x);
  if (use_blender) {
    x[0] -= 1;
    getRenderer()->GetActiveCamera()->SetPosition(x);
    getRenderer()->GetActiveCamera()->ComputeViewPlaneNormal();
    getRenderer()->GetActiveCamera()->SetViewUp(0,0,1);
  } else {
    x[0] -= 1;
    getRenderer()->GetActiveCamera()->SetPosition(x);
    getRenderer()->GetActiveCamera()->ComputeViewPlaneNormal();
    getRenderer()->GetActiveCamera()->SetViewUp(0,1,0);
  }
  getRenderer()->ResetCamera();
  getRenderWindow()->Render();
}

void GuiMainWindow::viewTop()
{
  bool use_blender;
  getSet("General","use Blender definition for front, top, etc.", true, use_blender);
  getRenderer()->ResetCamera();
  double x[3];
  getRenderer()->GetActiveCamera()->GetFocalPoint(x);
  if (use_blender) {
    x[2] += 1;
    getRenderer()->GetActiveCamera()->SetPosition(x);
    getRenderer()->GetActiveCamera()->ComputeViewPlaneNormal();
    getRenderer()->GetActiveCamera()->SetViewUp(0,1,0);
  } else {
    x[1] += 1;
    getRenderer()->GetActiveCamera()->SetPosition(x);
    getRenderer()->GetActiveCamera()->ComputeViewPlaneNormal();
    getRenderer()->GetActiveCamera()->SetViewUp(0,0,-1);
  }
  getRenderer()->ResetCamera();
  getRenderWindow()->Render();
}

void GuiMainWindow::viewBottom()
{
  bool use_blender;
  getSet("General","use Blender definition for front, top, etc.", true, use_blender);
  getRenderer()->ResetCamera();
  double x[3];
  getRenderer()->GetActiveCamera()->GetFocalPoint(x);
  if (use_blender) {
    x[2] -= 1;
    getRenderer()->GetActiveCamera()->SetPosition(x);
    getRenderer()->GetActiveCamera()->ComputeViewPlaneNormal();
    getRenderer()->GetActiveCamera()->SetViewUp(0,-1,0);
  } else {
    x[1] -= 1;
    getRenderer()->GetActiveCamera()->SetPosition(x);
    getRenderer()->GetActiveCamera()->ComputeViewPlaneNormal();
    getRenderer()->GetActiveCamera()->SetViewUp(0,0,-1);
  }
  getRenderer()->ResetCamera();
  getRenderWindow()->Render();
}

void GuiMainWindow::viewFront()
{
  bool use_blender;
  getSet("General","use Blender definition for front, top, etc.", true, use_blender);
  getRenderer()->ResetCamera();
  double x[3];
  getRenderer()->GetActiveCamera()->GetFocalPoint(x);
  if (use_blender) {
    x[1] -= 1;
    getRenderer()->GetActiveCamera()->SetPosition(x);
    getRenderer()->GetActiveCamera()->ComputeViewPlaneNormal();
    getRenderer()->GetActiveCamera()->SetViewUp(0,0,1);
  } else {
    x[2] += 1;
    getRenderer()->GetActiveCamera()->SetPosition(x);
    getRenderer()->GetActiveCamera()->ComputeViewPlaneNormal();
    getRenderer()->GetActiveCamera()->SetViewUp(0,1,0);
  }
  getRenderer()->ResetCamera();
  getRenderWindow()->Render();
}

void GuiMainWindow::viewBack()
{
  bool use_blender;
  getSet("General","use Blender definition for front, top, etc.", true, use_blender);
  getRenderer()->ResetCamera();
  double x[3];
  getRenderer()->GetActiveCamera()->GetFocalPoint(x);
  if (use_blender) {
    x[1] += 1;
    getRenderer()->GetActiveCamera()->SetPosition(x);
    getRenderer()->GetActiveCamera()->ComputeViewPlaneNormal();
    getRenderer()->GetActiveCamera()->SetViewUp(0,0,1);
  } else {
    x[2] -= 1;
    getRenderer()->GetActiveCamera()->SetPosition(x);
    getRenderer()->GetActiveCamera()->ComputeViewPlaneNormal();
    getRenderer()->GetActiveCamera()->SetViewUp(0,1,0);
  }
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
      //TriSurfaceProjection tmp03;
      SurfaceMesher tmp04;
      UpdateDesiredMeshDensity tmp05;
      InsertPoints tmp06;
      RemovePoints tmp07;
      LaplaceSmoother tmp08;
      SwapTriangles tmp09;
      SolverTools tmp10;
      BlenderReader tmp11;
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
  //Load the HTML code snippet with the list of contributions
  QFileInfo fileinfo;
  fileinfo.setFile(":/contributions.htm");
  QFile file(fileinfo.filePath());
  if (!file.exists()) {
    qDebug() << "ERROR: " << fileinfo.filePath() << " not found.";
    EG_BUG;
  }
  if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) {
    qDebug() << "ERROR:  Failed to open file " << fileinfo.filePath();
    EG_BUG;
  }
  QTextStream text_stream(&file);
  QString contributionsIncluded = text_stream.readAll();
  file.close();


  //Do the About box
  QMessageBox box(this);

  QString title="ENGRID";
  QString version = QString("version=") + ENGRID_VERSION + "<br/>branch=" + GIT_BRANCH + "<br/>commit=" + GIT_SHA1;
  version += "<br/>built on ";
  version += QString(__DATE__);
  version += " at ";
  version += QString(__TIME__);

  QString address = tr("ENGRID is being developed and maintained by:<br/>"
                       "enGits GmbH<br/>"
                       "Langenbachstrasse 3<br/>"
                       "79674 Todtnau<br/>"
                       "Germany<br/>");

  QString mainurl="<a href=\"http://engits.eu/\">http://engits.eu</a>";
  QString mail="<a href=\"mailto:info@engits.com\">info@engits.com</a>";
  QString gnuurl="<a href=\"http://www.gnu.org/licenses\">http://www.gnu.org/licenses</a>";
  QString license=tr("ENGRID is licenced under the GPL version 3.<br/>"
                     "(see ")+gnuurl+tr(" for details)<br/>");
  QString contributions=tr("Contributions:");

  box.setText(QString::fromLatin1("<center><img src=\":/icons/resources/icons/G.png\">"
                                  "<h3>%1</h3>"
                                  "<p>%2</p>"
                                  "<p>%3</p>"
                                  "<p>Homepage: %4</p>"
                                  "<p>E-mail: %5</p>"
                                  "<p>%6</p></center>"
                                  "<p>%7</p><blockquote>%8</blockquote>")
              .arg(title).arg(version).arg(address).arg(mainurl).arg(mail).arg(license)
              .arg(contributions).arg(contributionsIncluded));
  box.setWindowTitle(tr("about ENGRID"));
  box.setIcon(QMessageBox::NoIcon);
  box.exec();

}

///\todo Why not use bcs = m_AllBoundaryCodes; ?
void GuiMainWindow::getAllBoundaryCodes(QVector<int> &bcs)
{
  bcs.resize(m_AllBoundaryCodes.size());
  qCopy(m_AllBoundaryCodes.begin(), m_AllBoundaryCodes.end(), bcs.begin());
  qSort(bcs);
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
    QVector<int> bcs;
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

void GuiMainWindow::storeCadInterfaces(bool nosave)
{
  try {
    resetCadInterfaces();
    CgalTriCadInterface *cad = new CgalTriCadInterface(m_Grid);
    setUniversalCadInterface(cad);
    if (!nosave) {
      save();
      saveGrid(m_Grid, m_CurrentFilename + ".geo");
    }

  } catch (Error E) {
    E.display();
  }
}

void GuiMainWindow::setUniversalCadInterface(CadInterface *cad_interface)
{
  m_UniCadInterface = cad_interface;
  cad_interface->setForegroundGrid(m_Grid);
}

void GuiMainWindow::resetCadInterfaces()
{
  delete m_UniCadInterface;
  m_UniCadInterface = NULL;
  foreach (CadInterface* cad_interface, m_CadInterfaces) {
    delete cad_interface;
  }
  m_CadInterfaces.clear();
}

CadInterface *GuiMainWindow::getCadInterface(int bc, bool allow_null)
{
  QString bc_txt;
  bc_txt.setNum(bc);
  if (!m_CadInterfaces.contains(bc)) {
    bc = 0;
  }
  if (!m_CadInterfaces.contains(bc)) {
    if (m_UniCadInterface) {
      return m_UniCadInterface;
    }
    if (allow_null) {
      return NULL;
    }
    EG_ERR_RETURN("No surface projection found for boundary code " + bc_txt);
  }
  return m_CadInterfaces[bc];
}

bool GuiMainWindow::checkCadInterfaces()
{
  bool ok = true;
  if (!m_UniCadInterface) {
    foreach (int bc, m_AllBoundaryCodes) {
      if (!m_CadInterfaces.contains(bc)) {
        ok = false;
        break;
      }
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
#if QT_VERSION < 0x040500
  pts[0] = QInputDialog::getInteger(this, tr("id_node1"),tr("id_node1:"), 0, 0, m_Grid->GetNumberOfPoints(), 1, &ok1);
  pts[1] = QInputDialog::getInteger(this, tr("id_node2"),tr("id_node2:"), 0, 0, m_Grid->GetNumberOfPoints(), 1, &ok2);
  pts[2] = QInputDialog::getInteger(this, tr("id_node3"),tr("id_node3:"), 0, 0, m_Grid->GetNumberOfPoints(), 1, &ok3);
  vtkIdType id_cell = QInputDialog::getInteger(this, tr("copy cell data from id_cell"),tr("copy cell data from id_cell:"), 0, 0, m_Grid->GetNumberOfCells(), 1, &ok4);
#else
  pts[0] = QInputDialog::getInt(this, tr("id_node1"),tr("id_node1:"), 0, 0, m_Grid->GetNumberOfPoints(), 1, &ok1);
  pts[1] = QInputDialog::getInt(this, tr("id_node2"),tr("id_node2:"), 0, 0, m_Grid->GetNumberOfPoints(), 1, &ok2);
  pts[2] = QInputDialog::getInt(this, tr("id_node3"),tr("id_node3:"), 0, 0, m_Grid->GetNumberOfPoints(), 1, &ok3);
  vtkIdType id_cell = QInputDialog::getInt(this, tr("copy cell data from id_cell"),tr("copy cell data from id_cell:"), 0, 0, m_Grid->GetNumberOfCells(), 1, &ok4);
#endif
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
#if QT_VERSION < 0x040500
  vtkIdType id_node1 = QInputDialog::getInteger(this, tr("id_node1"),tr("id_node1:"), 0, 0, m_Grid->GetNumberOfPoints(), 1, &ok1);
  vtkIdType id_node2 = QInputDialog::getInteger(this, tr("id_node2"),tr("id_node2:"), 0, 0, m_Grid->GetNumberOfPoints(), 1, &ok2);
#else
  vtkIdType id_node1 = QInputDialog::getInt(this, tr("id_node1"),tr("id_node1:"), 0, 0, m_Grid->GetNumberOfPoints(), 1, &ok1);
  vtkIdType id_node2 = QInputDialog::getInt(this, tr("id_node2"),tr("id_node2:"), 0, 0, m_Grid->GetNumberOfPoints(), 1, &ok2);
#endif
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
      EG_GET_CELL(id_cell, m_Grid);
      for (int i = 0; i < num_pts; ++i) {
        ptIds->SetId(i, old2new_nodes[pts[i]]);
      }
      vtkIdType id_new_cell = new_grid->InsertNextCell(type_cell, ptIds);
      copyCellData(m_Grid, id_cell, new_grid, id_new_cell);
    }

    makeCopy(new_grid, m_Grid);
    m_Grid->Modified();
    qDebug()<<"The fusion is complete.";
  }

}

void GuiMainWindow::onEsc()
{
  setPickMode(true, true);
  pickCell(-1);
  m_CellPicker->Pick(-1e99,-1e99,0,m_Renderer);
  updateActors(true);
  updateStatusBar();
}

void GuiMainWindow::resetProgress(QString info_text, int p_max)
{
  m_StatusInfoLabel->setText(info_text);
  m_StatusProgressBar->setMaximum(p_max);
  m_StatusProgressBar->setValue(0);
  QApplication::processEvents();
}

void GuiMainWindow::setProgress(int p)
{
  m_StatusProgressBar->setValue(p);
  for (int i = 0; i < 3; ++i) {
    QApplication::processEvents();
  }
}

void GuiMainWindow::lock()
{
  m_Mutex.lock();
}

void GuiMainWindow::unlock()
{
  m_Mutex.unlock();
}

bool GuiMainWindow::tryLock()
{
  return m_Mutex.tryLock();
}

BoundaryCondition GuiMainWindow::getBC(BoundaryCondition BC)
{
  int last_code = 0;
  foreach (BoundaryCondition bc, m_bcmap.values()) {
    if (bc == BC) {
      return bc;
    }
    last_code = max(last_code, bc.getCode());
  }
  BC.setCode(last_code + 1);
  setBC(BC.getCode(), BC);
  return BC;
}
