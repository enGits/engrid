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
#ifndef mainwindow_H
#define mainwindow_H

class GuiMainWindow;

#include <QMainWindow>
#include <QSettings>
#include <QLabel>
#include <QSet>
#include <QFileSystemWatcher>
#include <QMutex>
#include <QTimer>
#include <QDockWidget>

#include <vtkUnstructuredGrid.h>
#include <vtkActor.h>
#include <vtkPolyDataMapper.h>
#include <vtkGeometryFilter.h>
#include <vtkCubeAxesActor2D.h>
#include <vtkCellPicker.h>
#include <vtkPointPicker.h>
#include <vtkSphereSource.h>
#include <vtkTextActor.h>
#include <vtkVectorText.h>
#include <vtkFollower.h>
#include <vtkScalarBarActor.h>
#include <vtkLookupTable.h>

#include "ui_guimainwindow.h"
#include "vtkEgBoundaryCodesFilter.h"
#include "vtkEgExtractVolumeCells.h"
#include "egvtkobject.h"
#include "boundarycondition.h"
#include "volumedefinition.h"
#include "checksurfaceintegrity.h"

#include "std_includes.h"
#include "guitransform.h"

/**
 * This is the main GUI class of enGrid.
 */
class GuiMainWindow : public QMainWindow, public EgVtkObject
{
  
  Q_OBJECT;
  
private: // attributes
  
  /** The user interface definition -- created by QtDesigner. */
  Ui::GuiMainWindow ui;
  

  vtkUnstructuredGrid *grid;           /// The current state of the grid that is being generated.
  unsigned long int    m_GridMTime;    /// last modification of the grid
  bool                 m_TetraExtr;    /// flag if tetras have been extracted yet
  bool                 m_PyraExtr;     /// flag if pyramids have been extracted yet
  bool                 m_PrismExtr;    /// flag if prisms have been extracted yet
  bool                 m_HexaExtr;     /// flag if hexes have been extracted yet
  bool                 m_VolExtr;      /// flag if any vol. cells have been extracted yet
  bool                 m_UsedClipping; /// flag if clipping had been used for last update
  QString              m_ClippingText; /// last clipping input from UI

  vtkPolyData *boundary_pd;
  vtkPolyData *tetras_pd;
  vtkPolyData *wedges_pd;
  vtkPolyData *pyras_pd;
  vtkPolyData *hexas_pd;
  vtkPolyData *volume_pd;
  
  /** The VTK renderer object, used for visualising the grid */
  vtkRenderer *renderer;
  
  vtkActor *surface_actor;
  vtkActor *surface_wire_actor;
  vtkActor *tetra_actor;
  vtkActor *wedge_actor;
  vtkActor *pyramid_actor;
  vtkActor *hexa_actor;
  vtkActor *volume_wire_actor;
  vtkProperty *backface_property;
  
  vtkLookupTable *lut;
  vtkScalarBarActor *iamlegend_actor;
  
  vtkPolyDataMapper *field_mapper;
  vtkPolyDataMapper *surface_mapper;
  vtkPolyDataMapper *surface_wire_mapper;
  vtkPolyDataMapper *tetra_mapper;
  vtkPolyDataMapper *pyramid_mapper;
  vtkPolyDataMapper *wedge_mapper;
  vtkPolyDataMapper *hexa_mapper;
  vtkPolyDataMapper *volume_wire_mapper;
  
  vtkEgExtractVolumeCells *extr_vol;
  vtkEgExtractVolumeCells *extr_tetras;
  vtkEgExtractVolumeCells *extr_pyramids;
  vtkEgExtractVolumeCells *extr_wedges;
  vtkEgExtractVolumeCells *extr_hexas;
  
  vtkGeometryFilter *volume_geometry;
  vtkGeometryFilter *tetra_geometry;
  vtkGeometryFilter *pyramid_geometry;
  vtkGeometryFilter *wedge_geometry;
  vtkGeometryFilter *hexa_geometry;
  
  /** VTK filter to extract the surface of the current grid. */
  vtkGeometryFilter *surface_filter;
  
  /** sphere to mark picked cell/points */
  vtkSphereSource *pick_sphere;
  
  /** Size to use for picker objects and annotations */
  double ReferenceSize;
  
  /** 2D Text actor to display node IDs */
  vector <vtkTextActor*> NodeText;
  
  /** 2D Text actor to display cell IDs */
  vector <vtkTextActor*> CellText;
  
  /** 3D Text actor to display node IDs */
  vector <vtkVectorText*> NodeText_VectorText;
  vector <vtkPolyDataMapper*> NodeText_PolyDataMapper;
  vector <vtkFollower*> NodeText_Follower;
  
  /** 3D Text actor to display cell IDs */
  vector <vtkVectorText*> CellText_VectorText;
  vector <vtkPolyDataMapper*> CellText_PolyDataMapper;
  vector <vtkFollower*> CellText_Follower;
  
  /** Picked point */
  static vtkIdType PickedPoint;
  
  /** Picked cell */
  static vtkIdType PickedCell;
  
  /** Boolean value specifying wether the VTK Interactor should be used or not to determine picked points/cells */
  static bool m_UseVTKInteractor;
  
/*  vtkTextActor* textActor[3];
  
  vtkVectorText* atext[3];
  vtkPolyDataMapper* textMapper[3];
  vtkFollower* textActor2[3];*/
  
  /** VTK mapper to map pick marker */
  vtkPolyDataMapper *pick_mapper;
  
  /** VTK actor to display pick marker */
  vtkActor *pick_actor;
  
  /** VTK actor to display the coordinate system */
  vtkCubeAxesActor2D *axes;
  
  /** The current file name of the grid. */
  QString current_filename;
  
  /** The current operation number. (used for undo/redo) */
  int current_operation;
  
  /** The last operation number. (used for undo/redo) */
  int last_operation;
  
  /** the log directory */
  QString m_LogDir;
  
  /** Status bar of the main window and application */
  QStatusBar *status_bar;
  
  /** Label for the information in the status bar */
  QLabel *status_label;
  
  /** A QList with all active boundary codes. */
  QSet<int> display_boundary_codes;
  
    /** A QList with all boundary codes. */
  QSet<int> all_boundary_codes;
  
  /** VTK filter to extract boundary elements with certain codes */
  vtkEgBoundaryCodesFilter *bcodes_filter;
  
  /** VTK CellPicker to pick cells for various user interactions */
  vtkCellPicker *CellPicker;

  /** VTK PointPicker to pick points for various user interactions */
  vtkPointPicker *PointPicker;
  
  /** flag to indicate that enGrid is busy with an operation */
  bool busy;
  
  /** log file to collect program output for display in the output window */
  QString log_file_name;
  
  /** number of lines that have been read from the log file */
  long int N_chars;
  
  FILE *system_stdout;
  static QMutex mutex;
  QTimer garbage_timer;
  QTimer log_timer;
  QDockWidget *dock_widget;
  
  /** mapping between numerical and symbolic boundary codes */
  QMap<int,BoundaryCondition> bcmap;

  /** all volume definitions */
  QMap<QString,VolumeDefinition> volmap;


  
private: // static attributes

  /**
   * Platform independant access to application settings.
   * For a UNIX system the user preferences will be stored in the file
   * folder ".config/enGits/enGrid.conf" in the user's home directory;
   * on Windows preferences will be stored in the registry.
   */
  static QSettings qset;
  
  /**
   * The current working directory of enGrid
   */
  static QString cwd;
  
  /** a static this pointer (somewhat ugly, but there is only one MainWindow) */
  static GuiMainWindow *THIS;
  
private: // methods
  
  /** Add VTK type information to the grid (useful for visualisation with ParaView). */
  void addVtkTypeInfo();
  
  /** callback for cell picking */
  static void pickCellCallBack(vtkObject *caller, unsigned long int eid, void *clientdata, void *calldata);

  /** callback for point picking */
  static void pickPointCallBack(vtkObject *caller, unsigned long int eid, void *clientdata, void *calldata);
  
public: // methods
  
  /**
   * The constructor connects the menu and toolbar actions and 
   * the VTK basics(i.e. renderer, actor, ...) will be set up.
   * Furthermore preferences will be read from qset.
   */
  GuiMainWindow();
  
  /**
   * Preferences will be written back.
   */
  virtual ~GuiMainWindow();
  
  /**
   * Get the VTK render window
   * @return the VTK render window
   */
  vtkRenderWindow* getRenderWindow();
  
  /**
   * Get the VTK renderer
   * @return the VTK renderer
   */
  vtkRenderer* getRenderer();
  
  /** 
   * Get the Qt-VTK interactor
   * @return the Qt-VTK interactor
   */
  QVTKInteractor* getInteractor();
  
  /**
   * Get a pointer to the current grid object
   * @return a pointer to the current vtkUnstructuredGrid object
   */
  vtkUnstructuredGrid* getGrid() { return grid; }
  
  void setBusy() { busy = true; updateStatusBar(); }
  void setIdle() { busy = false; updateStatusBar(); }
  
  /// Returns log directory
  QString getLogDir() { return m_LogDir; }
  
  /// Returns the path to the currently loaded file
  QString getFilePath();
  
public: // static methods
  
  /**
   * Get the current working directory.
   * @return the current working directory
   */
  static QString getCwd();
  
  /**
   * Set the current working directory
   * @param cwd the current working directory
   */
  static void setCwd(QString dir);
  
  /**
   * Get the currently picked cell.
   * @return the picked cell ID or -1 if no cell has been picked
   */
  static vtkIdType getPickedCell();

  /**
   * Get the currently picked point.
   * @return the picked point ID or -1 if no point has been picked
   */
  static vtkIdType getPickedPoint();
  
  /**
   * Access to the QSettings object/
   */
  static QSettings* settings() { return &qset; }
  
  BoundaryCondition getBC(int bc) { return bcmap[bc]; }
  VolumeDefinition  getVol(QString volname) { return volmap[volname]; }
  QList<VolumeDefinition> getAllVols();
  void setAllVols(QList<VolumeDefinition> vols);
  void createDefaultVol();
  
  static GuiMainWindow* pointer() { return THIS; }
  static void lock() { mutex.lock(); }
  static void unlock() { mutex.unlock(); }
  static bool tryLock() { return mutex.tryLock(); }
  void getAllBoundaryCodes(QSet<int> &bcs);
  void getDisplayBoundaryCodes(QSet<int> &bcs);
  vtkPointPicker* getPointPicker(){return(PointPicker);}
  vtkSphereSource* getPickSphere(){return(pick_sphere);}
  bool pickPoint(vtkIdType Point);
  bool pickCell(vtkIdType cellId);
  
  QString getFilename() { return(current_filename); }
  
public slots:
  void setUseVTKInteractor(int a_UseVTKInteractor);
  void setPickMode(bool a_UseVTKInteractor,bool a_CellPickerMode);
  
  /** Exit the application */
  void exit();
  
  /** Import an STL file (ASCII or binary) */
  void importSTL();
  
  /** Import a Gmsh grid from an ASCII file -- using version 1.0 of the Gmsh file format */
  void importGmsh1Ascii();
  
  /** Export a grid from to an ASCII Gmsh file -- using version 1.0 of the Gmsh file format */
  void exportGmsh1Ascii();
  
  /** Import a Gmsh grid from an ASCII file -- using version 2.0 of the Gmsh file format */
  void importGmsh2Ascii();
  
  /** Export a grid from to an ASCII Gmsh file -- using version 2.0 of the Gmsh file format */
  void exportGmsh2Ascii();
  
  /** Export a grid to neutral format for NETGEN */
  void exportNeutral();
  
  /** Update the VTK output */
  void updateActors();
  
  /** Scale to data */
  void ScaleToData();
  
  /** Move the camera in order to show everything on the screen */
  void zoomAll();
  void ZoomOnPickedObject();
  void DeselectAll();
  void PrintGrid() {cout<<"PrintGrid() called!"<<endl; cout_grid(cout,grid,true,true,true,true);}
  void Info();
  
  void openBC();
  void saveBC();
  void openBC(QString a_file);
  void saveBC(QString a_file);
  
  void Undo();
  void Redo();
  
  void ResetOperationCounter();
  
  ///@@@  TODO: Simplify available save/load functions
  /** Open an existing grid */
  void open();
  
  /** Save the current grid */
  void save();
  
  /** Save the current grid -- using a different file name */
  void saveAs();
  
  /** Save the current grid as a_filename */
  void QuickSave(QString a_filename);
  
  /** Load the current grid from a_filename */
  void QuickLoad(QString a_filename);
  
  /** Save the current grid as a_filename_a_operation */
  int QuickSave();
  
  /** Load a_filename_a_operation */
  void QuickLoad(int a_operation);
  
  /** Update the status bar */
  void updateStatusBar();
  
  /** Select the boundary codes to be displayed/hidden */
  void selectBoundaryCodes();
  
  /** Update the boundary code book keeping (e.g. after reading a mesh). */
  void updateBoundaryCodes(bool all_on);
  
  /** Normal extrusion of boundary elements (no validity check). */
  void normalExtrusion();
  
  /** Toggle the visibility of the axes annotation. */
  void setAxesVisibility();
  
  /** Toggle orthogonal viewing mode. */
  void setViewingMode();
  
  /** Toggle node ID viewing mode. */
  void ViewNodeIDs();
  
  /** Toggle cell ID viewing mode. */
  void ViewCellIDs();
  
  /** Change the orientation of all surface elements */
  void changeSurfaceOrientation();
  
  /** Check and, if required, change the orientation of all surface elements */
  void checkSurfaceOrientation();
  
  /** Eliminate edges in order to improve the aspect ratio of the cells */
  void improveAspectRatio();
  
  /** Write surface elements to an ASCII STL file. */
  void exportAsciiStl();
  
  /** Write surface elements to a binary STL file. */
  void exportBinaryStl();
  
  /** Edit boundary conditions (names and types) */
  void editBoundaryConditions();
  
  /** Edit settings */
  void configure();
  
  /** Display an about message */
  void about();
  
  void viewXP();
  void viewXM();
  void viewYP();
  void viewYM();
  void viewZP();
  void viewZM();
  
  void appendOutput(QString txt) { ui.textEditOutput->append(txt); }
  void clearOutput() { ui.textEditOutput->clear(); }
  void updateOutput();
  void periodicUpdate();
    
  // SLOTS for all standard operations should be defined below;
  // entries should look like this:
  // void callOperationName() { EG_STDSLOT(OperationName); };
  // The actual class name in this case, however, would be GuiOperationName.
  //
  // the following line can be used as a template:
  // void call() { EG_STDSLOT(); };
  // IMPORTANT: Using EG_STDSLOT sets gui to true, while EG_STDINTERSLOT does not (default is gui = false)
  // This is important to determine whether an operation is a GUI operation or not.
  // If it's a GUI operation, it locks everything.
  // Note: In practice, EG_STDINTERSLOT locks everything, while EG_STDSLOT prevents other operations, but doesn't lock the text output or prevent minimizing the window... Why?
  
  void callSmoothSurface() { EG_STDSLOT(GuiSmoothSurface); }
  void callCreateBoundaryLayer() { EG_STDSLOT(GuiCreateBoundaryLayer); }
  void callDivideBoundaryLayer() { EG_STDSLOT(GuiDivideBoundaryLayer); }
  void callDeleteVolumeGrid() { EG_STDSLOT(DeleteVolumeGrid); }
  void callDeleteTetras() { EG_STDSLOT(DeleteTetras); }
  void callCreateVolumeMesh() { EG_STDSLOT(GuiCreateVolumeMesh); }
  void callSmoothVolumeGrid() { EG_STDSLOT(SmoothVolumeGrid); }
  void callSetBoundaryCode()  { EG_STDINTERSLOT(GuiSetBoundaryCode); }
  void callDeleteBadAspectTris() { EG_STDINTERSLOT(GuiDeleteBadAspectTris); }
  void callDeletePickedCell() { EG_STDSLOT(DeletePickedCell); }
  void callDeletePickedPoint() { EG_STDINTERSLOT(DeletePickedPoint); }
  void callBoxSelect() { EG_STDINTERSLOT(BoxSelect); }
  void callCheckSurfaceIntegrity() { EG_STDINTERSLOT(CheckSurfaceIntegrity); }
  void callPick_cell_point() { EG_STDINTERSLOT(GuiPick); }
  void callTransform() { EG_STDINTERSLOT(GuiTransform); }
  
  void callFixSTL();
  
  void callFoamReader()       { EG_STDREADERSLOT(FoamReader); }
  void callFoamWriter()       { EG_STDINTERSLOT(FoamWriter); }
  void callSimpleFoamWriter() { EG_STDINTERSLOT(SimpleFoamWriter); }
  void callCgnsWriter()       { EG_STDINTERSLOT(CgnsWriter); }
  void callVtkReader()        { EG_STDREADERSLOT(VtkReader); }
  void callPolyDataReader()   { EG_STDREADERSLOT(PolyDataReader); }
  
};

#endif
