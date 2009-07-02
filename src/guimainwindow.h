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
#include <QDomDocument>

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
#include "surfaceprojection.h"

#include "std_includes.h"
#include "guitransform.h"

/**
 * This is the main GUI class of enGrid.
 */
class GuiMainWindow : public QMainWindow, public EgVtkObject
{
  
  Q_OBJECT;
  
private: // attributes
  
  QDomDocument         m_XmlDoc;        ///< XML document describing the complete case

  Ui::GuiMainWindow    ui;              ///< The user interface definition -- created by QtDesigner.
  vtkUnstructuredGrid *grid;            ///< The current state of the grid that is being generated.

  vtkRenderer *renderer; ///< The VTK renderer object, used for visualising the grid
  
  vtkActor* m_SurfaceActor;
  vtkActor* m_SurfaceWireActor;
  vtkActor* m_TetraActor;
  vtkActor* m_WedgeActor;
  vtkActor* m_PyramidActor;
  vtkActor* m_HexaActor;
  vtkActor* m_VolumeWireActor;

  vtkProperty*       backface_property;
  vtkLookupTable*    lut;
  vtkScalarBarActor* m_LegendActor;
  
  vtkPolyDataMapper* m_SurfaceMapper;
  vtkPolyDataMapper* m_SurfaceWireMapper;
  vtkPolyDataMapper* m_TetraMapper;
  vtkPolyDataMapper* m_PyramidMapper;
  vtkPolyDataMapper* m_WedgeMapper;
  vtkPolyDataMapper* m_HexaMapper;
  vtkPolyDataMapper* m_VolumeWireMapper;
  
  vtkEgExtractVolumeCells *m_ExtrVol;
  vtkEgExtractVolumeCells *m_ExtrTetras;
  vtkEgExtractVolumeCells *m_ExtrPyramids;
  vtkEgExtractVolumeCells *m_ExtrWedges;
  vtkEgExtractVolumeCells *m_ExtrHexes;
  
  vtkGeometryFilter *volume_geometry;
  vtkGeometryFilter *tetra_geometry;
  vtkGeometryFilter *pyramid_geometry;
  vtkGeometryFilter *wedge_geometry;
  vtkGeometryFilter *hexa_geometry;

  vtkIdType m_PickedPoint;      ///< Picked point
  vtkIdType m_PickedCell;       ///< Picked cell
  bool      m_UseVTKInteractor; ///< Boolean value specifying whether the VTK Interactor should be used or not

  static QMutex    mutex;

  vtkGeometryFilter* surface_filter;  ///< VTK filter to extract the surface of the current grid.
  double             m_ReferenceSize; ///< Size to use for picker objects and annotations

  vector <vtkTextActor*>      m_NodeText;               ///< 2D Text actor to display node IDs
  vector <vtkTextActor*>      m_CellText;               ///< 2D Text actor to display cell IDs
  vector <vtkVectorText*>     m_NodeTextVectorText;     ///< 3D Text actor to display node IDs
  vector <vtkPolyDataMapper*> m_NodeTextPolyDataMapper;
  vector <vtkFollower*>       m_NodeTextFollower;
  vector <vtkVectorText*>     m_CellTextVectorText;     ///< 3D Text actor to display cell IDs
  vector <vtkPolyDataMapper*> m_CellTextPolyDataMapper;
  vector <vtkFollower*>       m_CellTextFollower;

  vtkPolyDataMapper*        m_PickMapper;   ///< VTK mapper to map pick marker
  vtkActor*                 m_PickActor;    ///< VTK actor to display pick marker
  vtkSphereSource*          m_PickSphere;   ///< sphere to mark picked cell/points
  vtkCubeAxesActor2D*       m_Axes;           ///< VTK actor to display the coordinate system
  vtkEgBoundaryCodesFilter* m_BCodesFilter; ///< VTK filter to extract boundary elements with certain codes
  vtkCellPicker*            m_CellPicker;   ///< VTK CellPicker to pick cells for various user interactions
  vtkPointPicker*           m_PointPicker;  ///< VTK PointPicker to pick points for various user interactions

  QString      m_CurrentFilename;      ///< The current file name of the grid.
  int          current_operation;      ///< The current operation number. (used for undo/redo)
  int          last_operation;         ///< The last operation number. (used for undo/redo)
  QString      m_LogDir;               ///< the log directory
  QStatusBar*  status_bar;             ///< Status bar of the main window and application
  QLabel*      status_label;           ///< Label for the information in the status bar
  QSet<int>    m_DisplayBoundaryCodes; ///< A QList with all active boundary codes.
  QSet<int>    m_AllBoundaryCodes;     ///< A QList with all boundary codes.
  bool         busy;                   ///< flag to indicate that enGrid is busy with an operation
  QString      log_file_name;          ///< log file to collect program output for display in the output window
  long int     N_chars;                ///< number of lines that have been read from the log file
  FILE*        system_stdout;
  QTimer       garbage_timer;
  QTimer       log_timer;
  QDockWidget* dock_widget;

  QMap<int,BoundaryCondition>    bcmap;      ///< mapping between numerical and symbolic boundary codes
  QMap<QString,VolumeDefinition> volmap;     ///< all volume definitions

  QMap<int,SurfaceProjection*>   m_SurfProj; ///< all surface projectors for surface meshing
  
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
  
  void        setupVtk();
  void        addVtkTypeInfo(); ///< Add VTK type information to the grid (useful for visualisation with ParaView).
  static void pickCellCallBack(vtkObject *caller, unsigned long int eid, void *clientdata, void *calldata);
  static void pickPointCallBack(vtkObject *caller, unsigned long int eid, void *clientdata, void *calldata);
  void        updateSurfaceActors(bool forced);
  void        updateVolumeActors(bool forced);

private slots:

  void setClipX(const QString &txt);
  void setClipY(const QString &txt);
  void setClipZ(const QString &txt);
  void setClipNX(const QString &txt);
  void setClipNY(const QString &txt);
  void setClipNZ(const QString &txt);

  void openBC();
  void saveBC();
  void openGrid(QString file_name);
  void saveGrid(QString file_name);



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
  vtkIdType getPickedCell();

  /**
   * Get the currently picked point.
   * @return the picked point ID or -1 if no point has been picked
   */
  vtkIdType getPickedPoint();
  
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
  vtkPointPicker* getPointPicker() { return (m_PointPicker);}
  vtkSphereSource* getPickSphere() { return (m_PickSphere);}
  bool pickPoint(vtkIdType id_point);
  bool pickCell(vtkIdType id_cell);
  
  QString getFilename() { return(m_CurrentFilename); }

  SurfaceProjection* getSurfProj(int bc) { return m_SurfProj[bc]; }
  
public slots:

  void setUseVTKInteractor(int a_UseVTKInteractor);
  void setPickMode(bool a_UseVTKInteractor,bool a_CellPickerMode);
  
  void exit();                           ///< Exit the application
  void importSTL();                      ///< Import an STL file (ASCII or binary)
  void importGmsh1Ascii();               ///< Import a Gmsh grid from an ASCII file -- using version 1.0 of the Gmsh file format
  void exportGmsh1Ascii();               ///< Export a grid from to an ASCII Gmsh file -- using version 1.0 of the Gmsh file format
  void importGmsh2Ascii();               ///< Import a Gmsh grid from an ASCII file -- using version 2.0 of the Gmsh file format
  void exportGmsh2Ascii();               ///< Export a grid from to an ASCII Gmsh file -- using version 2.0 of the Gmsh file format
  void exportNeutral();                  ///< Export a grid to neutral format for NETGEN
  void updateActors(bool force = false); ///< Update the VTK output
  void forceUpdateActors();              ///< Force an update of the VTK output
  void scaleToData();                    ///< Scale to data
  void zoomAll();                        ///< Move the camera in order to show everything on the screen
  void zoomOnPickedObject();
  void deselectAll();
  void printGrid() {cout<<"PrintGrid() called!"<<endl; cout_grid(cout,grid,true,true,true,true);}
  void info();
  
  void undo();
  void redo();
  
  void resetOperationCounter();
  
  ///@@@  TODO: Simplify available save/load functions

  void saveXml();                        ///< Save the case in an XML file
  void open();                           ///< Open an existing case
  void save();                           ///< Save the current case
  void saveAs();                         ///< Save the current case -- using a different file name

  void quickSave(QString a_filename);    ///< Save the current grid as a_filename
  void quickLoad(QString a_filename);    ///< Load the current grid from a_filename
  int  quickSave();                      ///< Save the current grid as a_filename_a_operation
  void quickLoad(int a_operation);       ///< Load a_filename_a_operation
  void updateStatusBar();                ///< Update the status bar
  void selectBoundaryCodes();            ///< Select the boundary codes to be displayed/hidden
  void updateBoundaryCodes(bool all_on); ///< Update the boundary code book keeping (e.g. after reading a mesh).
  void normalExtrusion();                ///< Normal extrusion of boundary elements (no validity check).
  void setAxesVisibility();              ///< Toggle the visibility of the axes annotation.
  void setViewingMode();                 ///< Toggle orthogonal viewing mode.
  void viewNodeIDs();                    ///< Toggle node ID viewing mode.
  void viewCellIDs();                    ///< Toggle cell ID viewing mode.
  void changeSurfaceOrientation();       ///< Change the orientation of all surface elements
  void checkSurfaceOrientation();        ///< Check and, if required, change the orientation of all surface elements
  void improveAspectRatio();             ///< Eliminate edges in order to improve the aspect ratio of the cells
  void exportAsciiStl();                 ///< Write surface elements to an ASCII STL file.
  void exportBinaryStl();                ///< Write surface elements to a binary STL file.
  void editBoundaryConditions();         ///< Edit boundary conditions (names and types)
  void configure();                      ///< Edit settings
  void about();                          ///< Display an about message
  void markOutputLine();                 ///< Mark the current position in the output window

  QString getXmlSection(QString name);                   ///< Get a section from the XML case description
  void    setXmlSection(QString name, QString contents); ///< Set a section of the XML case description

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

  void storeSurfaceProjection();
    
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
