// 
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008-2012 enGits GmbH                                     +
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

#ifdef DLL_EXPORT
   #if defined(LIBENGRID_EXPORTS) || defined(libengrid_EXPORTS)
      #define LIBENGRID_DLL   __declspec(dllexport)
   #else
      #define LIBENGRID_DLL   __declspec(dllimport)
   #endif
   #define CLASS_LIBENGRID_DLL LIBENGRID_DLL
#else
   #define LIBENGRID_DLL
   #define CLASS_LIBENGRID_DLL
#endif

#include <stdio.h>

#include <QMainWindow>
#include <QSettings>
#include <QLabel>
#include <QSet>
#include <QFileSystemWatcher>
#include <QMutex>
#include <QTimer>
#include <QDockWidget>
#include <QDomDocument>
#include <QProgressBar>

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
#include "physicalboundarycondition.h"
#include "checksurfaceintegrity.h"

#include "surfaceprojection.h"

#include "openfoamcase.h"
#include "guitransform.h"
#include "openfoamtools.h"
#include "std_includes.h"
#include "fixcadgeometry.h"
#include "xmlhandler.h"

/**
 * This is the main GUI class of enGrid.
 */
class CLASS_LIBENGRID_DLL GuiMainWindow : public QMainWindow, public EgVtkObject
{

    Q_OBJECT;

  private: // attributes

    XmlHandler* m_XmlHandler;

    Ui::GuiMainWindow    ui;     ///< The user interface definition -- created by QtDesigner.
    vtkUnstructuredGrid *m_Grid; ///< The current state of the grid that is being generated.

    vtkRenderer *m_Renderer; ///< The VTK renderer object, used for visualising the grid

    vtkActor* m_SurfaceActor;
    vtkActor* m_SurfaceWireActor;
    vtkActor* m_TetraActor;
    vtkActor* m_WedgeActor;
    vtkActor* m_PyramidActor;
    vtkActor* m_HexaActor;
    vtkActor* m_PolyhedraActor;
    vtkActor* m_VolumeWireActor;

    vtkProperty*       m_BackfaceProperty;
    vtkLookupTable*    m_LookupTable;
    vtkScalarBarActor* m_LegendActor;

    vtkPolyDataMapper* m_SurfaceMapper;
    vtkPolyDataMapper* m_SurfaceWireMapper;
    vtkPolyDataMapper* m_TetraMapper;
    vtkPolyDataMapper* m_PyramidMapper;
    vtkPolyDataMapper* m_WedgeMapper;
    vtkPolyDataMapper* m_HexaMapper;
    vtkPolyDataMapper* m_PolyhedraMapper;
    vtkPolyDataMapper* m_VolumeWireMapper;

    double m_ColTetraR, m_ColTetraG, m_ColTetraB;
    double m_ColPyraR,  m_ColPyraG,  m_ColPyraB;
    double m_ColPrismR, m_ColPrismG, m_ColPrismB;
    double m_ColHexR,   m_ColHexG,   m_ColHexB;
    double m_ColAR,     m_ColAG,     m_ColAB;
    double m_ColBR,     m_ColBG,     m_ColBB;

    vtkEgExtractVolumeCells *m_ExtrVol;
    vtkEgExtractVolumeCells *m_ExtrTetras;
    vtkEgExtractVolumeCells *m_ExtrPyramids;
    vtkEgExtractVolumeCells *m_ExtrWedges;
    vtkEgExtractVolumeCells *m_ExtrHexes;
    vtkEgExtractVolumeCells *m_ExtrPolyhedra;

    vtkGeometryFilter *m_VolumeGeometry;
    vtkGeometryFilter *m_TetraGeometry;
    vtkGeometryFilter *m_PyramidGeometry;
    vtkGeometryFilter *m_WedgeGeometry;
    vtkGeometryFilter *m_HexaGeometry;
    vtkGeometryFilter *m_PolyhedraGeometry;

    vtkIdType m_PickedPoint;      ///< Picked point
    vtkIdType m_PickedCell;       ///< Picked cell
    bool      m_UseVTKInteractor; ///< Boolean value specifying whether the VTK Interactor should be used or not

    static QMutex    m_Mutex;

    vtkGeometryFilter* m_SurfaceFilter; ///< VTK filter to extract the surface of the current grid.
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
    vtkCubeAxesActor2D*       m_Axes;         ///< VTK actor to display the coordinate system
    vtkEgBoundaryCodesFilter* m_BCodesFilter; ///< VTK filter to extract boundary elements with certain codes
    vtkCellPicker*            m_CellPicker;   ///< VTK CellPicker to pick cells for various user interactions
    vtkPointPicker*           m_PointPicker;  ///< VTK PointPicker to pick points for various user interactions
    int                       m_PickedObject; ///< 0=none, 1=node, 2=cell

    QString       m_CurrentFilename;      ///< The current file name of the grid.
    int           m_CurrentOperation;     ///< The current operation number. (used for undo/redo)
    bool          m_undo_redo_enabled;    ///< if true, undo/redo operations will be usable.
    int           m_LastOperation;        ///< The last operation number. (used for undo/redo)
    QString       m_LogDir;               ///< the log directory
    QLabel*       m_StatusLabel;          ///< Label for the information in the status bar
    QLabel*       m_StatusInfoLabel;
    QProgressBar* m_StatusProgressBar;
    QSet<int>     m_DisplayBoundaryCodes; ///< A QList with all active boundary codes.
    QSet<int>     m_AllBoundaryCodes;     ///< A QList with all boundary codes.
    bool          m_Busy;                 ///< flag to indicate that enGrid is busy with an operation
    QString       m_LogFileName;          ///< log file to collect program output for display in the output window
    long int      m_N_chars;              ///< number of lines that have been read from the log file
#if defined( __linux__ ) //for Linux
    int           m_SystemStdout;
    int           m_LogFileStdout;
    fpos_t m_SystemStdout_pos;
    fpos_t m_LogFileStdout_pos;
#elif defined( _WIN32 ) //for Windows
    //Windows always uses CON
    FILE*        m_SystemStdout;
    FILE*        m_LogFileStdout;
#else
  #error "Please define the proper way to save/recover the stdout."
#endif

    QTimer       m_GarbageTimer;
    QTimer       m_LogTimer;

    QMap<int, BoundaryCondition>    m_bcmap;       ///< mapping between numerical and symbolic boundary codes
    QMap<QString, VolumeDefinition> m_VolMap;      ///< all volume definitions
    QMap<QString, PhysicalBoundaryCondition> m_PhysicalBoundaryConditionsMap;    ///< all physical boundary conditions definitions

    QMap<int, SurfaceProjection*>   m_SurfProj;    ///< all surface projectors for surface meshing
    SurfaceProjection              *m_UniSurfProj; ///< universal surface projection for all boundary conditions

    QMap<QAction*, Operation*> m_PluginOperations;
    QAction* m_EscAction;

    int m_SolverIndex;// deprecated
    OpenFOAMTools m_OpenFoamTools;

  // recent file list support
  private:
    QMap<QString,QDateTime> m_RecentFiles;
    QMenu* recentFileMenu() { return ui.menuOpen_recent; }
    void readRecentFiles();
    void writeRecentFiles();
    void addRecentFile(QString file_name, QDateTime date);
  private slots:
    void openRecent(QAction *action);

  public:
    void resetXmlDoc();

  private: // static attributes

    /**
     * Platform independant access to application settings.
     * For a UNIX system the user preferences will be stored in the file
     * folder ".config/enGits/enGrid.conf" in the user's home directory;
     * on Windows preferences will be stored in the registry.
     */
    static QSettings m_qset;

    /**
     * The current working directory of enGrid
     */
    static QString m_cwd;

    /**
     * Is the current case unsaved?
     */
    static bool m_UnSaved;

    /** a static this pointer (somewhat ugly, but there is only one MainWindow) */
    static GuiMainWindow* THIS;

  private: // methods

    void        setupVtk();
    static void pickCallBack( vtkObject *caller, unsigned long int eid, void *clientdata, void *calldata );
    void        updateSurfaceActors( bool forced );
    void        updateVolumeActors( bool forced );

  private slots:

    void setClipX( const QString &txt );
    void setClipY( const QString &txt );
    void setClipZ( const QString &txt );
    void setClipNX( const QString &txt );
    void setClipNY( const QString &txt );
    void setClipNZ( const QString &txt );

    void openBC();
    void saveBC();

    void openPhysicalBoundaryConditions();
    void savePhysicalBoundaryConditions();

    void openGrid(QString file_name);

    void pluginCalled();
    void onEsc();

  public: // methods
    GuiMainWindow();///< Default constructor.
    GuiMainWindow(QString file_name);///< Constructor which opens a file directly.

  private:
    /**
     * This function connects the menu and toolbar actions and
     * the VTK basics(i.e. renderer, actor, ...) will be set up.
     * Furthermore preferences will be read from qset.
     */
    void setupGuiMainWindow();
    bool m_open_last;

  public:
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
    vtkUnstructuredGrid* getGrid() { return m_Grid; }

    void setBusy() { m_Busy = true; updateStatusBar(); }
    void setIdle() { m_Busy = false; updateStatusBar(); }

    /// Returns log directory
    QString getLogDir() { return m_LogDir; }

    /// Returns the path to the currently loaded file
    QString getFilePath();

    /// Returns the index of the solver to use. The index corresponds to the position in solvers.txt .
    void setSolverIndex(int x) {m_SolverIndex = x;}
    int getSolverIndex() {return m_SolverIndex;}

  public: // static methods

    /**
     * Get the current working directory.
     * @return the current working directory
     */
    static QString getCwd();

    /**
     * Set the current working directory
     * @param dir the current working directory
     */
    static void setCwd( QString dir );

    /**
     * Set m_UnSaved.
     * @param unsaved Do you want to be asked where to save when clicking on save next time?
     */
    static void setUnsaved( bool unsaved );

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

    vtkIdType getPickedObject() { return m_PickedObject; }

    /**
     * Access to the QSettings object
     */
    static QSettings* settings() { return &m_qset; }

    BoundaryCondition getBC( int bc ) { return m_bcmap[bc]; }
    VolumeDefinition  getVol( QString volname ) { return m_VolMap[volname]; }
    void clearBCs() { m_bcmap.clear(); }
    void addBC(int bc, BoundaryCondition BC) { m_bcmap[bc] = BC; }

    QList<VolumeDefinition> getAllVols();
    void setAllVols( QList<VolumeDefinition> vols );
    void createDefaultVol();

    QList<PhysicalBoundaryCondition> getAllPhysicalBoundaryConditions();
    void setAllPhysicalBoundaryConditions (QList<PhysicalBoundaryCondition> physical_boundary_conditions);
    void setAllPhysicalBoundaryConditions (QMap<QString, PhysicalBoundaryCondition> physical_boundary_conditions);
    bool physicalTypeDefined(QString name) { return m_PhysicalBoundaryConditionsMap.contains(name); };
    PhysicalBoundaryCondition getPhysicalBoundaryCondition(QString name) { return m_PhysicalBoundaryConditionsMap[name]; }

    static GuiMainWindow* pointer() { return THIS; }
    static void lock() { m_Mutex.lock(); }
    static void unlock() { m_Mutex.unlock(); }
    static bool tryLock() { return m_Mutex.tryLock(); }
    void getAllBoundaryCodes(QVector<int> &bcs);
    QSet<int> getAllBoundaryCodes();
    void getDisplayBoundaryCodes(QSet<int> &bcs);
    vtkPointPicker* getPointPicker() { return ( m_PointPicker );}
    vtkSphereSource* getPickSphere() { return ( m_PickSphere );}
    bool pickPoint( vtkIdType id_point );
    bool pickCell( vtkIdType id_cell );

    QString getFilename() { return( m_CurrentFilename ); }
    void setFilename(QString filename) { m_CurrentFilename = filename; }

    SurfaceProjection* getSurfProj(int bc, bool allow_null = false);
    void setSurfProj(SurfaceProjection *surf_proj, int bc) { m_SurfProj[bc] = surf_proj; }
    void setUniversalSurfProj(SurfaceProjection *surf_proj);
    bool checkSurfProj();

    void setSystemOutput();
    void setLogFileOutput();

    void resetProgress(QString info_text, int p_max);
    void setProgress(int p);

  public slots:

    void setUseVTKInteractor( int a_UseVTKInteractor );
    void setPickMode( bool a_UseVTKInteractor, bool a_CellPickerMode );

    void exit();                           ///< Exit the application
    void importSTL();                      ///< Import an STL file (ASCII or binary)
    void importGmsh1Ascii();               ///< Import a Gmsh grid from an ASCII file -- using version 1.0 of the Gmsh file format
    void exportGmsh1Ascii();               ///< Export a grid from to an ASCII Gmsh file -- using version 1.0 of the Gmsh file format
    void importGmsh2Ascii();               ///< Import a Gmsh grid from an ASCII file -- using version 2.0 of the Gmsh file format
    void exportGmsh2Ascii();               ///< Export a grid from to an ASCII Gmsh file -- using version 2.0 of the Gmsh file format
    void exportNeutral();                  ///< Export a grid to neutral format for NETGEN
    void updateActors( bool force = false ); ///< Update the VTK output
    void forceUpdateActors();              ///< Force an update of the VTK output
    void scaleToData();                    ///< Scale to data
    void zoomAll();                        ///< Move the camera in order to show everything on the screen
    void zoomOnPickedObject();
    void deselectAll();
    void printGrid() {cout << "PrintGrid() called!" << endl; cout_grid( cout, m_Grid, true, true, true, true );}
    void info();

    void undo();
    void redo();

    void resetOperationCounter();

    void open();                           ///< Open an existing case
    void open( QString file_name, bool update_current_filename = true ); ///< Open case file_name
    void save();                           ///< Save the current case
    void saveAs();                         ///< Save the current case -- using a different file name
    QString saveAs( QString file_name, bool update_current_filename = true );   ///< Save the current case as file_name. Returns name under which file was saved (with missing .egc extension for example).

    int  quickSave();                      ///< Save the current grid as a_filename_a_operation
    void quickLoad( int a_operation );     ///< Load a_filename_a_operation
    void updateStatusBar();                ///< Update the status bar
    void selectBoundaryCodes();            ///< Select the boundary codes to be displayed/hidden
    void updateBoundaryCodes( bool all_on ); ///< Update the boundary code book keeping (e.g. after reading a mesh).
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
    void exportAsciiPly();                 ///< Write surface elements to an ASCII PLY file.
    void exportBinaryPly();                ///< Write surface elements to a binary PLY file.
    void editBoundaryConditions();         ///< Edit boundary conditions (names and types)
    void configure();                      ///< Edit settings
    void about();                          ///< Display an about message
    void markOutputLine();                 ///< Mark the current position in the output window

    QString getXmlSection( QString name );                 ///< Get a section from the XML case description
    void    setXmlSection( QString name, QString contents ); ///< Set a section of the XML case description

    void viewRight();
    void viewLeft();
    void viewTop();
    void viewBottom();
    void viewFront();
    void viewBack();

    void appendOutput( QString txt ) { ui.textEditOutput->append( txt ); }
    void clearOutput() { ui.textEditOutput->clear(); }
    void updateOutput();
    void periodicUpdate();

    void storeSurfaceProjection(bool nosave = false);
    void resetSurfaceProjection();    

    // SLOTS for all standard operations should be defined below;
    // entries should look like this:
    // void callOperationName() { EG_STDSLOT(OperationName); };
    // The actual class name in this case, however, would be GuiOperationName.
    //
    // the following line can be used as a template:
    // void call() { EG_STDSLOT(); };
    // IMPORTANT: Using EG_STDSLOT sets lock_gui to true, while EG_STDINTERSLOT does not (default is lock_gui = false)
    // This is important to determine whether an operation should try to lock the main mutex or not.
    // If lock_gui is true, the operation will try to lock the main mutex. If it fails (mutex locked by other operation), the operation is stopped.
    // SUMMARY:
    // EG_STDSLOT = background operation (There can not be more than one background operation!)
    // EG_STDINTERSLOT = foreground operation
    // Note: In practice, EG_STDINTERSLOT locks everything, while EG_STDSLOT prevents other operations, but doesn't lock the text output or prevent minimizing the window.

    void callCreateSurfaceMesh() { EG_STDINTERSLOT( GuiCreateSurfaceMesh ); }
    void callCreateBoundaryLayer() { EG_STDSLOT( GuiCreateBoundaryLayer ); }
    void callDivideBoundaryLayer() { EG_STDSLOT( GuiDivideBoundaryLayer ); }
    void callDeleteVolumeGrid() { EG_STDSLOT( DeleteVolumeGrid ); }
    void callDeleteTetras() { EG_STDSLOT( DeleteTetras ); }
    void callCreateVolumeMesh() { EG_STDSLOT( GuiCreateVolumeMesh ); }
    void callSmoothVolumeGrid() { EG_STDSLOT( SmoothVolumeGrid ); }
    void callSetBoundaryCode()  { EG_STDINTERSLOT( GuiSetBoundaryCode ); }
    void callDeleteBadAspectTris() { EG_STDINTERSLOT( GuiDeleteBadAspectTris ); }
    void callDeletePickedCell() { EG_STDSLOT( DeletePickedCell ); }
    void callMergeNodes();
    void callInsertNewCell();
    void callDeletePickedPoint();
    void callBoxSelect() { EG_STDINTERSLOT( BoxSelect ); }
    void callCheckSurfaceIntegrity() { EG_STDINTERSLOT( CheckSurfaceIntegrity ); }
    void callPick_cell_point() { EG_STDINTERSLOT( GuiPick ); }
    void callTransform() { EG_STDINTERSLOT( GuiTransform ); }
    void callUpdateSurfProj() { EG_STDINTERSLOT( UpdateSurfProj ); }
    void callImportOpenFoamCase() { EG_STDREADERSLOT(FoamReader); }
    void callMergeVolumes() { EG_STDSLOT(GuiMergeVolumes); }
    void callMirrorMesh() { EG_STDSLOT(GuiMirrorMesh); }
    void callOrthogonalityOptimiser() { EG_STDSLOT(OrthogonalityOptimiser); }
    void callCreateHexCore() { EG_STDSLOT( GuiCreateHexCore ); }

    void callFixSTL();

    void callFoamWriter()                 { EG_STDINTERSLOT( FoamWriter ); }
    void callSimpleFoamWriter()           { EG_STDINTERSLOT( SimpleFoamWriter ); }
    void callFoamCaseWriter()             { EG_STDINTERSLOT( OpenFOAMcase ); }
    void callCgnsWriter()                 { EG_STDINTERSLOT( CgnsWriter ); }
    void callVtkReader()                  { EG_STDREADERSLOT( VtkReader ); }
    void callBlenderReader()              { EG_STDREADERSLOT( BlenderReader ); }
    void callBlenderWriter()              { EG_STDREADERSLOT( BlenderWriter ); }
    void callPolyDataReader()             { EG_STDREADERSLOT( PolyDataReader ); }
    void callReducedPolyDataReader()      { EG_STDREADERSLOT( ReducedPolyDataReader ); }
    void callSeligAirfoilReader()         { EG_STDREADERSLOT( SeligAirfoilReader ); }
    void callBrlcadReader()               { EG_STDREADERSLOT( BrlcadReader ); }
    void callExportSu2()                  { EG_STDREADERSLOT( Su2Writer ); }
    void callExportDolfyn()               { EG_STDREADERSLOT( DolfynWriter ); }

    void callSurfaceMesher()              { EG_STDSLOT(GuiSurfaceMesher); }
    void callReduceSurfaceTriangulation() { EG_STDSLOT(ReduceSurfaceTriangulation); }
    void callEliminateSmallBranches()     { EG_STDSLOT(EliminateSmallBranches); }
    void callSmoothAndSwapSurface()       { EG_STDSLOT(SmoothAndSwapSurface); }
    void callSharpenEdges()               { EG_STDSLOT(SharpenEdges); }
    void callCheckForOverlap()            { EG_STDSLOT(CheckForOverlap); }

    void callFixCADGeometry()             { EG_STDSLOT(FixCadGeometry); }

};

#endif
