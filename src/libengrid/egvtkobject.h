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
#ifndef egvtkobject_H
#define egvtkobject_H

class EgVtkObject;

class BezierTriangle;

#include "engrid.h"
#include "utilities.h"
#include "boundarycondition.h"

#include <vtkUnstructuredGrid.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkLongArray.h>
#include <vtkDoubleArray.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkCellLocator.h>
#include <vtkIdList.h>

#include <QSettings>
#include <QSet>
#include <QVector>

#define EG_SIMPLE_VERTEX         0
#define EG_FEATURE_EDGE_VERTEX   1
#define EG_BOUNDARY_EDGE_VERTEX  2
#define EG_FEATURE_CORNER_VERTEX 3
#define EG_FIXED_VERTEX          4

class EgVtkObject
{

public: // data-types

  typedef const QVector<vtkIdType>&     l2g_t;
  typedef const QVector<int>&           g2l_t;
  typedef const QVector<QVector<int> >& l2l_t;
  
private: // methods
  
  void addToC2C
    (
      vtkIdType               id_cell,
      QVector<int>           &_cells,
      QVector<QVector<int> > &c2c,
      int                     j,
      vtkIdList              *nds,
      vtkIdList              *cls,
      vtkUnstructuredGrid    *grid
    );
  
  void addToN2N
    (
      QVector<QSet<int> > &n2n,
      int                  n1,
      int                  n2
    );
  
  void createNodeField(vtkUnstructuredGrid *grid, QString field_name, QString type_name, int Nnodes, bool overwrite = false);
  void createCellField(vtkUnstructuredGrid *grid, QString field_name, QString type_name, int Ncells, bool overwrite = false);

protected: // attributes
  
  QSet<int> m_BoundaryCodes;
  static int DebugLevel;
  
protected: // methods
  
  /**
   * if key=value pair not found in settings file, write it + read key value from settings file and assign it to variable
   * Version for int variables
   */
  int getSet(QString group, QString key, int value, int& variable);
  
  /**
   * if key=value pair not found in settings file, write it + read key value from settings file and assign it to variable
   * Version for double variables
   */
  double getSet(QString group, QString key, double value, double& variable);
  
  /**
   * if key=value pair not found in settings file, write it + read key value from settings file and assign it to variable
   * Version for bool variables
   */
  bool getSet(QString group, QString key, bool value, bool& variable);
  
  /**
   * if key=value pair not found in settings file, write it + read key value from settings file and assign it to variable
   * Version for string variables
   */
  QString getSet(QString group, QString key, QString value, QString& variable);

  /**
   * if key=value pair not found in settings file, write it + read key value from settings file and assign it to variable
   * Version for string variables.
   */
  QString getSet(QString group, QString key, QString value, QString& variable, int type);
    
  /**
   * Update the cell index array.
   */
  void UpdateCellIndex(vtkUnstructuredGrid *grid);
  
  /**
   * Update the point index array.
   */
  void UpdateNodeIndex(vtkUnstructuredGrid *grid);
  
  /**
   * Compute normal vectors on nodes and cells of a subset of a grid.
   * The paramters nodes and cells must be consistent; this means the nodes
   * represent exactly (not more, not less) the nodes forming the cells.
   * @param cell_normals On return, this will contain the cell normals (same order as cells)
   * @param node_normals On return, this will contain the cell normals (same order as cells)
   * @param cells        The cells to compute the normals of
   * @param nodes        The nodes to compute the normals of
   * @param grid         The grid to operate on
   */
  void computeNormals(QVector<vec3_t> &cell_normals, QVector<vec3_t> &node_normals, QVector<vtkIdType> &cells, QVector<vtkIdType>  &nodes, vtkUnstructuredGrid *grid);
      
  /**
   * Create a mapping from global node indices to the indeces of a subset of nodes.
   * @param nodes  The subset of nodes.
   * @param _nodes On return, this will contain the mapping.
   * @param grid   The grid to operate on.
   */
  void createNodeMapping(QVector<vtkIdType> &nodes, QVector<int> &_nodes, vtkUnstructuredGrid *grid);
  
  /**
   * Create a mapping from global cell indices to the indices of a subset of cells.
   * @param cells  The subset of cells.
   * @param _cells On return, this will contain the mapping.
   * @param grid   The grid to operate on.
   */
  void createCellMapping(QVector<vtkIdType> &cells, QVector<int> &_cells, vtkUnstructuredGrid *grid);
  
  /**
   * Create a node to boundary condition ("cell_code") mapping.
   * Only non-zero boundary conditions will be considered.
   * @param bcs  On return, this will hold the codes of all boundary elements that are
   *             attached to a node.
   * @param grid The grid to operate on.
   */
  void createNodeToBcMapping(QVector<QSet<int> > &bcs, vtkUnstructuredGrid *grid);
  
  /**
   * Create a node to cell structure for a given set of cells and nodes.
   * This creates a vector of sets which might have performance issues.
   * @param cells  the subset of cells
   * @param nodes  the subset of nodes
   * @param _nodes the reverse mapping for the nodes
   * @param n2c    On return, this will hold the node to cell structure
   * @param grid   The grid to operate on
   */
  void createNodeToCell(QVector<vtkIdType> &cells, QVector<vtkIdType> &nodes, QVector<int> &_nodes, QVector<QSet<int> >  &n2c, vtkUnstructuredGrid *grid);
  
  /**
   * Create a node to cell structure for a given set of cells and nodes.
   * This creates a vector of vectors.
   * @param cells  the subset of cells
   * @param nodes  the subset of nodes
   * @param _nodes the reverse mapping for the nodes
   * @param n2c    On return, this will hold the node to cell structure
   * @param grid   The grid to operate on
   */
  void createNodeToCell(QVector<vtkIdType> &cells, QVector<vtkIdType> &nodes, QVector<int> &_nodes, QVector<QVector<int> > &n2c, vtkUnstructuredGrid *grid);

  /**
   * Create a node to node structure for a given set of cells and nodes.
   * This creates a vector of sets which might have performance issues.
   * @param cells  the subset of cells
   * @param nodes  the subset of nodes
   * @param _nodes the reverse mapping for the nodes
   * @param n2n    On return, this will hold the node to node structure
   * @param grid   The grid to operate on
   */
  void createNodeToNode(QVector<vtkIdType> &cells, QVector<vtkIdType> &nodes, QVector<int> &_nodes, QVector<QSet<int> > &n2n, vtkUnstructuredGrid *grid);
  
  /**
   * Create a node to node structure for a given set of cells and nodes.
   * This creates a vector of vectors.
   * @param cells  the subset of cells
   * @param nodes  the subset of nodes
   * @param _nodes the reverse mapping for the nodes
   * @param n2n    On return, this will hold the node to node structure
   * @param grid   The grid to operate on
   */
  void createNodeToNode(QVector<vtkIdType> &cells, QVector<vtkIdType> &nodes, QVector<int> &_nodes, QVector<QVector<int> > &n2n, vtkUnstructuredGrid *grid);

  /**
   * Extract the nodes which are part of a given set of cells.
   * @param cells the subset of cells
   * @param nodes On return, this will contain the nodes that correspond to the subset of cells
   * @param grid  The grid to operate on
   */
  template <class C>
  void getNodesFromCells(const C &cells, QVector<vtkIdType> &nodes, vtkUnstructuredGrid *grid);
  
  /**
   * Check if a cell is a volume cell.
   * @param cellId The id fof the cell in question
   * @param grid   The grid to operate on
   * @return true if the cell represents a volume and false if not
   */
  bool isVolume(vtkIdType id_cell, vtkUnstructuredGrid *grid);
  
  
  /**
   * Check if a cell is a surface cell.
   * @param cellId The id fof the cell in question
   * @param grid   The grid to operate on
   * @return true if the cell represents a surface and false if not
   */
  bool isSurface(vtkIdType id_cell, vtkUnstructuredGrid *grid);
  
  /**
   * Get all volume cells of a grid.
   * @param cells On return this will hold the Ids of all volume cells.
   * @param grid  The grid to operate on.
   */
  void getAllVolumeCells(QVector<vtkIdType> &cells, vtkUnstructuredGrid *grid);
  
  /**
   * Get all cells of a grid.
   * @param cells On return this will hold the Ids of all cells.
   * @param grid  The grid to operate on.
   */
  void getAllCells(QVector<vtkIdType> &cells, vtkUnstructuredGrid *grid);
  
  /**
   * Get all cells of a grid and a specific type.
   * @param type  The type of the cells (e.g. VTK_TETRA, VTK_TRIANGLE, etc.)
   * @param cells On return this will hold the Ids of all cells.
   * @param grid  The grid to operate on.
   */
  void getAllCellsOfType(vtkIdType type, QVector<vtkIdType> &cells, vtkUnstructuredGrid *grid);
  
  /**
   * Get all surface cells of a grid.
   * @param cells On return this will hold the Ids of all surface cells.
   * @param grid  The grid to operate on.
   */
  void getAllSurfaceCells(QVector<vtkIdType>  &cells, vtkUnstructuredGrid *grid);
  
  /**
   * Get all surface cells of a grid with a specific boundary condition.
   * @param bcs   The set of boundary conditions
   * @param cells On return this will hold the Ids of the surface cells.
   * @param grid  The grid to operate on.
   */
  void getSurfaceCells(QSet<int> &bcs, QVector<vtkIdType> &cells, vtkUnstructuredGrid *grid);
  
  /**
   * Create a cell neighbourship list for a subset grid. 
   * This has been implemented using VTK's vtkCellLinks structures.
   * @param cells the subset of cells
   * @param c2c   On return this will hold the neighbourship list
   * @param grid  The grid to operate on.
   */
  void createCellToCell(QVector<vtkIdType> &cells, QVector<QVector<int> > &c2c, vtkUnstructuredGrid *grid);
  
  /**
   * Insert a subset of a grid into a vtkPolyData structure.
   * This is can be used in order to make use of many readily available 
   * operations within VTK; one example is smoothing.
   * Cell index and node index arrays will be created and passed to the
   * poly data structure. Thus any purely geometric effect (no topology change)
   * can be directly reintroduced into the vtkUnstructuredGrid.
   * @param cells the subset of cells
   * @param pdata the vtkPolyData to add the nodes and cells to
   * @param grid  The grid to operate on.
   */
  void addToPolyData(QVector<vtkIdType>  &cells, vtkPolyData *pdata, vtkUnstructuredGrid *grid);
  
  /**
   * Copy the attributes from an input to an output cell.
   * @param old_grid the input grid
   * @param oldId    the existing input cell
   * @param new_grid the output grid
   * @param newId    the new output cell
   */
  void copyCellData(vtkUnstructuredGrid *old_grid, vtkIdType oldId, vtkUnstructuredGrid *new_grid, vtkIdType newId);
  
  /**
   * Copy the attributes from an input to an output node.
   * @param old_grid the input grid
   * @param oldId    the existing input node
   * @param new_grid the output grid
   * @param newId    the new output node
   */
  void copyNodeData(vtkUnstructuredGrid *old_grid, vtkIdType oldId, vtkUnstructuredGrid *new_grid, vtkIdType newId);
  
  /**
   * Create the basic fields on a given grid.
   * Care should be taken with the overwrite parameter; if it is set to <i>false</i>
   * and the field does not have the correct size it can lead to <i>segmentation faults</i>.
   * @param Ncells the number of output cells
   * @param Nnodes the number of output nodes
   * @param overwrite f set to true existing fields will be re-created
   */
  void createBasicFields(vtkUnstructuredGrid *grid, vtkIdType Ncells, vtkIdType Nnodes, bool overwrite = false);
  
  /**
   * Create the basic cell fields on a given grid.
   * Care should be taken with the overwrite parameter; if it is set to <i>false</i>
   * and the field does not have the correct size it can lead to <i>segmentation faults</i>.
   * @param Ncells the number of output cells
   * @param overwrite f set to true existing fields will be re-created
   */
  void createBasicCellFields(vtkUnstructuredGrid *grid, vtkIdType Ncells, bool overwrite = false);
  
  /**
   * Create the basic node fields on a given grid.
   * Care should be taken with the overwrite parameter; if it is set to <i>false</i>
   * and the field does not have the correct size it can lead to <i>segmentation faults</i>.
   * @param Nnodes the number of output nodes
   * @param overwrite f set to true existing fields will be re-created
   */
  void createBasicNodeFields(vtkUnstructuredGrid *grid, vtkIdType Nnodes, bool overwrite = false);
  
  /**
   * Allocate memory for a grid. This method will also create the basic
   * attribute fields (e.g. "cell_code").
   * @param grid   the grid for which to allocate memory
   * @param Ncells the number of output cells
   * @param Nnodes the number of output nodes
   * @param create_fields flag to determine if node and cell data shall be created
   */
  void allocateGrid(vtkUnstructuredGrid *grid, vtkIdType Ncells, vtkIdType Nnodes, bool create_fields = true);
  
  /**
   * Get the names of all node (point) attributes (fields) of a VTK grid
   * @param field_names On return this vector will contain the names of all fields
   * @param grid the grid
   */
  void getAllNodeDataNames(QVector<QString> &field_names, vtkUnstructuredGrid *grid);

  /**
   * Get the names of all cell attributes (fields) of a VTK grid
   * @param field_names On return this vector will contain the names of all fields
   * @param grid the grid
   */
  void getAllCellDataNames(QVector<QString> &field_names, vtkUnstructuredGrid *grid);

  /**
   * Compute the intersection of two Q containers.
   * This will return a set.
   * @param set1 the first container
   * @param set2 the second container
   * @param inters on return this will hold the intersection
   */
  template <class C1, class C2>
  void qcontIntersection(const C1& c1, const C2& c2, QSet<typename C1::value_type> &inters);

  /**
   * Compute the intersection of two Q containers.
   * This will return a vector.
   * @param set1 the first container
   * @param set2 the second container
   * @param inters on return this will hold the intersection
   */
  template <class C1, class C2>
  void qcontIntersection(const C1& c1, const C2& c2, QVector<typename C1::value_type> &inters);

  /**
   * Compute the centre of a cell
   * @param grid the grid to use
   * @param id_cell the cell of which to compute the centre
   * @return the centre of the cell
   */
  vec3_t cellCentre(vtkUnstructuredGrid *grid, vtkIdType id_cell);

  /**
   * Compute the angle between two faces.
   * @param grid the grid to use
   * @param id_face1 index of the first face
   * @param id_face2 index of the second face
   * @return the angle between the faces (M_PI for a flat surface)
   */
  double faceAngle(vtkUnstructuredGrid* grid, vtkIdType id_face1, vtkIdType id_face2);
  
  /**
   * Get the cells of a grid that are not part of a given set of cells.
   * @param grid the grid to use
   * @param cells the given set of cells
   * @param rest_cells on return this will hold the rest of all cells of the grid (not part of cells)
   */
  void getRestCells(vtkUnstructuredGrid *grid, const QVector<vtkIdType> &cells, QVector<vtkIdType> &rest_cells);
  
  /**
   * Find the corresponding volume cell of a surface cell
   * @param grid the grid to use
   * @param id_surf the id of the surface cell
   * @param n2n the node to cell structure for this grid
   * @return the id of the corresponding volume cell (or -1 if not found)
   */
  vtkIdType findVolumeCell(vtkUnstructuredGrid *grid, vtkIdType id_surf, g2l_t _nodes, l2g_t cells, g2l_t _cells, l2l_t n2c);

  /**
   * @brief Copy a cell from one grid to another.
   * @param src the source grid
   * @param id_cell the cell index within the source grid
   * @param dst the destination grid
   * @return the index of the cell in the destination grid
   */
  vtkIdType copyCell(vtkUnstructuredGrid* src, vtkIdType id_cell, vtkUnstructuredGrid* dst, vtkIdType offset = 0);

  /**
   * @brief Copy a cell from one grid to another and translate node indices.
   * @param src the source grid
   * @param id_cell the cell index within the source grid
   * @param dst the destination grid
   * @param src2dst a generic container containing the mapping from src to destination index
   * @return the index of the cell in the destination grid
   */
  template <typename C>
  vtkIdType copyCell(vtkUnstructuredGrid* src, vtkIdType id_cell, vtkUnstructuredGrid* dst, const C& src2dst);

  /**
   * Copy "src" grid to "dst" grid. Allocate "dst" so that it fits the data of "src".
   * @param src a pointer to the source grid
   * @param dst a pointer to the destination grid
   */
  void makeCopy(vtkUnstructuredGrid *src, vtkUnstructuredGrid *dst, bool copy_data = true);

  /**
   * Copy a part of "src" grid to "dst" grid. Allocate "dst" so that it fits the data to be copied.
   * @param src a pointer to the source grid
   * @param dst a pointer to the destination grid
   * @param cells a container with the cells to be copied
   */
  template <class C>
  void makeCopy(vtkUnstructuredGrid *src, vtkUnstructuredGrid *dst, const C &cells);

  /**
   * Copy "src" grid to "dst" grid. DO NOT allocate "dst" so that it fits the data of "src".
   * Allocation is left for the user to do.
   * @param src a pointer to the source grid
   * @param dst a pointer to the destination grid
   */
  void makeCopyNoAlloc(vtkUnstructuredGrid *src, vtkUnstructuredGrid *dst);

  /**
   * Change the orientation of a face.
   * @param grid the grid to use
   * @param id_face the id of the face to change
   */
  void reorientateFace(vtkUnstructuredGrid *grid, vtkIdType id_face);

  /**
   * Reset face orientation to original orientation.
   * @param grid the grid with the faces
   */
  void resetOrientation(vtkUnstructuredGrid *grid);
  
  void createIndices(vtkUnstructuredGrid *grid);
  
  /**
   * Get the boundary condition of a boundary code.
   * @param bc the boundary code
   * @return the boundary condition
   */
  BoundaryCondition getBC(int bc);
  
  /**
   * Save the subgrid defined by cls from grid.
   * @param grid The source grid
   * @param cls The cells to extract
   * @param file_name The file to save to
   */
  template <class C>
  void writeCells(vtkUnstructuredGrid *grid, const C &cls, QString file_name);
  
  /**
   * Get the SubGrid defined by cls from grid. The function takes care of allocation for SubGrid.
   * @param grid The source grid
   * @param cls The cells to extract
   * @param SubGrid The SubGrid to create 
   */
  template <class C>
  void getSubGrid(vtkUnstructuredGrid *grid, const C &cls, vtkUnstructuredGrid *SubGrid);

  
  /**
   * Save grid to file filename.
   * @param grid The source grid
   * @param filename Name of the file to save to
   */
  void writeGrid(vtkUnstructuredGrid *grid, QString filename);
  
  /**
   * Get a file name without extension.
   * @param file_name the full name (with extension)
   * @return the name without the extension
   */
  QString stripFromExtension(QString file_name);

  /**
   * Get the extension of a file name
   * @param file_name the full name (with extension)
   * @return the extension
   */
  QString getExtension(QString file_name);

  /**
   * Get a face of a cell
   * @param grid the unstructured grid
   * @param id_cell the index of the cell
   * @param i_face the index of the face within the cell
   * @param ids on return this will contain the nodes of the face
   */
  void getFaceOfCell(vtkUnstructuredGrid *grid, vtkIdType id_cell, int i_face, QVector<vtkIdType> &ids);

  /**
   * Get the normal of a face of a volume cell.
   * @param grid the unstructured grid
   * @param id_cell the index of the cell
   * @param i_face the index of the face within the cell
   * @return the normal vector (absolute value corresponds to the area)
   */
  vec3_t getNormalOfCell(vtkUnstructuredGrid *grid, vtkIdType id_cell, int i_face);

  /**
   * Get the centre of a face of a volume cell.
   * @param grid the unstructured grid
   * @param id_cell the index of the cell
   * @param i_face the index of the face within the cell
   * @return the centre
   */
  vec3_t getCentreOfCellFace(vtkUnstructuredGrid *grid, vtkIdType id_cell, int i_face);

  /**
   * Get an edge of a face/cell
   * @param grid the unstructured grid
   * @param id_cell the index of the cell
   * @param i_edge the index of the edge within the cell
   * @param ids on return this will contain the nodes of the edge
   */
  void getEdgeOfCell(vtkUnstructuredGrid *grid, vtkIdType id_cell, int i_edge, QVector<vtkIdType> &ids);

  /**
   * Get all boundary codes fo a grid.
   * @param grid the grid to extract the boundaru codes from
   * @return a set with all boundary codes
   */
  QSet<int> getAllBoundaryCodes(vtkUnstructuredGrid *grid);

  /**
   * Check if a cell (face) contains a given node
   * @param grid the unstructured grid
   * @param id_cell the id of the cell to investigate
   * @param id_node the id of the required node
   * @return true if id_cell contains id_node
   */
  bool cellContainsNode(vtkUnstructuredGrid *grid, vtkIdType id_cell, vtkIdType id_node);

  /**
   * Get all node indices which are shared by two cells.
   * These cells can be surface or volume cells; also a combination
   * of a volume and a surface cell is possible.
   * @param grid the grid to use
   * @param id_cell1 index of the first cell
   * @param id_cell2 index of the second cell
   * @param cont a generic Qt container which will hold the shared node indices on return
   */
  template <typename C>
  void sharedNodesOfCells(vtkUnstructuredGrid* grid, vtkIdType id_cell1, vtkIdType id_cell2, C& cont);

  /**
   * @brief Get all nodes of a cell from a vtkUnstructuredGrid
   * This methods collects all nodes of a cell in a generic Qt container.
   * It can be used uniformly for VTK_POLYHEDRON cells and standard cells.
   * @param grid the grid to use
   * @param id_cell the cell index
   * @param cont a generic Qt container which will hold node indices on return
   */
  template <typename C>
  void getPointsOfCell(vtkUnstructuredGrid* grid, vtkIdType id_cell, C& cont);

  template <class C> void createPolyData(const C &x, vtkPolyData *poly_data, bool closed_loop = false);
  void createPolyDataC2C(vtkPolyData *poly_data, QVector<QVector<vtkIdType> > &c2c);
  void createPolyDataN2C(vtkPolyData *poly_data, QVector<QSet<vtkIdType> > &n2c);
  void createPolyDataN2N(vtkPolyData *poly_data, QVector<QSet<vtkIdType> > &n2n);
  template <class C> double convexRatio(const C &x, vec3_t n_plane, bool closed_loop = false);

  template <class C> void invertQContainer(C &cont);


public: // methods
  
  EgVtkObject() { DebugLevel = 0; }

  void setBoundaryCodes(const QSet<int> &bcs);
  QSet<int> getBoundaryCodes();
  void setDebugLevel(int a_DebugLevel) { DebugLevel = a_DebugLevel; }
  
  bool saveGrid( vtkUnstructuredGrid* a_grid, QString file_name );

  vtkIdType addGrid(vtkUnstructuredGrid *main_grid, vtkUnstructuredGrid *grid_to_add, vtkIdType offset);

  template <typename C>
  void vtkIdList2QContainer(vtkIdList *id_list, C &cont);

  void checkGridConsitency(vtkUnstructuredGrid *grid);

  template <typename C1, typename C2>
  void triangulatePolygon(vtkUnstructuredGrid *grid, const C1 &polygon, C2 &triangles);
  
private:

  void addVtkTypeInfo(vtkUnstructuredGrid* a_grid); ///< Add VTK type information to the grid (useful for visualisation with ParaView).

};

//End of class EgVtkObject

template <class C1, class C2>
void EgVtkObject::qcontIntersection(const C1& c1, const C2& c2, QSet<typename C1::value_type> &inters)
{
  inters.clear();
  foreach (typename C1::value_type t1, c1) {
    foreach (typename C2::value_type t2, c2) {
      if (t1 == t2) {
        inters.insert(t1);
      }
    }
  }
}

template <class C1, class C2>
void EgVtkObject::qcontIntersection(const C1& c1, const C2& c2, QVector<typename C1::value_type> &inters)
{
  QSet<typename C1::value_type> inters_set;
  qcontIntersection(c1, c2, inters_set);
  inters.resize(inters_set.size());
  qCopy(inters_set.begin(), inters_set.end(), inters.begin());
}

template <class C>
void EgVtkObject::getSubGrid(vtkUnstructuredGrid *grid, const C &cls, vtkUnstructuredGrid *sub_grid)
{
  createIndices(grid);
  QVector<vtkIdType> cells;
  QVector<vtkIdType> nodes;
  cells.resize(cls.size());
  qCopy(cls.begin(), cls.end(), cells.begin());
  getNodesFromCells(cells, nodes, grid);
  allocateGrid(sub_grid, cells.size(), nodes.size());
  vtkIdType id_new_node = 0;
  QVector<vtkIdType> old2new(grid->GetNumberOfPoints(), -1);
  foreach (vtkIdType id_node, nodes) {
    vec3_t x;
    grid->GetPoint(id_node, x.data());
    sub_grid->GetPoints()->SetPoint(id_new_node, x.data());
    old2new[id_node] = id_new_node;
    copyNodeData(grid, id_node, sub_grid, id_new_node);
    ++id_new_node;
  }
  foreach (vtkIdType id_cell, cells) {
    vtkIdType id_new_cell = copyCell(grid, id_cell, sub_grid, old2new);
    copyCellData(grid, id_cell, sub_grid, id_new_cell);
  }
}

template <class C>
void EgVtkObject::writeCells(vtkUnstructuredGrid *grid, const C &cls, QString file_name)
{
  qDebug()<<"Saving cells from grid as "<<file_name;
  
  EG_VTKSP(vtkUnstructuredGrid,SubGrid);
  getSubGrid(grid,cls,SubGrid);
  
  EG_VTKSP(vtkXMLUnstructuredGridWriter,vtu);
  vtu->SetFileName(qPrintable(file_name));
  vtu->SetDataModeToBinary();
  vtu->SetInputData(SubGrid);
  vtu->Write();
}

template <class C>
void EgVtkObject::getNodesFromCells(const C& cells, QVector<vtkIdType>  &nodes, vtkUnstructuredGrid *grid)
{
  QSet<vtkIdType> ex_nodes;
  vtkIdType id_cell;
  foreach(id_cell, cells) {
    QList<vtkIdType> pts;
    getPointsOfCell(grid, id_cell, pts);
    foreach (vtkIdType id_node, pts) {
      if (id_node >= grid->GetNumberOfPoints()) {
        EG_BUG;
      }
      ex_nodes.insert(id_node);
    }
  }
  nodes.resize(ex_nodes.size());
  {
    int j = 0;
    vtkIdType i;
    foreach(i,ex_nodes) {
      nodes[j] = i;
      ++j;
    }
  }
}

inline vtkIdType EgVtkObject::copyCell(vtkUnstructuredGrid *src, vtkIdType id_cell, vtkUnstructuredGrid *dst, vtkIdType offset)
{
  EG_VTKSP(vtkIdList, stream);
  vtkIdType type_cell = src->GetCellType(id_cell);
  vtkIdType id_new_cell = -1;
  if (type_cell == VTK_POLYHEDRON) {
    src->GetFaceStream(id_cell, stream);
    int id = 1;
    for (int i = 0; i < stream->GetId(0); ++i) {
      int num_pts = stream->GetId(id);
      ++id;
      for (int j = 0; j < num_pts; ++j) {
        stream->SetId(id, stream->GetId(id) + offset);
        ++id;
      }
    }
    id_new_cell = dst->InsertNextCell(type_cell, stream);
  } else {
    src->GetCellPoints(id_cell, stream);
    for (int i = 0; i < stream->GetNumberOfIds(); ++i) {
      stream->SetId(i, stream->GetId(i) + offset);
    }
    id_new_cell = dst->InsertNextCell(type_cell, stream);
  }
  return id_new_cell;
}

template <typename C>
vtkIdType EgVtkObject::copyCell(vtkUnstructuredGrid *src, vtkIdType id_cell, vtkUnstructuredGrid *dst, const C& src2dst)
{
  EG_VTKSP(vtkIdList, stream);
  vtkIdType type_cell = src->GetCellType(id_cell);
  vtkIdType id_new_cell = -1;
  if (type_cell == VTK_POLYHEDRON) {
    src->GetFaceStream(id_cell, stream);
    int id = 1;
    for (int i = 0; i < stream->GetId(0); ++i) {
      int num_pts = stream->GetId(id);
      ++id;
      for (int j = 0; j < num_pts; ++j) {
        stream->SetId(id, src2dst[stream->GetId(id)]);
        ++id;
      }
    }
    id_new_cell = dst->InsertNextCell(type_cell, stream);
  } else {
    src->GetCellPoints(id_cell, stream);
    for (int i = 0; i < stream->GetNumberOfIds(); ++i) {
      stream->SetId(i, src2dst[stream->GetId(i)]);
    }
    id_new_cell = dst->InsertNextCell(type_cell, stream);
  }
  return id_new_cell;
}

template <class C>
void EgVtkObject::makeCopy(vtkUnstructuredGrid *src, vtkUnstructuredGrid *dst, const C& cells)
{
  QVector<vtkIdType> nodes;
  getNodesFromCells(cells, nodes, src);
  allocateGrid(dst, cells.size(), nodes.size());
  vtkIdType id_new_node = 0;
  QVector<vtkIdType> old2new(src->GetNumberOfPoints(), -1);
  foreach (vtkIdType id_node, nodes) {
    vec3_t x;
    src->GetPoints()->GetPoint(id_node, x.data());
    dst->GetPoints()->SetPoint(id_new_node, x.data());
    copyNodeData(src, id_node, dst, id_new_node);
    old2new[id_node] = id_new_node;
    ++id_new_node;
  }
  foreach (vtkIdType id_cell, cells) {
    vtkIdType id_new_cell = copyCell(src, id_cell, dst, old2new);
    copyCellData(src, id_cell, dst, id_new_cell);
  }
}

template <class C>
void EgVtkObject::createPolyData(const C &x, vtkPolyData *poly_data, bool closed_loop)
{
  int N = x.size();
  if (closed_loop) {
    --N;
  }
  EG_VTKSP(vtkPoints, points);
  points->SetNumberOfPoints(N);
  QVector<vtkIdType> pts(N);
  for (int i = 0; i < N; ++i) {
    points->SetPoint(i, x[i][0], x[i][1], x[i][2]);
    pts[i] = i;
  }
  poly_data->Allocate(1);
  poly_data->SetPoints(points);
  poly_data->InsertNextCell(VTK_POLYGON, N, pts.data());
}

template <class C>
double EgVtkObject::convexRatio(const C &x, vec3_t n_plane, bool closed_loop)
{
  double L_max = -1e99;
  double L_min =  1e99;
  int N = x.size();
  if (closed_loop) {
    --N;
  }
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) {
      int p1 = j;
      int p2 = j+1;
      if (j == N - 1 && !closed_loop) {
        p2 = 0;
      }
      if (i != p1 && i != p2) {
        vec3_t n = n_plane.cross(x[p2] - x[p1]);
        n.normalise();
        double L = (x[i] - x[j])*n;
        L_max = max(L, L_max);
        L_min = min(L, L_min);
      }
    }
  }
  return L_min/L_max;
}

template <typename C>
void EgVtkObject::sharedNodesOfCells(vtkUnstructuredGrid* grid, vtkIdType id_cell1, vtkIdType id_cell2, C& cont)
{
  cont.clear();
  EG_GET_CELL(id_cell1, grid);
  for (int i = 0; i < num_pts; ++i) {
    if (cellContainsNode(grid, id_cell2, pts[i])) {
      cont << pts[i];
    }
  }
}

template <typename C>
void EgVtkObject::getPointsOfCell(vtkUnstructuredGrid* grid, vtkIdType id_cell, C& cont)
{
  cont.clear();
  EG_VTKSP(vtkIdList, stream);
  vtkIdType type_cell = grid->GetCellType(id_cell);
  if (type_cell == VTK_POLYHEDRON) {
    grid->GetFaceStream(id_cell, stream);
    vtkIdType id = 1;
    for (int i = 0; i < stream->GetId(0); ++i) {
      int num_pts = stream->GetId(id++);
      for (int j = 0; j < num_pts; ++j) {
        cont << stream->GetId(id);
        ++id;
      }
    }
  } else {
    grid->GetCellPoints(id_cell, stream);
    for (vtkIdType i = 0; i < stream->GetNumberOfIds(); ++i) {
      cont << stream->GetId(i);
    }
  }
}

template <typename C>
void EgVtkObject::vtkIdList2QContainer(vtkIdList *id_list, C &cont)
{
  cont.clear();
  for (int i = 0; i < id_list->GetNumberOfIds(); ++i) {
    cont << id_list->GetId(i);
  }
}

template <typename C1, typename C2>
void EgVtkObject::triangulatePolygon(vtkUnstructuredGrid *grid, const C1 &polygon, C2 &triangles)
{
  int N = polygon.size();
  if (N < 3) {
    EG_BUG;
  }
  QVector<vtkIdType> poly(N+4);
  for (int i = 2; i <= N+1; ++i) {
    poly[i] = polygon[i-2];
  }
  poly[0]   = poly[N];
  poly[1]   = poly[N+1];
  poly[N+2] = poly[2];
  poly[N+3] = poly[3];

  int    i_best     = 0;
  double best_angle = M_PI;

  for (int i = 2; i <= N+1; ++i) {
    double angle = 0;
    vec3_t a, b, c;
    grid->GetPoint(poly[i], a.data());
    for (int j = 0; j <= N+1; ++j) {
      if (j < i-1 || j > i) {
        grid->GetPoint(poly[j],   b.data());
        grid->GetPoint(poly[j+1], c.data());
        angle = max(angle, GeometryTools::angle(b-a, c-a));
        angle = max(angle, GeometryTools::angle(a-b, c-b));
        angle = max(angle, GeometryTools::angle(a-c, b-c));
      }
    }
    if (angle < best_angle) {
      i_best = i;
      best_angle = angle;
    }
  }
  triangles.resize(N-2);
  int i1 = 0;
  int i2 = 1;
  if (i1 == i_best || i1 == i_best - N) {
    ++i1;
    ++i2;
  }
  for (int i = 0; i < triangles.size(); ++i) {
    triangles[i].resize(3);
    if (i2 == i_best || i2 == i_best - N) {
      i1 += 2;
      i2 += 2;
    }
    triangles[i][0] = poly[i1];
    triangles[i][1] = poly[i2];
    triangles[i][2] = poly[i_best];
    ++i1;
    ++i2;
  }
}

template <typename C>
void EgVtkObject::invertQContainer(C &cont)
{
  QList<typename C::value_type> original;
  foreach (typename C::value_type item, cont) {
    original << item;
  }
  cont.clear();
  while (original.size() > 0) {
    cont << original.last();
    original.pop_back();
  }
}

#endif
