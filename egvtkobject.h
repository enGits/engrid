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
#ifndef egvtkobject_H
#define egvtkobject_H

class EgVtkObject;

#include "engrid.h"
#include "boundarycondition.h"

#include <vtkUnstructuredGrid.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkLongArray.h>
#include <vtkDoubleArray.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <QSettings>

class EgVtkObject
{
  
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
  
protected: // attributes
  
  QSet<int> boundary_codes;
      
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
  void computeNormals
    (
      QVector<vec3_t>     &cell_normals,
      QVector<vec3_t>     &node_normals,
      QVector<vtkIdType>  &cells,
      QVector<vtkIdType>  &nodes,
      vtkUnstructuredGrid *grid
    );
      
  /**
   * Create a mapping from global node indices to the indeces of a subset of nodes.
   * @param nodes  The subset of nodes.
   * @param _nodes On return, this will contain the mapping.
   * @param grid   The grid to operate on.
   */
  void createNodeMapping
    (
      QVector<vtkIdType>  &nodes,
      QVector<int>        &_nodes,
      vtkUnstructuredGrid *grid
    );
  
  /**
   * Create a mapping from global cell indices to the indeces of a subset of cells.
   * @param cells  The subset of cells.
   * @param _cells On return, this will contain the mapping.
   * @param grid   The grid to operate on.
   */
  void createCellMapping
    (
      QVector<vtkIdType>  &cells,
      QVector<int>        &_cells,
      vtkUnstructuredGrid *grid
    );
  
  /**
   * Create a node to boundary condition ("cell_code") mapping.
   * Only non-zero boundary conditions will be considered.
   * @param bcs  On return, this will hold the codes of all boundary elements that are
   *             attached to a node.
   * @param grid The grid to operate on.
   */
  void createNodeToBcMapping
    (
      QVector<QSet<int> > &bcs,
      vtkUnstructuredGrid *grid
    );
  
  /**
   * Create a node to cell structure for a given set of cells and nodes.
   * @param cells  the subset of cells
   * @param nodes  the subset of nodes
   * @param _nodes the reverse mapping for the nodes
   * @param n2c    On return, this will hold the node to cell structure
   * @param grid   The grid to operate on
   */
  void createNodeToCell
    (
      QVector<vtkIdType>  &cells,
      QVector<vtkIdType>  &nodes,
      QVector<int>        &_nodes,
      QVector<QSet<int> > &n2c,
      vtkUnstructuredGrid *grid
    );
  
  /**
   * Create a node to node structure for a given set of cells and nodes.
   * @param cells  the subset of cells
   * @param nodes  the subset of nodes
   * @param _nodes the reverse mapping for the nodes
   * @param n2c    On return, this will hold the node to node structure
   * @param grid   The grid to operate on
   */
  void createNodeToNode
    (
      QVector<vtkIdType>  &cells,
      QVector<vtkIdType>  &nodes,
      QVector<int>        &_nodes,
      QVector<QSet<int> > &n2n,
      vtkUnstructuredGrid *grid
    );
  
  /**
   * Extract the nodes which are part of a given set of cells.
   * @param cells the subset of cells
   * @param nodes On return, this will contain the nodes that correspond to the subset of cells
   * @param grid  The grid to operate on
   */
  void getNodesFromCells
    (
      QVector<vtkIdType>  &cells,
      QVector<vtkIdType>  &nodes,
      vtkUnstructuredGrid *grid
    );
  
  /**
   * Check if a cell is a volume cell.
   * @param cellId The id fof the cell in question
   * @param grid   The grid to operate on
   * @return true if the cell represents a volume and false if not
   */
  bool isVolume
    (
      vtkUnstructuredGrid *grid,
      vtkIdType            cellId
    );
  bool isVolume(vtkIdType id_cell, vtkUnstructuredGrid *grid) { return isVolume(grid, id_cell); };
  
  
  /**
   * Check if a cell is a surface cell.
   * @param cellId The id fof the cell in question
   * @param grid   The grid to operate on
   * @return true if the cell represents a surface and false if not
   */
  bool isSurface
    (
      vtkUnstructuredGrid *grid,
      vtkIdType            id_cell
    );
  bool isSurface(vtkIdType id_cell, vtkUnstructuredGrid *grid) { return isSurface(grid, id_cell); };
  
  /**
   * Get all volume cells of a grid.
   * @param cells On return this will hold the Ids of all volume cells.
   * @param grid  The grid to operate on.
   */
  void getAllVolumeCells
    (
      QVector<vtkIdType>  &cells,
      vtkUnstructuredGrid *grid
    );
  
  /**
   * Get all cells of a grid.
   * @param cells On return this will hold the Ids of all cells.
   * @param grid  The grid to operate on.
   */
  void getAllCells
    (
      QVector<vtkIdType>  &cells,
      vtkUnstructuredGrid *grid
    );
  
  /**
   * Get all cells of a grid and a specific type.
   * @param type  The type of the cells (e.g. VTK_TETRA, VTK_TRIANGLE, etc.)
   * @param cells On return this will hold the Ids of all cells.
   * @param grid  The grid to operate on.
   */
  void getAllCellsOfType
    (
      vtkIdType            type,
      QVector<vtkIdType>  &cells,
      vtkUnstructuredGrid *grid
    );
  
  /**
   * Get all surface cells of a grid.
   * @param cells On return this will hold the Ids of all surface cells.
   * @param grid  The grid to operate on.
   */
  void getAllSurfaceCells
    (
      QVector<vtkIdType>  &cells,
      vtkUnstructuredGrid *grid
    );
  
  /**
   * Get all surface cells of a grid with a specific boundary condition.
   * @param bcs   The set of boundary conditions
   * @param cells On return this will hold the Ids of the surface cells.
   * @param grid  The grid to operate on.
   */
  void getSurfaceCells
    (
      QSet<int>           &bcs,
      QVector<vtkIdType>  &cells,
      vtkUnstructuredGrid *grid
    );
  
  /**
   * Get all surface nodes of a grid with a specific boundary condition.
   * @param bcs   The set of boundary conditions
   * @param nodes On return this will hold the Ids of the surface nodes.
   * @param grid  The grid to operate on.
   */
  void getSurfaceNodes
    (
      QSet<int>           &bcs,
      QSet <vtkIdType> &SelectedNodes,
      vtkUnstructuredGrid *grid
    );
    
  /**
   * Create a cell neighbourship list for a subset grid. 
   * This has been implemented using VTK's vtkCellLinks structures.
   * @param cells the subset of cells
   * @param c2c   On return this will hold the neighbourship list
   * @param grid  The grid to operate on.
   */
  void createCellToCell
    (
      QVector<vtkIdType>     &cells,
      QVector<QVector<int> > &c2c,
      vtkUnstructuredGrid    *grid
    );
  
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
  void addToPolyData
    (
      QVector<vtkIdType>  &cells,
      vtkPolyData         *pdata,
      vtkUnstructuredGrid *grid
    );
  
  /**
   * Copy the attributes from an input to an output cell.
   * @param old_grid the input grid
   * @param oldId    the existing input cell
   * @param new_grid the output grid
   * @param newId    the new output cell
   */
  void copyCellData
    (
      vtkUnstructuredGrid *old_grid, 
      vtkIdType            oldId, 
      vtkUnstructuredGrid *new_grid,
      vtkIdType            newId
    );
  
  /**
   * Copy the attributes from an input to an output node.
   * @param old_grid the input grid
   * @param oldId    the existing input node
   * @param new_grid the output grid
   * @param newId    the new output node
   */
  void copyNodeData
    (
      vtkUnstructuredGrid *old_grid, 
      vtkIdType            oldId, 
      vtkUnstructuredGrid *new_grid,
      vtkIdType            newId
    );
  
  /**
   * Create the basic fields on a given grid.
   * Care should be taken with the overwrite parameter; if it is set to <i>false</i>
   * and the field does not have the correct size it can lead to <i>segmentation faults</i>.
   * @param Ncells the number of output cells
   * @param Nnodes the number of output nodes
   * @param overwrite f set to true existing fields will be re-created
   */
  void createBasicFields
    (
      vtkUnstructuredGrid *grid,
      vtkIdType            Ncells, 
      vtkIdType            Nnodes,
      bool                 overwrite = true
    );
  
  /**
   * Create the basic cell fields on a given grid.
   * Care should be taken with the overwrite parameter; if it is set to <i>false</i>
   * and the field does not have the correct size it can lead to <i>segmentation faults</i>.
   * @param Ncells the number of output cells
   * @param overwrite f set to true existing fields will be re-created
   */
  void createBasicCellFields
    (
      vtkUnstructuredGrid *grid,
      vtkIdType            Ncells, 
      bool                 overwrite = true
    );
  
  /**
   * Create the basic node fields on a given grid.
   * Care should be taken with the overwrite parameter; if it is set to <i>false</i>
   * and the field does not have the correct size it can lead to <i>segmentation faults</i>.
   * @param Nnodes the number of output nodes
   * @param overwrite f set to true existing fields will be re-created
   */
  void createBasicNodeFields
    (
      vtkUnstructuredGrid *grid,
      vtkIdType            Nnodes,
      bool                 overwrite = true
    );
  
  /**
   * Allocate memory for a grid. This method will also create the basic
   * attribute fields (e.g. "cell_code").
   * @param grid   the grid for which to allocate memory
   * @param Ncells the number of output cells
   * @param Nnodes the number of output nodes
   */
  void allocateGrid
    (
      vtkUnstructuredGrid *grid,
      vtkIdType            Ncells,
      vtkIdType            Nnodes
    );
  
  /**
   * Compute the intersection of two QSets.
   * @param set1 the first set
   * @param set2 the second set
   * @param inters On return this will hold the intersection
   */
  template <class T>
  void setIntersection(const QSet<T> &set1,
                       const QSet<T> &set2,
                       QSet<T>       &inters);
    
  /**
   * Compute the intersection of two QVectors.
   * @param set1 the first vector
   * @param set2 the second vector
   * @param inters On return this will hold the intersection
   */
  template <class T>
  void vectorIntersection(const QVector<T> &set1,
                          const QVector<T> &set2,
                          QSet<T>          &inters);
  
  /**
   * Compute the centre of a cell
   * @param grid the grid to use
   * @param id_cell the cell of which to compute the centre
   * @return the centre of the cell
   */
  vec3_t cellCentre(vtkUnstructuredGrid *grid, vtkIdType id_cell);
  
  /**
   * Get the cells of a grid that are not part of a given set of cells.
   * @param grid the grid to use
   * @param cells the given set of cells
   * @param rest_cells on return this will hold the rest of all cells of the grid (not part of cells)
   */
  void getRestCells(vtkUnstructuredGrid      *grid, 
                    const QVector<vtkIdType> &cells,
                    QVector<vtkIdType>       &rest_cells);
  
  /**
   * Find the corresponding volume cell of a surface cell
   * @param grid the grid to use
   * @param id_surf the id of the surface cell
   * @param n2n the node to cell structure for this grid
   * @return the id of the corresponding volume cell (or -1 if not found)
   */
  int findVolumeCell
    (
      vtkUnstructuredGrid     *grid,
      vtkIdType                id_surf,
      const QVector<int>      _nodes,      
      const QVector<vtkIdType> cells,      
      const QVector<int>      _cells,      
      QVector<QSet<int> >     &n2c
    );

  /**
   * Copy "src" grid to "dst" grid. Allocate "dst" so that it fits the data of "src".
   */
  void makeCopy(vtkUnstructuredGrid *src, vtkUnstructuredGrid *dst);
  /**
   * Copy "src" grid to "dst" grid. DO NOT allocate "dst" so that it fits the data of "src".
   * Allocation is left for the user to do.
   */
  void makeCopyNoAlloc(vtkUnstructuredGrid *src, vtkUnstructuredGrid *dst);
  
  void createIndices(vtkUnstructuredGrid *grid);
  
  /**
   * Get the boundary condition of a boundary code.
   * @param bc the boundary code
   * @return the boundary condition
   */
  BoundaryCondition getBC(int bc);
  
  template <class T>
  void writeCells(vtkUnstructuredGrid *grid, const T &cls, QString file_name);
  
public: // methods
  
  void setBoundaryCodes(const QSet<int> &bcs);
  
};


template <class T>
void EgVtkObject::setIntersection(const QSet<T> &set1,
                                  const QSet<T> &set2,
                                  QSet<T>       &inters)
{
  inters.clear();
  foreach (T t1, set1) {
    if (set2.contains(t1)) {
      inters.insert(t1);
    };
  };
};

template <class T>
void EgVtkObject::vectorIntersection(const QVector<T> &set1,
                                     const QVector<T> &set2,
                                     QSet<T>          &inters)
{
  inters.clear();
  foreach (T t1, set1) {
    if (set2.has(t1)) {
      inters.insert(t1);
    };
  };
};

template <class T>
void EgVtkObject::writeCells(vtkUnstructuredGrid *grid, const T &cls, QString file_name)
{
  createIndices(grid);
  EG_VTKSP(vtkUnstructuredGrid,tmp_grid);
  QVector<vtkIdType> cells;
  QVector<vtkIdType> nodes;
  cells.resize(cls.size());
  qCopy(cls.begin(), cls.end(), cells.begin());
  getNodesFromCells(cells, nodes, grid);
  allocateGrid(tmp_grid, cells.size(), nodes.size());
  vtkIdType id_new_node = 0;
  QVector<vtkIdType> old2new(grid->GetNumberOfPoints(), -1);
  foreach (vtkIdType id_node, nodes) {
    vec3_t x;
    grid->GetPoint(id_node, x.data());
    tmp_grid->GetPoints()->SetPoint(id_new_node, x.data());
    old2new[id_node] = id_new_node;
    copyNodeData(grid, id_node, tmp_grid, id_new_node);
    ++id_new_node;
  };
  foreach (vtkIdType id_cell, cells) {
    vtkIdType N_pts, *pts;
    vtkIdType type_cell = grid->GetCellType(id_cell);
    grid->GetCellPoints(id_cell, N_pts, pts);
    QVector<vtkIdType> new_pts(N_pts);
    for (int i_pts = 0; i_pts < N_pts; ++i_pts) {
      new_pts[i_pts] = old2new[pts[i_pts]];
    };
    vtkIdType id_new_cell = tmp_grid->InsertNextCell(type_cell, N_pts, new_pts.data());
    copyCellData(grid, id_cell, tmp_grid, id_new_cell);
  };
  EG_VTKSP(vtkXMLUnstructuredGridWriter,vtu);
  vtu->SetFileName(file_name.toAscii().data());
  vtu->SetDataModeToBinary();
  vtu->SetInput(tmp_grid);
  vtu->Write();
};

  /**
   * Utility function that allows printing selected data from an vtkUnstructuredGrid to any ostream (includes ofstream objects)
   * @param stream ostream object to print to
   * @param grid vtkUnstructuredGrid you want to print data from
   * @param npoints print number of points in the grid
   * @param ncells print number of cells in the grid
   * @param points print points in the grid
   * @param cells print cells in the grid
   */
int cout_grid(ostream &stream, vtkUnstructuredGrid *grid, bool npoints=true, bool ncells=true, bool points=false, bool cells=false);

///////////////////////////////////////////
int addPoint(vtkUnstructuredGrid* a_grid,vtkIdType index,vec3_t x);
int addCell(vtkUnstructuredGrid* a_grid, vtkIdType A, vtkIdType B, vtkIdType C, int bc);
int getShortestSide(vtkIdType a_id_cell,vtkUnstructuredGrid* a_grid);
int getLongestSide(vtkIdType a_id_cell,vtkUnstructuredGrid* a_grid);
QSet <int> complementary_bcs(QSet <int> &bcs, vtkUnstructuredGrid *a_grid, QVector <vtkIdType> &a_cells);
QString cell2str(vtkIdType id_cell,vtkUnstructuredGrid* grid);
Qt::CheckState int2CheckState(int a);
int CheckState2int(Qt::CheckState a);

///////////////////////////////////////////
template <class T>
ostream &operator<<(ostream &out, QVector<T> & vector)
{
  int N=vector.size();
  out<<"[";
  for (int i = 0; i < N; ++i) {
    out<<vector.at(i);
    if(i!=N-1) out<<",";
  }
  out<<"]";
  return(out);
}

template <class T>
ostream &operator<<(ostream &out, QSet<T> const & set )
{
  out << "[ ";
  foreach (T value, set) out << value << " ";
  out << "]";
  return(out);
}

template <class T>
ostream &operator<<(ostream &out, QVector<QSet<T> > & vector)
{
  int N=vector.size();
  out<<"[";
  for (int i = 0; i < N; ++i) {
    QSet<T> set=vector.at(i);
    out<<set;
    if(i!=N-1) out<<",";
  }
  out<<"]";
  return(out);
}

template <class T>
ostream &operator<<(ostream &out, QVector<QVector<T> > & vector)
{
  int N=vector.size();
  out<<"[";
  for (int i = 0; i < N; ++i) {
    QVector<T> subvector=vector.at(i);
    out<<subvector;
    if(i!=N-1) out<<",";
  }
  out<<"]";
  return(out);
}
///////////////////////////////////////////
// ///////////////////////////////////////////
// /* Here is how we we get QTextStreams that look like iostreams */
// QTextStream Qcin;
// QTextStream Qcout;
// QTextStream Qcerr;
// ///////////////////////////////////////////

pair<vtkIdType,vtkIdType> OrderedPair(vtkIdType a, vtkIdType b);

#endif
