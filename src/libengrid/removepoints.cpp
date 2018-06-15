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
#include "removepoints.h"

#include "checksurfaceintegrity.h"

#include "geometrytools.h"
using namespace GeometryTools;

#include <cmath>
using namespace std;

#include <iostream>
using namespace std;

RemovePoints::RemovePoints() : SurfaceOperation() {
  setQuickSave(true);
  getSet("surface meshing", "point removal threshold", 2, m_Threshold);
  m_ProtectFeatureEdges = false;
  m_PerformGeometricChecks = true;
  m_UpdatePSP = false;
}

void RemovePoints::markFeatureEdges()
{
  EG_VTKDCN(vtkDoubleArray, characteristic_length_desired, m_Grid, "node_meshdensity_desired");
  m_IsFeatureNode.fill(false, m_Part.getNumberOfNodes());
  if(m_ProtectFeatureEdges) {
    for (int i_nodes = 0; i_nodes < m_Part.getNumberOfNodes(); ++i_nodes) {
      if (!m_IsFeatureNode[i_nodes]) {
        vtkIdType id_node1 = m_Part.globalNode(i_nodes);
        for (int j = 0; j < m_Part.n2nLSize(i_nodes); ++j) {
          vtkIdType id_node2 = m_Part.n2nLG(i_nodes, j);
          QSet<vtkIdType> edge_cells;
          int N = getEdgeCells(id_node1, id_node2, edge_cells);
          if (N != 2) {
            m_IsFeatureNode[i_nodes] = true;
            m_IsFeatureNode[m_Part.localNode(id_node2)] = true;
          } else {
            QSet<vtkIdType>::iterator iter = edge_cells.begin();
            vtkIdType id_cell1 = *iter;
            ++iter;
            vtkIdType id_cell2 = *iter;
            vec3_t n1 = cellNormal(m_Grid, id_cell1);
            vec3_t n2 = cellNormal(m_Grid, id_cell2);
            if (angle(n1, n2) >= m_FeatureAngle) {
              m_IsFeatureNode[i_nodes] = true;
              m_IsFeatureNode[m_Part.localNode(id_node2)] = true;
            }
          }
        }
      }
    }
    int N = 0;
    for (int i = 0; i < m_IsFeatureNode.size(); ++i) {
      if(m_IsFeatureNode[i]) {
        ++N;
      }
    }
    cout << N << " nodes on feature edges (angle >= " << GeometryTools::rad2deg(m_FeatureAngle) << "deg)" << endl;
  }
}

void RemovePoints::fixNodes(const QVector<bool> &fixnodes)
{
  if (fixnodes.size() != m_Grid->GetNumberOfPoints()) {
    EG_BUG;
  }
  m_Fixed = fixnodes;
}

bool RemovePoints::checkEdge(vtkIdType id_node1, vtkIdType id_node2)
{
  EG_VTKDCN(vtkDoubleArray, cl, m_Grid, "node_meshdensity_desired");
  double cl1 = cl->GetValue(id_node1);
  double cl2 = cl->GetValue(id_node2);
  vec3_t x1, x2;
  m_Grid->GetPoint(id_node1, x1.data());
  m_Grid->GetPoint(id_node2, x2.data());
  double L = (x1 - x2).abs();
  double cl_crit = max(cl1, cl2) / m_Threshold;
  if (L < cl_crit) {
    return true;
  }
  return false;
}

void RemovePoints::operate()
{
  if (m_Fixed.size() != m_Grid->GetNumberOfPoints()) {
    m_Fixed.fill(false, m_Grid->GetNumberOfPoints());
  }

  MeshPartition full_partition(m_Grid, true);
  full_partition.setAllCells();

  int N1 = m_Grid->GetNumberOfPoints();

  QVector<vtkIdType> selected_cells;
  getSurfaceCells(m_BoundaryCodes, selected_cells, m_Grid);
  QVector<vtkIdType> selected_nodes;
  getNodesFromCells(selected_cells, selected_nodes, m_Grid);

  setAllSurfaceCells();
  computeNormals();
  updateNodeInfo();

  EG_VTKDCN(vtkCharArray, node_type, m_Grid, "node_type");
  EG_VTKDCC(vtkIntArray, cell_code, m_Grid, "cell_code");
  EG_VTKDCN(vtkDoubleArray, characteristic_length_desired, m_Grid, "node_meshdensity_desired");

  // global values
  QVector <vtkIdType> all_deadcells;
  QVector <vtkIdType> all_mutatedcells;
  int num_newpoints = 0;
  int num_newcells = 0;

  QVector<bool> marked_nodes(m_Grid->GetNumberOfPoints(), false);

  QVector <vtkIdType> deadnode_vector;
  QVector <vtkIdType> snappoint_vector;

  foreach (vtkIdType id_node, selected_nodes) {
    if (node_type->GetValue(id_node) != EG_FIXED_VERTEX && !m_Fixed[id_node]) {

      // check that node is only surrounded by triangles
      bool tri_only = true;
      for (int j = 0; j < m_Part.n2cGSize(id_node); ++j) {
        vtkIdType id_cell = m_Part.n2cGG(id_node, j);
        if(m_Grid->GetCellType(id_cell) != VTK_TRIANGLE) {
          tri_only = false;
          break;
        }
      }

      if (!marked_nodes[id_node] && tri_only) {
        for (int j = 0; j < m_Part.n2nGSize(id_node); ++j) {
          vtkIdType id_neigh = m_Part.n2nGG(id_node, j);
          if (checkEdge(id_node, id_neigh)) {
            QVector <vtkIdType> dead_cells;
            QVector <vtkIdType> mutated_cells;
            int l_num_newpoints = 0;
            int l_num_newcells = 0;
            if (isSnapPoint(id_node, id_neigh, dead_cells, mutated_cells, l_num_newpoints, l_num_newcells, marked_nodes)) {

              // add deadnode/snappoint pair
              deadnode_vector.push_back(id_node);
              snappoint_vector.push_back(id_neigh);

              // update global values
              num_newpoints += l_num_newpoints;
              num_newcells  += l_num_newcells;
              all_deadcells += dead_cells;
              all_mutatedcells += mutated_cells;

              // mark neighbour nodes
              for (int k = 0; k < m_Part.n2nGSize(id_node); ++k) {
                marked_nodes[m_Part.n2nGG(id_node, k)] = true;
              }
              for (int k = 0; k < m_Part.n2nGSize(id_neigh); ++k) {
                marked_nodes[m_Part.n2nGG(id_neigh, k)] = true;
              }

              break;
            }
          }
        }
      }
    }
  }

  //delete
  if (num_newpoints != -deadnode_vector.size()) {
    EG_BUG;
  }
  if (num_newcells != -all_deadcells.size()) {
    EG_BUG;
  }

  deleteSetOfPoints(deadnode_vector, snappoint_vector, all_deadcells, all_mutatedcells);

  int N2 = m_Grid->GetNumberOfPoints();
  m_NumRemoved = N1 - N2;
}

/// \todo finish this function and optimize it.
bool RemovePoints::checkForDestroyedVolumes(vtkIdType id_node1, vtkIdType id_node2, int& N_common_points)
{
  if (id_node1 == id_node2) {
    EG_BUG;
  }

  l2l_t  n2n   = getPartN2N();
  g2l_t _nodes = getPartLocalNodes();
  l2g_t nodes  = getPartNodes();

  QVector<int> node1_neighbours = n2n[_nodes[id_node1]];
  QVector<int> node2_neighbours = n2n[_nodes[id_node2]];
  QVector<int> common_points;
  qcontIntersection(node1_neighbours, node2_neighbours, common_points);
  // set N_common_points
  N_common_points = common_points.size();

  // TEST 0: TOPOLOGICAL: DeadNode, PSP and any common point must belong to a cell.
  for(int i = 0; i < N_common_points; i++) {
    int i_common_point_1 = common_points[i];
    vtkIdType id_common_point_1 = nodes[i_common_point_1];
    if(!isCell(id_node1, id_node2, id_common_point_1)) {
      if(DebugLevel > 100) {
        qDebug() << "test 0 failed";
        qDebug() << "id_node1, id_node2, id_common_point_1=" << id_node1 << id_node2 << id_common_point_1;
      }
      return true;
    }
    // TEST 1: TOPOLOGICAL: Moving DeadNode to PSP must not lay any cell on another cell.
    //                      => For any pair of common points (cp1,cp2), (cp1,cp2,DeadNode)+(cp1,cp2,PSP)
    //                         must not be cells at the same time!
    for(int j = i + 1; j < common_points.size(); j++) {
      int i_common_point_2 = common_points[j];
      vtkIdType id_common_point_2 = nodes[i_common_point_2];
      if(isCell(id_common_point_1, id_common_point_2, id_node1) && isCell(id_common_point_1, id_common_point_2, id_node2)) {
        if(DebugLevel > 100) {
          qDebug() << "test 1 failed";
          qDebug() << "id_common_point_1, id_common_point_2, id_node1=" << id_common_point_1 << id_common_point_2 << id_node1;
          qDebug() << "id_common_point_1, id_common_point_2, id_node2=" << id_common_point_1 << id_common_point_2 << id_node2;
        }
        return true;
      }
    }
  }

  QSet<vtkIdType> all_faces;
  for (int i = 0; i < m_Part.n2nGSize(id_node1); ++i) {
    for (int j = 0; j < m_Part.n2cGSize(m_Part.n2nGG(id_node1, i)); ++j) {
      all_faces.insert(m_Part.n2cGG(m_Part.n2nGG(id_node1, i), j));
    }
  }
  for (int i = 0; i < m_Part.n2nGSize(id_node2); ++i) {
    for (int j = 0; j < m_Part.n2cGSize(m_Part.n2nGG(id_node2, i)); ++j) {
      all_faces.insert(m_Part.n2cGG(m_Part.n2nGG(id_node2, i), j));
    }
  }
  QSet<vtkIdType> near_faces;
  for (int i = 0; i < m_Part.n2cGSize(id_node1); ++i) {
    near_faces.insert(m_Part.n2cGG(id_node1, i));
  }
  for (int i = 0; i < m_Part.n2cGSize(id_node2); ++i) {
    near_faces.insert(m_Part.n2cGG(id_node2, i));
  }
  QSet<vtkIdType> far_faces = all_faces - near_faces;
  bool tetra = true;
  foreach (vtkIdType id_cell, far_faces) {
    vtkIdType N_pts, *pts;
    m_Grid->GetCellPoints(id_cell, N_pts, pts);
    for (int i = 0; i < N_pts; ++i) {
      if (!m_Part.hasNeighNode(pts[i], id_node1) && !m_Part.hasNeighNode(pts[i], id_node2)) {
        tetra = false;
        break;
      }
    }
    if (!tetra) {
      break;
    }
  }
  if (tetra) {
    return true;
  }



  //FIX THIS!!!!
  /*
  // check if DeadNode, PSP and common points form a tetrahedron.
  if ( n2n[_nodes[intersection1]].contains( _nodes[intersection2] ) ) { //if there's an edge between intersection1 and intersection2
    //check if (node1,intersection1,intersection2) and (node2,intersection1,intersection2) are defined as cells!
    QVector<int> S1 = n2c[_nodes[intersection1]];
    QVector<int> S2 = n2c[_nodes[intersection2]];
    QVector<int> Si;
    qcontIntersection( S1, S2, Si );
    int counter = 0;
    foreach( int i_cell, Si ) {
      vtkIdType num_pts, *pts;
      m_Grid->GetCellPoints( cells[i_cell], num_pts, pts );
      for ( int i = 0; i < num_pts; ++i ) {
        if ( pts[i] == id_node1 || pts[i] == id_node2 ) counter++;
      }
    }
    if ( counter >= 2 ) {
      IsTetra = true;
    }
  }
*/
  return false;
}

int RemovePoints::numberOfCommonPoints(vtkIdType id_node1, vtkIdType id_node2, bool& IsTetra)
{
  l2l_t  n2n   = getPartN2N();
  l2l_t  n2c   = getPartN2C();
  g2l_t _nodes = getPartLocalNodes();
  l2g_t nodes  = getPartNodes();
  l2g_t cells = getPartCells();

  QVector<int> node1_neighbours = n2n[_nodes[id_node1]];
  QVector<int> node2_neighbours = n2n[_nodes[id_node2]];
  QVector<int> intersection;
  qcontIntersection(node1_neighbours, node2_neighbours, intersection);
  int N = intersection.size();
  IsTetra = false;
  if(N == 2) {
    vtkIdType intersection1 = nodes[intersection[0]];
    vtkIdType intersection2 = nodes[intersection[1]];

    // test if id_node1, id_node2 and intersection* form a cell
    QVector <vtkIdType> EdgeCells_1i;
    QVector <vtkIdType> EdgeCells_2i;
    QVector <vtkIdType> inter;
    int Ncells;

    // intersection1
    Ncells = getEdgeCells(id_node1, intersection1, EdgeCells_1i);
    if(N != 2) {
      qWarning() << "Ncells=" << Ncells;
      EG_BUG;
    }
    Ncells = getEdgeCells(id_node2, intersection1, EdgeCells_2i);
    if(Ncells != 2) {
      qWarning() << "Ncells=" << Ncells;
      EG_BUG;
    }
    qcontIntersection(EdgeCells_1i, EdgeCells_2i, inter);
    if(inter.size() <= 0) EG_BUG;   // (id_node1, id_node2, intersection1) is not a cell

    // intersection2
    Ncells = getEdgeCells(id_node1, intersection2, EdgeCells_1i);
    if(Ncells != 2) {
      qWarning() << "Ncells=" << Ncells;
      EG_BUG;
    }
    Ncells = getEdgeCells(id_node2, intersection2, EdgeCells_2i);
    if(Ncells != 2) {
      qWarning() << "Ncells=" << Ncells;
      EG_BUG;
    }
    qcontIntersection(EdgeCells_1i, EdgeCells_2i, inter);
    if(inter.size() <= 0) EG_BUG;   // (id_node1, id_node2, intersection2) is not a cell

    // check if DeadNode, PSP and common points form a tetrahedron.
    if(n2n[_nodes[intersection1]].contains(_nodes[intersection2])) {      //if there's an edge between intersection1 and intersection2
      //check if (node1,intersection1,intersection2) and (node2,intersection1,intersection2) are defined as cells!
      QVector<int> S1 = n2c[_nodes[intersection1]];
      QVector<int> S2 = n2c[_nodes[intersection2]];
      QVector<int> Si;
      qcontIntersection(S1, S2, Si);
      int counter = 0;
      foreach(int i_cell, Si) {
        vtkIdType num_pts, *pts;
        m_Grid->GetCellPoints(cells[i_cell], num_pts, pts);
        for(int i = 0; i < num_pts; ++i) {
          if(pts[i] == id_node1 || pts[i] == id_node2) counter++;
        }
      }
      if(counter >= 2) {
        IsTetra = true;
      }
    }
  }
  return(N);
}

bool RemovePoints::flippedCell2(vtkIdType id_node, vec3_t x_new) {
  /*
  for (int i = 0; i < m_Part.n2cGSize(id_node); ++i) {

    vtkIdType id_cell = m_Part.n2cGG(id_node, i);

    vtkIdType N_pts, *pts;
    m_Grid->GetCellPoints(id_cell, N_pts, pts);
    if(N_pts!=3) EG_BUG;

    int i_pts=0;
    for(i_pts=0; i_pts<N_pts; i_pts++) {
      if(pts[i_pts]==id_node) break;
    }
    if(pts[i_pts]!=id_node) EG_BUG;

    vec3_t x1, x2, x_old;
    m_Grid->GetPoint(pts[(i_pts+1)%N_pts],x1.data());
    m_Grid->GetPoint(pts[(i_pts+2)%N_pts],x2.data());

    vec3_t old_cell_normal = GeometryTools::triNormal(x_old, x1, x2);
    vec3_t new_cell_normal = GeometryTools::triNormal(x_new, x1, x2);

    if(old_cell_normal.abs2()==0) EG_BUG;
    if(old_cell_normal.abs2()==0) EG_BUG;

    GeometryTools::cellNormal(m_Grid, );
    cell_normals.normalise();

    vtkIdType *pts;
    vtkIdType npts;
    vec3_t n(0,0,0);
    grid->GetCellPoints(i, npts, pts);
    if (npts == 3) {
      return triNormal(grid,pts[0],pts[1],pts[2]);
    } else if (npts == 4) {
      return quadNormal(grid,pts[0],pts[1],pts[2],pts[3]);
    } else {
      EG_BUG;
    }
    return n;

  }
  */
  return true;
}

/// \todo adapt for multiple volumes
bool RemovePoints::flippedCell(vtkIdType id_node, vec3_t x_new, vtkIdType id_cell) {
  if( m_Grid->GetCellType(id_cell) == VTK_WEDGE ) EG_BUG;

  vec3_t x_old;
  m_Grid->GetPoint(id_node, x_old.data());

  vec3_t n(0, 0, 0);
  bool move = true;
  QVector<vec3_t> cell_normals(m_Part.n2cGSize(id_node));
  double A_max = 0;
  for(int i = 0; i < m_Part.n2cGSize(id_node); ++i) {
    double A = fabs(GeometryTools::cellVA(m_Grid, m_Part.n2cGG(id_node, i)));
    A_max = max(A, A_max);
    cell_normals[i] = GeometryTools::cellNormal(m_Grid, m_Part.n2cGG(id_node, i));
    cell_normals[i].normalise();
  }
  int N = 0;
  for(int i = 0; i < m_Part.n2cGSize(id_node); ++i) {
    double A = fabs(GeometryTools::cellVA(m_Grid, m_Part.n2cGG(id_node, i)));
    if(A > 0.01 * A_max) {
      n += cell_normals[i];
      ++N;
    }
  }
  if(N == 0) {
    move = false;
  } else {
    n.normalise();
    double L_max = 0;
    for(int i = 0; i < m_Part.n2nGSize(id_node); ++i) {
      vec3_t xn;
      m_Grid->GetPoint(m_Part.n2nGG(id_node, i), xn.data());
      double L = (xn - x_old).abs();
      L_max = max(L, L_max);
    }
    vec3_t x_summit = x_old + L_max * n;
    vec3_t x[3];
    vtkIdType N_pts, *pts;
    m_Grid->GetCellPoints(id_cell, N_pts, pts);
    if(N_pts != 3) {
      EG_BUG;
    }
    for(int j = 0; j < N_pts; ++j) {
      m_Grid->GetPoint(pts[j], x[j].data());
    }
    if(GeometryTools::tetraVol(x[0], x[1], x[2], x_summit, false) <= 0) {
      move = false;
    }
  }

  return !move;
}

/** This function tries to find a valid snappoint for DeadNode and returns its ID if it finds one, otherwise it returns -1.
    If a valid snappoint is found, the corresponding dead and mutated cells are returned via DeadCells and MutatedCells.

 DEFINITIONS:
 Normal cell: nothing has changed
 Dead cell: the cell does not exist anymore
 Mutated cell: the cell's form has changed

 Basic algorithm:\n
 foreach(potential snap point of DeadNode) {\n
   bool IsValidSnapPoint = true;\n
   some tests; if any fails: IsValidSnapPoint = false; continue;\n
   // reset output variables\n
   num_newpoints = -1;\n
   num_newcells = 0;\n
   DeadCells.clear();\n
   MutatedCells.clear();\n
   foreach(neighbour cell of DeadNode) {\n
     more tests; if any fails: IsValidSnapPoint = false; continue;\n
     fill DeadCells + MutatedCells;\n
   }\n
   even more tests; if any fails: IsValidSnapPoint = false; continue;\n
   if(IsValidSnapPoint) {\n
     SnapPoint = PSP;\n
     break;\n
   }\n
 }\n

 \todo Clean up this function
 */
bool RemovePoints::isSnapPoint(vtkIdType id_node1, vtkIdType id_node2,
                               QVector<vtkIdType>& dead_cells, QVector<vtkIdType>& mutated_cells,
                               int& num_newpoints, int& num_newcells,
                               const QVector<bool>& marked_nodes)
{
  EG_VTKDCN(vtkCharArray, node_type, m_Grid, "node_type");
  if (node_type->GetValue(id_node1) == EG_FIXED_VERTEX) {
    EG_BUG;
  }

  QVector <vtkIdType> psp_vector = getPotentialSnapPoints(id_node1);
  if (!psp_vector.contains(id_node2)) {
    return false;
  }

  vec3_t x1;
  m_Grid->GetPoint(id_node1, x1.data());

  // TEST 1 : TOPOLOGICAL : Is the node already marked?
  if (marked_nodes[id_node2]) {
    return false;
  }

  // TEST 2: TOPOLOGICAL: do not cut off feature corners
  {
    QSet<vtkIdType> common_nodes, n2n2;
    m_Part.getGlobalN2N(id_node1, common_nodes);
    m_Part.getGlobalN2N(id_node2, n2n2);
    common_nodes.intersect(n2n2);
    foreach (vtkIdType id_neigh, common_nodes) {
      if (node_type->GetValue(id_neigh) == EG_FEATURE_CORNER_VERTEX) {
        return false;
      }
    }
  }

  // TEST 3: TOPOLOGICAL: Moving id_node1 to id_node2 must not lay any cell on another cell.
  int num_common_points = 0;
  if (checkForDestroyedVolumes(id_node1, id_node2, num_common_points)) {
    return false;
  }

  // TEST 4: normal irregularity
  if (normalIrregularity(id_node1) > normalIrregularity(id_node2)) {
    //return false;
  }

  //count number of points and cells to remove + analyse cell transformations
  num_newpoints = -1;
  num_newcells = 0;
  dead_cells.clear();
  mutated_cells.clear();
  //foreach (int i_cell, n2c[_nodes[DeadNode]]) {   //loop through potentially dead cells
  for (int i = 0; i < m_Part.n2cGSize(id_node1); ++i) {
    vtkIdType id_cell = m_Part.n2cGG(id_node1, i);

    EG_GET_CELL_AND_TYPE(id_cell, m_Grid);
    if (type_cell == VTK_WEDGE) {
      EG_BUG;
    }
    if(num_pts != 3) {
      return false;
    }

    bool contains_snap_point = false;
    bool invincible = false; // a point with only one cell is declared invincible.
    for (int j = 0; j < num_pts; ++j) {
      if (pts[j] == id_node2) {
        contains_snap_point = true;
      }
      if (pts[j] != id_node1 && pts[j] != id_node2 && m_Part.n2cGSize(pts[j]) <= 1) {
        invincible = true;
      }
    }

    if (contains_snap_point) {    // potential dead cell
      if (invincible) {
        // TEST 4: TOPOLOGICAL: Check that empty lines aren't left behind when a cell is killed
        return false;
      } else {
        dead_cells.push_back(id_cell);
        num_newcells -= 1;
      }
    } else { // if the cell does not contain the SnapPoint (potential mutated cell)

      vtkIdType old_triangle[3];
      vtkIdType new_triangle[3];

      for (int j = 0; j < num_pts; ++j) {
        old_triangle[j] = pts[j];
        new_triangle[j] = ((pts[j] == id_node1) ? id_node2 : pts[j]);
      }
      vec3_t n_old = triNormal(m_Grid, old_triangle[0], old_triangle[1], old_triangle[2]);
      vec3_t n_new = triNormal(m_Grid, new_triangle[0], new_triangle[1], new_triangle[2]);
      double A_old = n_old.abs();
      double A_new = n_new.abs();
      n_old.normalise();
      n_new.normalise();

      // TEST 5: GEOMETRICAL: area + inversion check
      if (m_PerformGeometricChecks) {
        //if(Old_N * New_N < 0.1 || New_N * New_N < Old_N * Old_N * 1. / 100.) {
        if (n_old * n_new < 0.2 || A_new < 0.1*A_old) {
          return false;
        }
      }

      mutated_cells.push_back(id_cell);
    } // end of if the cell does not contain the SnapPoint (potential mutated cell)
  }

  // TEST 6: TOPOLOGICAL: survivor check
  if (m_Grid->GetNumberOfCells() + num_newcells <= 0) {
    return false;
  }

  return true;
}

bool RemovePoints::deleteSetOfPoints(const QVector<vtkIdType>& deadnode_vector,
                                     const QVector<vtkIdType>& snappoint_vector,
                                     const QVector<vtkIdType>& all_deadcells,
                                     const QVector<vtkIdType>& all_mutatedcells) {
  int initial_num_points = m_Grid->GetNumberOfPoints();

  //src grid info
  int num_points = m_Grid->GetNumberOfPoints();
  int num_cells = m_Grid->GetNumberOfCells();

  int num_newcells = -all_deadcells.size();
  int num_newpoints = -deadnode_vector.size();

  //  if ( num_newcells != 2*num_newpoints ) {
  //    EG_BUG;
  //  }

  //allocate
  EG_VTKSP(vtkUnstructuredGrid, dst);
  allocateGrid(dst, num_cells + num_newcells, num_points + num_newpoints);

  //vector used to redefine the new point IDs
  QVector <vtkIdType> OffSet(num_points);

  //copy undead points
  QVector<bool> is_deadnode(m_Grid->GetNumberOfPoints(), false);
  QVector<int> glob2dead(m_Grid->GetNumberOfPoints(), -1);
  for(int i_deadnodes = 0; i_deadnodes < deadnode_vector.size(); ++i_deadnodes) {
    vtkIdType id_node = deadnode_vector[i_deadnodes];
    if(id_node > m_Grid->GetNumberOfPoints()) {
      EG_BUG;
    }
    is_deadnode[id_node] = true;
    glob2dead[id_node] = i_deadnodes;
  }
  vtkIdType dst_id_node = 0;
  for(vtkIdType src_id_node = 0; src_id_node < num_points; ++src_id_node) { //loop through src points
    OffSet[src_id_node] = src_id_node - dst_id_node;
    if(!is_deadnode[src_id_node]) {  //if the node isn't dead, copy it
      vec3_t x;
      m_Grid->GetPoints()->GetPoint(src_id_node, x.data());
      dst->GetPoints()->SetPoint(dst_id_node, x.data());
      copyNodeData(m_Grid, src_id_node, dst, dst_id_node);
      dst_id_node++;
    } else {
      if(DebugLevel > 0) {
        cout << "dead node encountered: src_id_node=" << src_id_node << " dst_id_node=" << dst_id_node << endl;
      }
    }
  }

  //Copy undead cells

  // Fill is_deadcell
  QVector<bool> is_deadcell(m_Grid->GetNumberOfCells(), false);
  foreach(vtkIdType id_cell, all_deadcells) {
    if( m_Grid->GetCellType(id_cell) == VTK_WEDGE ) EG_BUG;
    is_deadcell[id_cell] = true;
  }

  // Fill is_mutatedcell
  QVector<bool> is_mutatedcell(m_Grid->GetNumberOfCells(), false);
  foreach(vtkIdType id_cell, all_mutatedcells) {
    if( m_Grid->GetCellType(id_cell) == VTK_WEDGE ) EG_BUG;
    is_mutatedcell[id_cell] = true;
  }

  for(vtkIdType id_cell = 0; id_cell < m_Grid->GetNumberOfCells(); ++id_cell) { //loop through src cells
    //    if( m_Grid->GetCellType(id_cell) == VTK_WEDGE ) continue;
    //    if(isVolume(id_cell, m_Grid)) continue;

    if(!is_deadcell[id_cell]) {  //if the cell isn't dead
      vtkIdType src_num_pts, *src_pts;
      m_Grid->GetCellPoints(id_cell, src_num_pts, src_pts);
      vtkIdType type_cell = m_Grid->GetCellType(id_cell);

      vtkIdType dst_num_pts = src_num_pts;
      QVector<vtkIdType> dst_pts(dst_num_pts);

      if(is_mutatedcell[id_cell]) {  //mutated cell
        if(dst_num_pts != 3) {
          // Not fully supported yet
          qWarning() << "all_mutatedcells=" << all_mutatedcells;
          qWarning() << "A non-triangle cell was mutated!";
          EG_BUG;
        }
        int num_deadnode = 0;
        for(int i = 0; i < src_num_pts; i++) {
          int DeadIndex = glob2dead[src_pts[i]];
          if(DeadIndex != -1) { // It is a dead node.
            dst_pts[i] = snappoint_vector[DeadIndex] - OffSet[snappoint_vector[DeadIndex]]; // dead node
            num_deadnode++;
          } else {
            dst_pts[i] = src_pts[i] - OffSet[src_pts[i]]; // not a dead node
          }
        }
        if(num_deadnode != 1) {
          qWarning() << "FATAL ERROR: Mutated cell has more than one dead node!";
          qWarning() << "num_deadnode=" << num_deadnode;
          qWarning() << "type_cell=" << type_cell << " VTK_TRIANGLE=" << VTK_TRIANGLE << " VTK_QUAD=" << VTK_QUAD;
          EG_BUG;
        }
      } else { //normal cell

        if(DebugLevel > 10) {
          cout << "processing normal cell " << id_cell << endl;
        }

        if(isVolume(id_cell, m_Grid)) {
          int num_deadnode = 0;
          for(int i = 0; i < src_num_pts; i++) {
            int DeadIndex = glob2dead[src_pts[i]];
            if(DeadIndex != -1) { // It is a dead node.
              dst_pts[i] = snappoint_vector[DeadIndex] - OffSet[snappoint_vector[DeadIndex]]; // dead node
              num_deadnode++;
            } else {
              dst_pts[i] = src_pts[i] - OffSet[src_pts[i]]; // not a dead node
            }
          }
          if(num_deadnode > 1) {
            qWarning() << "FATAL ERROR: Mutated cell has more than one dead node!";
            qWarning() << "num_deadnode=" << num_deadnode;
            qWarning() << "type_cell=" << type_cell << " VTK_TRIANGLE=" << VTK_TRIANGLE << " VTK_QUAD=" << VTK_QUAD;
            for(int k = 0; k < src_num_pts; k++) {
              int DeadIndex = glob2dead[src_pts[k]];
              qWarning()<<"k="<<k<<" DeadIndex="<<"glob2dead["<<src_pts[k]<<"]="<<DeadIndex;
            }
            EG_BUG;
          }
        }
        else {
          for(int j = 0; j < src_num_pts; j++) {
            if(is_deadnode[src_pts[j]]) {
              qWarning() << "FATAL ERROR: Normal cell contains a dead node!";
              qWarning() << "is_deadnode=" << is_deadnode;
              qWarning() << "src_pts[" << j << "]=" << src_pts[j];
              qWarning() << "type_cell=" << type_cell << " VTK_TRIANGLE=" << VTK_TRIANGLE << " VTK_QUAD=" << VTK_QUAD;
              saveGrid(m_Grid, "crash");
              EG_BUG;
            }
            dst_pts[j] = src_pts[j] - OffSet[src_pts[j]];
          }
        }

      } // end of normal cell processing

      // copy the cell
      //\todo adapt type_cell in the case of mutilated cells!
      vtkIdType id_new_cell = dst->InsertNextCell(type_cell, dst_num_pts, dst_pts.data());
      copyCellData(m_Grid, id_cell, dst, id_new_cell);
    } //end of undead cell processing
  } //end of loop through src cells

  // update m_Grid
  makeCopy(dst, m_Grid);

  int final_num_points = m_Grid->GetNumberOfPoints();
  if(initial_num_points - final_num_points != deadnode_vector.size()) {
    EG_BUG;
  }

  return(true);
}
//End of DeleteSetOfPoints
