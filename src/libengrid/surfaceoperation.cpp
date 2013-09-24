// 
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008-2013 enGits GmbH                                      +
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
#include "surfaceoperation.h"

#include "guimainwindow.h"

#include <vtkCharArray.h>
#include <vtkMath.h>
#include <vtkCellArray.h>
#include <vtkPolygon.h>

#include "geometrytools.h"
#include "meshqualityfaceorientation.h"

using namespace GeometryTools;

SurfaceOperation::SurfaceOperation() : Operation()
{
  //default values for determining node types and for smoothing operations
  getSet("surface meshing", "edge angle to determine fixed vertices", 180, m_EdgeAngle);
  getSet("surface meshing", "feature angle", 20, m_FeatureAngle);
  getSet("surface meshing", "boundary codes define features", true, m_BCodeFeatureDefinition);
  getSet("surface meshing", "number of steps to protect node types", 2, m_TypeProtectionCount);
  getSet("surface meshing", "threshold for face orientation quality", 0.1, m_FaceOrientationThreshold);
  m_FeatureAngle = GeometryTools::deg2rad(m_FeatureAngle);
  m_EdgeAngle = GeometryTools::deg2rad(m_EdgeAngle);
  setEdgeAngle(m_EdgeAngle);
  m_StretchingFactor = 0;
  m_UniformSnapPoints = false;
  m_StrictFeatureSnap = !m_BCodeFeatureDefinition;
}

void SurfaceOperation::operate()
{

}

ostream& operator<<(ostream &out, stencil_t S)
{
  out << "S.id_cell = " << S.id_cell << " ";
  out << "S.id_node = " << S.id_node << " ";
  out << "S.sameBC = " << S.sameBC << " ";
  out << "S.type = " << S.type_cell << " ";
  out << "S.p1 = " << S.p1 << " ";
  out << "S.p2 = " << S.p2 << " ";
  return(out);
}

stencil_t SurfaceOperation::getStencil(vtkIdType id_cell1, int j1)
{
  stencil_t S;
  {
    vtkIdType N_pts, *pts;
    m_Grid->GetCellPoints(id_cell1, N_pts, pts);
    S.p1 = pts[j1];
    S.p2 = pts[0];
    if (j1 < N_pts - 1) {
      S.p2 = pts[j1 + 1];
    }
  }
  QSet<vtkIdType> cells_p1;
  for (int i = 0; i < m_Part.n2cGSize(S.p1); ++i) {
    vtkIdType id_cell = m_Part.n2cGG(S.p1, i);
    if (id_cell != id_cell1) {
      cells_p1.insert(id_cell);
    }
  }
  QSet<vtkIdType> cells_p2;
  for (int i = 0; i < m_Part.n2cGSize(S.p2); ++i) {
    vtkIdType id_cell = m_Part.n2cGG(S.p2, i);
    if (id_cell != id_cell1) {
      cells_p2.insert(id_cell);
    }
  }
  QSet<vtkIdType> cells = cells_p1.intersect(cells_p2);
  EG_VTKDCC(vtkIntArray, cell_code, m_Grid, "cell_code");
  S.sameBC = true;
  S.id_cell.resize(1);
  S.id_cell[0] = id_cell1;
  foreach (vtkIdType id_cell, cells) {
    if (isSurface(id_cell, m_Grid)) {
      S.id_cell.push_back(id_cell);
      if (cell_code->GetValue(id_cell) != cell_code->GetValue(id_cell1)) {
        S.sameBC = false;
      }
    }
  }
  S.id_node.resize(S.id_cell.size());
  S.type_cell.resize(S.id_cell.size());
  for (int i = 0; i < S.id_cell.size(); ++i) {
    vtkIdType N_pts, *pts;
    m_Grid->GetCellPoints(S.id_cell[i], N_pts, pts);
    S.type_cell[i] = m_Grid->GetCellType(S.id_cell[i]);
    for (int j = 0; j < N_pts; ++j) {
      if (pts[j] != S.p1 && pts[j] != S.p2) {
        S.id_node[i] = pts[j];
        break;
      }
    }
  }
  return S;
}

void SurfaceOperation::setBCodesFeatureDefinition(bool flag)
{
  m_BCodeFeatureDefinition = flag;
  if (m_BCodeFeatureDefinition) {
    m_FeatureAngle = deg2rad(180);
  }
}

int SurfaceOperation::UpdateCurrentMeshDensity()
{
  if ( DebugLevel > 0 ) {
    cout << "===UpdateMeshDensity START===" << endl;
  }
  QVector<vtkIdType> cells;
  getAllSurfaceCells( cells, m_Grid );
  EG_VTKDCC( vtkIntArray, cell_code, m_Grid, "cell_code" );
  EG_VTKDCN( vtkDoubleArray, characteristic_length_desired, m_Grid, "node_meshdensity_desired" );
  setGrid( m_Grid );
  setCells( cells );
  if ( DebugLevel > 5 ) {
    cout << "cells.size()=" << cells.size() << endl;
  }
  EG_VTKDCN( vtkDoubleArray, node_meshdensity_current, m_Grid, "node_meshdensity_current" );
  l2g_t nodes = getPartNodes();
  foreach( vtkIdType node, nodes ) {
    node_meshdensity_current->SetValue( node, CurrentMeshDensity( node ) );
  }
  if ( DebugLevel > 0 ) {
    cout << "===UpdateMeshDensity END===" << endl;
  }
  return( 0 ); ///\todo what for???
}

void SurfaceOperation::readVMD()
{
  QString buffer = GuiMainWindow::pointer()->getXmlSection("engrid/surface/table").replace("\n", " ");
  int row_count = 0;
  int column_count = 0;
  m_VMDvector.clear();

  if(!buffer.isEmpty()) {
    QTextStream in(&buffer, QIODevice::ReadOnly);
    in >> row_count >> column_count;
    QVector<int> tmp_bcs;
    GuiMainWindow::pointer()->getAllBoundaryCodes(tmp_bcs);
    if (column_count == tmp_bcs.size() + 3) {
      m_VMDvector.fill(VertexMeshDensity(), row_count);
      for (int i = 0; i < row_count; ++i) {
        int row, column;
        QString formula;
        foreach (int bc, tmp_bcs) {
          in >> row >> column >> formula;
          m_VMDvector[row].BCmap[bc] = formula.toInt();
        }
        in >> row >> column >> formula;
        m_VMDvector[row].type = Str2VertexType(formula);
        in >> row >> column >> formula;
        if (formula == "{{{empty}}}") {
          formula = "";
        }
        m_VMDvector[i].setNodes(formula);
        in >> row >> column >> formula;
        m_VMDvector[i].density = formula.toDouble();
      }
    } else {
      EG_ERR_RETURN(QObject::tr("Mismatch of number of boundary codes!"));
    }
  }
}

void SurfaceOperation::updateNodeInfo()
{
  setAllCells();
  readVMD();
  l2g_t nodes = getPartNodes();
  computeNormals();
  EG_VTKDCN(vtkCharArray, node_type, m_Grid, "node_type");
  EG_VTKDCN(vtkIntArray, node_type_counter, m_Grid, "node_type_counter");
  EG_VTKDCC(vtkIntArray, cell_code, m_Grid, "cell_code");
  foreach (vtkIdType id_node, nodes) {

    char old_type = node_type->GetValue(id_node);
    char new_type = getNodeType(id_node, true);

    if (old_type == EG_FIXED_VERTEX && new_type != EG_FIXED_VERTEX) {
      //EG_BUG;
    }

    if ((old_type != EG_FEATURE_CORNER_VERTEX && old_type != EG_FEATURE_EDGE_VERTEX) || m_BCodeFeatureDefinition) {
      node_type->SetValue(id_node, new_type);
    } else {

      bool done = false;

      // EG_FEATURE_CORNER_VERTEX -> any other type
      if (old_type == EG_FEATURE_CORNER_VERTEX && new_type != EG_FEATURE_CORNER_VERTEX) {
        if (node_type_counter->GetValue(id_node) > m_TypeProtectionCount) {
          node_type->SetValue(id_node, new_type);
        } else {
          node_type_counter->SetValue(id_node, node_type_counter->GetValue(id_node) + 1);
        }
        done = true;
      }

      // EG_FEATURE_EDGE_VERTEX -> EG_SIMPLE_VERTEX
      if (old_type == EG_FEATURE_EDGE_VERTEX && new_type == EG_SIMPLE_VERTEX) {
        if (node_type_counter->GetValue(id_node) > m_TypeProtectionCount) {
          node_type->SetValue(id_node, new_type);
        } else {
          node_type_counter->SetValue(id_node, node_type_counter->GetValue(id_node) + 1);
        }
        done = true;
      }

      if (!done) {
        node_type->SetValue(id_node, new_type);
        node_type_counter->SetValue(id_node, 0);
      }

    }

    //density index from table
    EG_VTKDCN(vtkIntArray, node_specified_density, m_Grid, "node_specified_density");

    VertexMeshDensity nodeVMD = getVMD(id_node);
    int idx = nodeVMD.findSmallestVMD(m_VMDvector);
    node_specified_density->SetValue(id_node, idx);
  }

  // mesh quality
  if (!m_BCodeFeatureDefinition) {
    MeshQualityFaceOrientation mesh_quality;
    mesh_quality();
    EG_VTKDCN(vtkDoubleArray, node_mesh_quality, m_Grid, "node_mesh_quality");
    EG_FORALL_NODES(id_node, m_Grid) {
      if (node_mesh_quality->GetValue(id_node) < m_FaceOrientationThreshold) {
        node_type->SetValue(id_node, EG_SIMPLE_VERTEX);
      }
    }
  }

  /*
  // look for feature edge triple stars
  QList<vtkIdType> new_corners;
  EG_FORALL_NODES(id_node1, m_Grid) {
    if (node_type->GetValue(id_node1) == EG_FEATURE_EDGE_VERTEX) {
      int N = 0;
      for (int i = 0; i < m_Part.n2nGSize(id_node1); ++i) {
        vtkIdType id_node2 = m_Part.n2cGG(id_node1, i);
        if (node_type->GetValue(id_node2) == EG_FEATURE_EDGE_VERTEX) {
          ++N;
        }
      }
      if (N == 3) {
        new_corners.append(id_node1);
      }
    }
  }
  foreach (vtkIdType id_node, new_corners) {
    node_type->SetValue(id_node, EG_FEATURE_CORNER_VERTEX);
  }
  */

  updatePotentialSnapPoints();
}

bool SurfaceOperation::checkSnapPointPairForBcMatch(vtkIdType id_node1, vtkIdType id_node2)
{
  QSet<int> bcs1, bcs2;
  for (int i = 0; i < m_Part.n2bcGSize(id_node1); ++i) {
    bcs1.insert(m_Part.n2bcG(id_node1, i));
  }
  for (int i = 0; i < m_Part.n2bcGSize(id_node2); ++i) {
    bcs2.insert(m_Part.n2bcG(id_node2, i));
  }
  if (bcs2.contains(bcs1)) {
    return true;
  }
  return false;
}

void SurfaceOperation::updatePotentialSnapPoints()
{
  setAllSurfaceCells();  
  l2g_t nodes  = getPartNodes();

  m_PotentialSnapPoints.resize(m_Grid->GetNumberOfPoints());

  if (m_UniformSnapPoints) {
    m_PotentialSnapPoints.resize(m_Grid->GetNumberOfPoints());
    foreach( vtkIdType id_node, nodes ) {
      for (int i = 0; i < m_Part.n2nGSize(id_node); ++i) {
        m_PotentialSnapPoints[id_node].append(m_Part.n2nGG(id_node, i));
      }
    }
    return;
  }

  EG_VTKDCN(vtkCharArray, node_type, m_Grid, "node_type");

  foreach( vtkIdType id_node1, nodes ) {
    m_PotentialSnapPoints[id_node1].clear();
    char type1 = node_type->GetValue(id_node1);    
    if (type1 != EG_FIXED_VERTEX) { // fixed vertices do not have any snap-points
      QSet<vtkIdType> exclude_nodes;
      if (type1 == EG_FEATURE_EDGE_VERTEX || type1 == EG_BOUNDARY_EDGE_VERTEX) {
        for (int i = 0; i < m_Part.n2nGSize(id_node1); ++i) {
          vtkIdType id_node2 = m_Part.n2nGG(id_node1, i);
          char type2 = node_type->GetValue(id_node2);
          if (type2 == EG_FEATURE_CORNER_VERTEX || type2 == EG_FIXED_VERTEX) {
            for (int j = 0; j < m_Part.n2nGSize(id_node2); ++j) {
              vtkIdType id_node3 = m_Part.n2nGG(id_node2, j);
              exclude_nodes.insert(id_node3);
            }
          }
        }
      }
      for (int i = 0; i < m_Part.n2nGSize(id_node1); ++i) {
        vtkIdType id_node2 = m_Part.n2nGG(id_node1, i);
        char type2 = node_type->GetValue(id_node2);
        if (m_StrictFeatureSnap) {
          if (   (type1 == EG_SIMPLE_VERTEX)
                 || (type1 == EG_FEATURE_EDGE_VERTEX && (type2 == EG_FEATURE_EDGE_VERTEX || type2 == EG_FEATURE_CORNER_VERTEX))
                 || (type1 == EG_FEATURE_CORNER_VERTEX && type2 == EG_FEATURE_CORNER_VERTEX)
                 || (type1 == EG_BOUNDARY_EDGE_VERTEX && (type2 == EG_BOUNDARY_EDGE_VERTEX || type2 == EG_FIXED_VERTEX)))
          {
            if (!exclude_nodes.contains(id_node2)) {
              if (checkSnapPointPairForBcMatch(id_node1, id_node2)) {
                m_PotentialSnapPoints[id_node1].append(id_node2);
              }
            }
          }
        } else {
          if (   (type1 == EG_SIMPLE_VERTEX)
                 || (type1 == EG_FEATURE_EDGE_VERTEX)
                 || (type1 == EG_BOUNDARY_EDGE_VERTEX && (type2 == EG_BOUNDARY_EDGE_VERTEX || type2 == EG_FIXED_VERTEX)))
          {
            if (checkSnapPointPairForBcMatch(id_node1, id_node2)) {
              m_PotentialSnapPoints[id_node1].append(id_node2);
            }
          }
        }
      }

      // make sure feature edge vertices have at least two snap points ...
      /*
    if (type1 == EG_FEATURE_EDGE_VERTEX) {
      if (m_PotentialSnapPoints[id_node1].size() < 2) {
        m_PotentialSnapPoints[id_node1].clear();
      }
    }
    */
    }
  }
}

char SurfaceOperation::getNodeType(vtkIdType id_node, bool fix_unselected)
{
  if (m_Part.n2nGSize(id_node) == 0) {
    return EG_FIXED_VERTEX;
  }

  //char type = EG_SIMPLE_VERTEX;
  EG_VTKDCC(vtkIntArray, cell_code, m_Grid, "cell_code");

  for (int i = 0; i < m_Part.n2cGSize(id_node); ++i) {
    vtkIdType id_cell = m_Part.n2cGG(id_node, i);

    // fix all vertices that are part of a volume element, quad, or polygon
    if (m_Grid->GetCellType(id_cell) != VTK_TRIANGLE) {
      return EG_FIXED_VERTEX;
    }

    int bc = cell_code->GetValue(id_cell);

    // fix nodes which belong to faces with unselected boundary codes
    if (!m_BoundaryCodes.contains(bc) && fix_unselected) {
      return EG_FIXED_VERTEX;
    }

  }

  int num_bcs = m_Part.n2bcGSize(id_node);
  if (num_bcs >= 3 || num_bcs == 0) {
    return EG_FIXED_VERTEX;
  }
  if (num_bcs == 2) {
    if (m_BCodeFeatureDefinition) {
      QList<vtkIdType> neigh;
      for (int i = 0; i < m_Part.n2nGSize(id_node); ++i) {
        vtkIdType id_neigh = m_Part.n2nGG(id_node, i);
        if (getEdgeType(id_node, id_neigh) == EG_BOUNDARY_EDGE_VERTEX) {
          neigh << id_neigh;
        }
      }
      if (neigh.size() == 2) {
        if (M_PI - GeometryTools::angle(m_Grid, neigh[0], id_node, neigh[1]) >= m_EdgeAngle) {
          return EG_FIXED_VERTEX;
        }
      }
    }
    return EG_BOUNDARY_EDGE_VERTEX;
  }

  CadInterface* cad_interface = GuiMainWindow::pointer()->getCadInterface(m_Part.n2bcG(id_node,0), true);
  if (cad_interface && !m_BCodeFeatureDefinition) {

    vec3_t x, x0;
    m_Grid->GetPoint(id_node, x0.data());

    // compute average edge length
    double L = 0;
    for (int i = 0; i < m_Part.n2nGSize(id_node); ++i) {
      vec3_t xn;
      m_Grid->GetPoint(m_Part.n2nGG(id_node, i), xn.data());
      L += (x0 - xn).abs();
    }
    L /= m_Part.n2nGSize(id_node);
    L *= 0.1;
    bool convex = isConvexNode(id_node);

    x0 = cad_interface->projectNode(id_node, x0, m_NodeNormal[id_node]);
    //x0 = cad_interface->snapNode(id_node, x0);
    if (convex) {
      x = x0 - L*m_NodeNormal[id_node];
    } else {
      x = x0 + L*m_NodeNormal[id_node];
    }
    vec3_t n = cad_interface->getLastNormal();
    if (GeometryTools::angle(n, m_NodeNormal[id_node]) > 0.5*m_FeatureAngle && !cad_interface->failed()) {
      double d = L/tan(0.5*m_FeatureAngle);
      static const int num_steps = 36;
      double D_alpha = 2*M_PI/num_steps;
      vec3_t v;

      v = GeometryTools::orthogonalVector(m_NodeNormal[id_node]);
      int num_miss = 0;
      for (int i = 0; i < num_steps; ++i) {
        v = GeometryTools::rotate(v, m_NodeNormal[id_node], D_alpha);
        vec3_t xp = cad_interface->projectNode(id_node, x, v, true);
        if (cad_interface->failed()) {
          ++num_miss;
        } else {
          double l = (x - xp).abs();
          if (l >= d) {
            ++num_miss;
          }
        }
      }
      if (num_miss == 0) {
        return EG_FEATURE_CORNER_VERTEX;
      }

      v = GeometryTools::orthogonalVector(n);
      int num_hit = 0;
      for (int i = 0; i < num_steps; ++i) {
        v = GeometryTools::rotate(v, n, D_alpha);
        vec3_t xp = cad_interface->projectNode(id_node, x, v, true);
        if (!cad_interface->failed()) {
          double l = (x - xp).abs();
          if (l < d) {
            ++num_hit;
          }
        }
      }
      if (num_hit > 0) {
        return EG_FEATURE_EDGE_VERTEX;
      }
    }
  }

  return EG_SIMPLE_VERTEX;
}

int SurfaceOperation::getEdgeCells(vtkIdType id_node1, vtkIdType id_node2, QVector <vtkIdType> &EdgeCells)
{
  g2l_t _nodes = getPartLocalNodes();
  l2g_t cells  = getPartCells();
  l2l_t n2c    = getPartN2C();

  QSet<vtkIdType> S1;
  foreach (int i, n2c[_nodes[id_node1]]) {
    S1.insert(cells[i]);
  }

  QSet<vtkIdType> S2;
  foreach( int i, n2c[_nodes[id_node2]] ) {
    S2.insert(cells[i]);
  }

  S2.intersect(S1);
  EdgeCells = set2Vector(S2, false);
  return EdgeCells.size();
}

int SurfaceOperation::getEdgeCells( vtkIdType id_node1, vtkIdType id_node2, QSet <vtkIdType> &EdgeCells )
{
  g2l_t _nodes = getPartLocalNodes();
  l2g_t cells  = getPartCells();
  l2l_t n2c    = getPartN2C();

  QSet<vtkIdType> S1;
  foreach( int i, n2c[_nodes[id_node1]] ) {
    S1.insert( cells[i] );
  }

  QSet<vtkIdType> S2;
  foreach( int i, n2c[_nodes[id_node2]] ) {
    S2.insert( cells[i] );
  }

  EdgeCells = S2.intersect( S1 );
  return EdgeCells.size();
}

char SurfaceOperation::getEdgeType(vtkIdType a_node1, vtkIdType a_node2, bool fix_unselected)
{
  double cos_feature_angle = cos(this->m_FeatureAngle);
  bool feature_edges_disabled = m_FeatureAngle >= M_PI;

  // compute number of cells around edge [a_node,p2] and put them into neighbour_cells
  QVector <vtkIdType> neighbour_cells;
  int numNei = getEdgeCells(a_node1, a_node2, neighbour_cells) - 1;

  // set default value
  char edge = EG_SIMPLE_VERTEX;

  if (numNei == 0) {
    edge = EG_BOUNDARY_EDGE_VERTEX;
  } else if (numNei >= 2) {
    edge = EG_BOUNDARY_EDGE_VERTEX;
  } else if (numNei == 1) {

    // check angle between cell1 and cell2 against FeatureAngle
    double cosa = cosAngle(m_Grid, neighbour_cells[0], neighbour_cells[1]);
    if (cosa <= cos_feature_angle && !feature_edges_disabled) {
      edge = EG_FEATURE_EDGE_VERTEX;
    }

    // check the boundary codes
    EG_VTKDCC( vtkIntArray, cell_code, m_Grid, "cell_code" );
    int cell_code_0 = cell_code->GetValue( neighbour_cells[0] );
    int cell_code_1 = cell_code->GetValue( neighbour_cells[1] );
    if (cell_code_0 !=  cell_code_1) {
      edge = EG_BOUNDARY_EDGE_VERTEX;
    }

    if (fix_unselected) {
      if (!m_BoundaryCodes.contains(cell_code_0) || !m_BoundaryCodes.contains(cell_code_1)) {
        // does not make sense, but should make the points of the edge fixed
        edge = EG_FIXED_VERTEX;
      }
    }
  }

  return edge;
}

VertexMeshDensity SurfaceOperation::getVMD( vtkIdType id_node )
{
  g2l_t _nodes = getPartLocalNodes();
  l2g_t  cells = getPartCells();
  l2l_t  n2c   = getPartN2C();

  EG_VTKDCN( vtkCharArray, node_type, m_Grid, "node_type" );
  EG_VTKDCC( vtkIntArray, cell_code, m_Grid, "cell_code" );

  VertexMeshDensity VMD;
  VMD.type = node_type->GetValue( id_node );
  VMD.density = 0;
  VMD.CurrentNode = id_node;

  foreach( int i_cell, n2c[_nodes[id_node]] ) {
    vtkIdType id_cell = cells[i_cell];
    VMD.BCmap[cell_code->GetValue( id_cell )] = 2;
  }
  return( VMD );
}

//////////////////////////////////////////////
double SurfaceOperation::currentVertexAvgDist( vtkIdType id_node )
{
  l2g_t  nodes = getPartNodes();
  g2l_t _nodes = getPartLocalNodes();
  l2l_t  n2n   = getPartN2N();

  double total_dist = 0;
  double avg_dist = 0;
  int N = n2n[_nodes[id_node]].size();
  vec3_t C;
  m_Grid->GetPoint( id_node, C.data() );
  foreach( int i_node_neighbour, n2n[_nodes[id_node]] ) {
    vtkIdType id_node_neighbour = nodes[i_node_neighbour];
    vec3_t M;
    m_Grid->GetPoint( id_node_neighbour, M.data() );
    total_dist += ( M - C ).abs();
  }
  avg_dist = total_dist / ( double )N;
  return( avg_dist );
}

double SurfaceOperation::CurrentMeshDensity( vtkIdType id_node )
{
  return 1.0 / currentVertexAvgDist( id_node );
}

///\todo change meshdensity fields to edgelength fields since this is what is mostly used?

/// desired edge length for id_node
double SurfaceOperation::desiredEdgeLength( vtkIdType id_node )
{
  EG_VTKDCN( vtkDoubleArray, characteristic_length_desired, m_Grid, "node_meshdensity_desired" );
  return( 1.0 / characteristic_length_desired->GetValue( id_node ) );
}

/// mean desired edge length for id_cell
double SurfaceOperation::meanDesiredEdgeLength( vtkIdType id_cell )
{
  vtkIdType num_pts, *pts;
  m_Grid->GetCellPoints( id_cell, num_pts, pts );
  int total = 0;
  for ( int i = 0; i < num_pts; i++ ) {
    total += desiredEdgeLength( pts[i] );
  }
  return total / ( double )num_pts;
}

QVector <vtkIdType> SurfaceOperation::getPotentialSnapPoints( vtkIdType id_node )
{
  if ((id_node < 0) || (id_node >= m_PotentialSnapPoints.size())) {
    // UpdatePotentialSnapPoints should probably be called before using this function.
    EG_BUG;
  }
  return m_PotentialSnapPoints[id_node];
}

bool SurfaceOperation::isCell(vtkIdType id_node1, vtkIdType id_node2, vtkIdType id_node3)
{
  QVector <vtkIdType> EdgeCells_12;
  QVector <vtkIdType> EdgeCells_13;
  QVector <vtkIdType> inter;
  
  getEdgeCells( id_node1, id_node2, EdgeCells_12 );
  getEdgeCells( id_node1, id_node3, EdgeCells_13 );
  qcontIntersection( EdgeCells_12, EdgeCells_13, inter );
  if(inter.size()>1) {
    qWarning()<<"(id_node1, id_node2, id_node3)="<<"("<<id_node1<<", "<<id_node2<<", "<<id_node3<<")";
    qWarning()<<"EdgeCells_12="<<EdgeCells_12;
    qWarning()<<"EdgeCells_13="<<EdgeCells_13;
    qWarning()<<"inter="<<inter;
    writeGrid(m_Grid, "abort");
    EG_BUG;// multiple cells in the same place
  }
  if(DebugLevel>100) {
    qDebug()<<"(id_node1, id_node2, id_node3)="<<"("<<id_node1<<", "<<id_node2<<", "<<id_node3<<")";
    qDebug()<<"EdgeCells_12="<<EdgeCells_12;
    qDebug()<<"EdgeCells_13="<<EdgeCells_13;
    qDebug()<<"inter="<<inter;
  }
  return(inter.size()>0);
}

void SurfaceOperation::computeNormals()
{
  EG_VTKDCC(vtkIntArray, cell_code, m_Grid, "cell_code");
  m_NodeNormal.fill(vec3_t(0,0,0), m_Grid->GetNumberOfPoints());
  for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
    if (m_Part.localNode(id_node) != -1) {
      QSet<int> bcs;
      for (int i = 0; i < m_Part.n2cGSize(id_node); ++i) {
        vtkIdType id_cell = m_Part.n2cGG(id_node, i);
        if (isSurface(id_cell, m_Grid)) {
          int bc = cell_code->GetValue(id_cell);
          if (m_BoundaryCodes.contains(bc)) {
            bcs.insert(bc);
          }
        }
      }
      int num_bcs = bcs.size();
      QVector<vec3_t> normal(num_bcs, vec3_t(0,0,0));
      QMap<int,int> bcmap;
      int i_bc = 0;
      foreach (int bc, bcs) {
        bcmap[bc] = i_bc;
        ++i_bc;
      }
      for (int i = 0; i < m_Part.n2cGSize(id_node); ++i) {
        vtkIdType id_cell = m_Part.n2cGG(id_node, i);
        if (isSurface(id_cell, m_Grid)) {
          int bc = cell_code->GetValue(id_cell);
          if (m_BoundaryCodes.contains(bc)) {
            vtkIdType N_pts, *pts;
            m_Grid->GetCellPoints(id_cell, N_pts, pts);
            vec3_t a, b, c;
            for (int j = 0; j < N_pts; ++j) {
              if (pts[j] == id_node) {
                m_Grid->GetPoint(pts[j], a.data());
                if (j > 0) {
                  m_Grid->GetPoint(pts[j-1], b.data());
                } else {
                  m_Grid->GetPoint(pts[N_pts-1], b.data());
                }
                if (j < N_pts - 1) {
                  m_Grid->GetPoint(pts[j+1], c.data());
                } else {
                  m_Grid->GetPoint(pts[0], c.data());
                }
              }
            }
            vec3_t u = b - a;
            vec3_t v = c - a;
            double alpha = GeometryTools::angle(u, v);
            vec3_t n = u.cross(v);
            n.normalise();
            if (checkVector(n)) {
              normal[bcmap[bc]] -= alpha*n;
            }
          }
        }
      }
      for (int i = 0; i < num_bcs; ++i) {
        normal[i].normalise();
      }
      if (num_bcs > 0) {
        if (num_bcs > 1) {
          if (num_bcs == 3) {
            for (int i = 0; i < num_bcs; ++i) {
              for (int j = i + 1; j < num_bcs; ++j) {
                vec3_t n = normal[i] + normal[j];
                n.normalise();
                m_NodeNormal[id_node] += n;
              }
            }
          } else {
            for (int i = 0; i < num_bcs; ++i) {
              m_NodeNormal[id_node] += normal[i];
            }
          }
        } else {
          m_NodeNormal[id_node] = normal[0];
        }
        m_NodeNormal[id_node].normalise();
      }
    }
  }
}

double SurfaceOperation::normalIrregularity(vtkIdType id_node)
{
  double nirr = 0;
  QVector<vec3_t> nc(m_Part.n2cGSize(id_node));
  for (int i = 0; i < nc.size(); ++i) {
    nc[i] = GeometryTools::cellNormal(m_Grid, m_Part.n2cGG(id_node, i));
    nc[i].normalise();
  }
  for (int i = 0; i < nc.size(); ++i) {
    for (int j = i + 1; j < nc.size(); ++j) {
      nirr += 1.0 - nc[i]*nc[j];
    }
  }
  return nirr;
}

bool SurfaceOperation::isConvexNode(vtkIdType id_node)
{
  int N = m_Part.n2nGSize(id_node);
  if (N == 0) {
    return false;
  }
  vec3_t x1, x2(0,0,0);
  m_Grid->GetPoint(id_node, x1.data());
  for (int i = 0; i < N; ++i) {
    vec3_t x;
    m_Grid->GetPoint(m_Part.n2nGG(id_node, i), x.data());
    x2 += x;
  }
  x2 *= 1.0/N;
  if ((x1 - x2)*m_NodeNormal[id_node] > 0) {
    return true;
  }
  return false;
}

double SurfaceOperation::getSurfaceDeviation(vtkIdType id_node)
{
  vec3_t x;
  m_Grid->GetPoint(id_node, x.data());
  double d = 0;
  for (int i = 0; i < m_Part.n2nGSize(id_node); ++i) {
    vec3_t xe;
    m_Grid->GetPoint(m_Part.n2nGG(id_node, i), xe.data());
    d += (x - xe)*m_NodeNormal[id_node];
  }
  if (m_Part.n2nGSize(id_node) > 0) {
    d *= 1.0/m_Part.n2nGSize(id_node);
  }
  return d;
}

bool SurfaceOperation::isFeatureNode(vtkIdType id_node)
{
  EG_VTKDCN(vtkCharArray, node_type, m_Grid, "node_type");
  char type = node_type->GetValue(id_node);
  if (type == EG_FEATURE_CORNER_VERTEX || type == EG_FEATURE_EDGE_VERTEX) {
    return true;
  }
  return false;
}
