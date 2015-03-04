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
  getSet("surface meshing", "feature angle", 30, m_FeatureAngle);
  getSet("surface meshing", "threshold for face orientation quality", 0.1, m_FaceOrientationThreshold);
  m_FeatureAngle = GeometryTools::deg2rad(m_FeatureAngle);
  m_EdgeAngle = GeometryTools::deg2rad(m_EdgeAngle);
  setEdgeAngle(m_EdgeAngle);
  m_StretchingFactor = 0;
  m_UniformSnapPoints = false;

  //EG_STOPDATE("2015-03-01");
  m_StrictFeatureSnap = true;
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
    {
      m_VMDvector.fill(VertexMeshDensity(), row_count);
      for (int i = 0; i < row_count; ++i) {
        int row, column;
        QString formula;
        for (int j = 0; j < tmp_bcs.size(); ++j) {
          int bc = tmp_bcs[j];
          if (j < column_count - 3) {
            in >> row >> column >> formula;
          } else {
            formula = "1";
          }
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
    }
  }
}

/*
void SurfaceOperation::updateNodeInfo()
{
  setAllCells();
  readVMD();
  l2g_t nodes = getPartNodes();
  computeNormals();
  EG_VTKDCN(vtkCharArray, node_type, m_Grid, "node_type");
  EG_VTKDCN(vtkIntArray, node_type_counter, m_Grid, "node_type_counter");
  EG_VTKDCC(vtkIntArray, cell_code, m_Grid, "cell_code");

  bool num_non_simple = 0;
  foreach (vtkIdType id_node, nodes) {
    if (node_type->GetValue(id_node) != EG_SIMPLE_VERTEX) {
      ++num_non_simple;
    }
  }

  foreach (vtkIdType id_node, nodes) {

    char old_type = node_type->GetValue(id_node);

    if (!m_AngleFeatureDefinition) {
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
    } else {
      if (old_type == EG_SIMPLE_VERTEX) {
        char new_type = getNodeType(id_node);
        if (new_type == EG_FEATURE_CORNER_VERTEX || new_type == EG_FEATURE_EDGE_VERTEX ) {
          if (num_non_simple == 0) {
            node_type->SetValue(id_node, new_type);
          }
        } else {
          if (num_non_simple == 0) {
            node_type->SetValue(id_node, new_type);
          }
        }
      }
    }

    //density index from table
    EG_VTKDCN(vtkIntArray, node_specified_density, m_Grid, "node_specified_density");

    VertexMeshDensity nodeVMD = getVMD(id_node);
    int idx = nodeVMD.findSmallestVMD(m_VMDvector);
    node_specified_density->SetValue(id_node, idx);
  }

  // mesh quality
//  if (!m_BCodeFeatureDefinition) {
//    MeshQualityFaceOrientation mesh_quality;
//    mesh_quality();
//    EG_VTKDCN(vtkDoubleArray, node_mesh_quality, m_Grid, "node_mesh_quality");
//    EG_FORALL_NODES(id_node, m_Grid) {
//      if (node_mesh_quality->GetValue(id_node) < m_FaceOrientationThreshold) {
//        node_type->SetValue(id_node, EG_SIMPLE_VERTEX);
//      }
//    }
//  }

  updatePotentialSnapPoints();
}
*/

void SurfaceOperation::updateNodeInfo()
{
  setAllCells();
  readVMD();
  l2g_t nodes = getPartNodes();
  computeNormals();
  EG_VTKDCN(vtkCharArray, node_type, m_Grid, "node_type");
  EG_VTKDCN(vtkIntArray, node_type_counter, m_Grid, "node_type_counter");
  EG_VTKDCC(vtkIntArray, cell_code, m_Grid, "cell_code");

  // check if there are non simple vertices
  bool all_simple = true;
  foreach (vtkIdType id_node, nodes) {
    if (node_type->GetValue(id_node) != EG_SIMPLE_VERTEX) {
      all_simple = false;
      break;
    }
  }

  foreach (vtkIdType id_node, nodes) {

    if (all_simple) {
      node_type->SetValue(id_node, getNodeType(id_node, true));
    }

    //density index from table
    EG_VTKDCN(vtkIntArray, node_specified_density, m_Grid, "node_specified_density");

    VertexMeshDensity nodeVMD = getVMD(id_node);
    int idx = nodeVMD.findSmallestVMD(m_VMDvector);
    node_specified_density->SetValue(id_node, idx);
  }

  /*
  foreach (vtkIdType id_node, nodes) {

    char old_type = node_type->GetValue(id_node);
    char new_type = getNodeType(id_node, true);

    if (new_type == EG_FIXED_VERTEX) {
      if (old_type == EG_SIMPLE_VERTEX) {
        node_type->SetValue(id_node, new_type);
      }
    }

    if (new_type > old_type) {
      node_type->SetValue(id_node, new_type);
    } else {
    }

    //density index from table
    EG_VTKDCN(vtkIntArray, node_specified_density, m_Grid, "node_specified_density");

    VertexMeshDensity nodeVMD = getVMD(id_node);
    int idx = nodeVMD.findSmallestVMD(m_VMDvector);
    node_specified_density->SetValue(id_node, idx);
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
        /*
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
        */
        if (isFeatureNode(id_node1)) {
          char type2 = node_type->GetValue(id_node2);
          if (type1 == EG_FEATURE_EDGE_VERTEX && type2 >= type1) {
            QVector<vtkIdType> edge_faces;
            m_Part.getEdgeFaces(id_node1, id_node2, edge_faces);
            if (edge_faces.size() == 2) {
              vec3_t n1 = cellNormal(m_Grid, edge_faces[0]);
              vec3_t n2 = cellNormal(m_Grid, edge_faces[1]);
              if (GeometryTools::angle(n1, n2) >= m_FeatureAngle) {
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

double SurfaceOperation::edgeAngle(vtkIdType id_node1, vtkIdType id_node2)
{
  QVector<vtkIdType> faces;
  m_Part.getEdgeFaces(id_node1, id_node2, faces);
  if (faces.size() == 2) {
    vec3_t n1 = cellNormal(m_Grid, faces[0]);
    vec3_t n2 = cellNormal(m_Grid, faces[1]);
    return GeometryTools::angle(n1, n2);
  }
  return 0;
}

/*
char SurfaceOperation::getNodeType(vtkIdType id_node, bool fix_unselected)
{
  char type = EG_SIMPLE_VERTEX;

  if (m_Part.n2nGSize(id_node) == 0) {
    return EG_FIXED_VERTEX;
  }

  EG_VTKDCC(vtkIntArray, cell_code, m_Grid, "cell_code");

  for (int i = 0; i < m_Part.n2cGSize(id_node); ++i) {
    vtkIdType id_cell = m_Part.n2cGG(id_node, i);

    // fix all vertices that are part of a volume element, quad, or polygon
    if (m_Grid->GetCellType(id_cell) != VTK_TRIANGLE) {
      type = max(type, char(EG_FIXED_VERTEX));
    }

    int bc = cell_code->GetValue(id_cell);

    // fix nodes which belong to faces with unselected boundary codes
    if (!m_BoundaryCodes.contains(bc) && fix_unselected) {
      type = max(type, char(EG_FIXED_VERTEX));
    }

  }

  int num_bcs = m_Part.n2bcGSize(id_node);
  if (num_bcs >= 3 || num_bcs == 0) {
    type = max(type, char(EG_FIXED_VERTEX));
  }
  if (num_bcs == 2) {
    QList<vtkIdType> neigh;
    for (int i = 0; i < m_Part.n2nGSize(id_node); ++i) {
      vtkIdType id_neigh = m_Part.n2nGG(id_node, i);
      if (getEdgeType(id_node, id_neigh) == EG_BOUNDARY_EDGE_VERTEX) {
        neigh << id_neigh;
      }
    }
    if (neigh.size() == 2) {
      if (M_PI - GeometryTools::angle(m_Grid, neigh[0], id_node, neigh[1]) >= m_EdgeAngle) {
        type = max(type, char(EG_FIXED_VERTEX));
      }
    }
    type = max(type, char(EG_BOUNDARY_EDGE_VERTEX));
  }

  // features are defined via angles
  QList<vtkIdType> feature_neigh;
  for (int i_neigh = 0; i_neigh < m_Part.n2nGSize(id_node); ++i_neigh) {
    vtkIdType id_neigh = m_Part.n2nGG(id_node, i_neigh);
    QVector<vtkIdType> id_face;
    m_Part.getEdgeFaces(id_node, id_neigh, id_face);
    if (id_face.size() == 2) {
      vec3_t n1    = cellNormal(m_Grid, id_face[0]);
      vec3_t n2    = cellNormal(m_Grid, id_face[1]);
      double alpha = GeometryTools::angle(n1, n2);
      if (alpha >= m_FeatureAngle) {
        feature_neigh.append(id_neigh);
      }
    }
  }
  if (feature_neigh.size() == 0) {
    type = max(type, char(EG_SIMPLE_VERTEX));
  } else if (feature_neigh.size() == 2) {

    vec3_t x1, x2, xc;
    m_Grid->GetPoint(id_node, xc.data());
    m_Grid->GetPoint(feature_neigh[0], x1.data());
    m_Grid->GetPoint(feature_neigh[1], x2.data());
    double alpha1 = edgeAngle(id_node, feature_neigh[0]);
    double alpha2 = edgeAngle(id_node, feature_neigh[1]);
    if (alpha1 >= m_FeatureAngle && alpha2 >= m_FeatureAngle) {
      if (GeometryTools::angle(xc - x1, x2 - xc) > m_EdgeAngle) {
        type = max(type, char(EG_FIXED_VERTEX));
      }
    }
    type = max(type, char(EG_FEATURE_EDGE_VERTEX));
  } else {
    type = max(type, char(EG_FIXED_VERTEX));
    //return EG_FIXED_VERTEX; //EG_FEATURE_CORNER_VERTEX;
  }

  return type;
}
*/

char SurfaceOperation::getNodeType(vtkIdType id_node, bool fix_unselected)
{
  l2g_t  nodes = getPartNodes();
  g2l_t _nodes = getPartLocalNodes();
  l2l_t  n2n   = getPartN2N();

  // fix all vertices that are part of a volume element or a quad
  for (int i = 0; i < m_Part.n2cGSize(id_node); ++i) {
    if (m_Grid->GetCellType(m_Part.n2cGG(id_node, i)) != VTK_TRIANGLE) {
      return EG_FIXED_VERTEX;
    }
  }

  //initialize default value
  char type = EG_SIMPLE_VERTEX;

  //loop through edges around id_node

  QVector <vtkIdType> edges;

  double CosEdgeAngle = cos(this->m_EdgeAngle);

  foreach (int i_node2, n2n[_nodes[id_node]]) {
    vtkIdType id_node2 = nodes[i_node2];
    //-----------------------
    //determine edge type
    char edge = getEdgeType(id_node2, id_node, fix_unselected);

    //-----------------------
    //determine node type pre-processing (count nb of complex edges if the node is complex, otherwise, just count the nb of edges)
    if (edge && type == EG_SIMPLE_VERTEX) {
      edges.clear();
      edges.push_back( id_node2 );
      type = edge;
    } else if (( edge && type == EG_BOUNDARY_EDGE_VERTEX) ||
               ( edge && type == EG_FEATURE_EDGE_VERTEX) ||
               (!edge && type == EG_SIMPLE_VERTEX)) {
      edges.push_back(id_node2);
      if (type && edge == EG_BOUNDARY_EDGE_VERTEX) {
        type = EG_BOUNDARY_EDGE_VERTEX;//VTK_BOUNDARY_EDGE_VERTEX has priority over VTK_FEATURE_EDGE_VERTEX
      }
    }
  }
  //-----------------------
  //determine node type post-processing
  if (type == EG_FEATURE_EDGE_VERTEX || type == EG_BOUNDARY_EDGE_VERTEX) { //see how many edges; if two, what the angle is

    if (edges.size() == 2) { //check angle between edges
      double x1[3], x2[3], x3[3], l1[3], l2[3];
      m_Grid->GetPoint( edges[0], x1 );
      m_Grid->GetPoint( id_node, x2 );
      m_Grid->GetPoint( edges[1], x3 );
      for ( int k = 0; k < 3; k++ ) {
        l1[k] = x2[k] - x1[k];
        l2[k] = x3[k] - x2[k];
      }
      if ( vtkMath::Normalize( l1 ) >= 0.0 &&
           vtkMath::Normalize( l2 ) >= 0.0 &&
           vtkMath::Dot( l1, l2 ) < CosEdgeAngle ) {
             type = EG_FIXED_VERTEX;
      }
    } else {
      type = EG_FIXED_VERTEX;
    }

    //if along edge
  } //if edge vertex

  return type;
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

char SurfaceOperation::getEdgeType(vtkIdType id_node1, vtkIdType id_node2, bool fix_unselected)
{
  // compute number of cells around edge [a_node,p2] and put them into neighbour_cells
  QVector <vtkIdType> neighbour_cells;
  int num_neigh = getEdgeCells(id_node1, id_node2, neighbour_cells) - 1;

  // set default value
  char edge = EG_SIMPLE_VERTEX;

  if (num_neigh == 0) {
    edge = EG_BOUNDARY_EDGE_VERTEX;
  } else if (num_neigh >= 2) {
    edge = EG_BOUNDARY_EDGE_VERTEX;
  } else if (num_neigh == 1) {

    /*
    if (m_Part.isFeatureEdge(id_node1, id_node2, m_FeatureAngle)) {
      edge = EG_FEATURE_EDGE_VERTEX;
    }
    */
    vec3_t n1 = cellNormal(m_Grid, neighbour_cells[0]);
    vec3_t n2 = cellNormal(m_Grid, neighbour_cells[1]);
    if (GeometryTools::angle(n1, n2) > m_FeatureAngle) {
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
