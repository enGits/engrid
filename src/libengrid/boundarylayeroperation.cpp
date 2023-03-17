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
//
#include "boundarylayeroperation.h"
#include "optimisenormalvector.h"
#include "guimainwindow.h"
#include "pointfinder.h"
#include "cgaltricadinterface.h"
#include "geometrysmoother.h"
#include "vtkEgPolyDataToUnstructuredGridFilter.h"

#include <vtkDataSetSurfaceFilter.h>
#include <vtkIdList.h>
#include <vtkLoopSubdivisionFilter.h>
#include <vtkLinearSubdivisionFilter.h>
#include <vtkButterflySubdivisionFilter.h>
#include <vtkSmoothPolyDataFilter.h>
#include <vtkType.h>
#include <vtkWindowedSincPolyDataFilter.h>
#include <vtkSmartPointer.h>

BoundaryLayerOperation::BoundaryLayerOperation()
{
  m_ShellGrid = vtkUnstructuredGrid::New();
}

BoundaryLayerOperation::~BoundaryLayerOperation()
{
  m_ShellGrid->Delete();
}

void BoundaryLayerOperation::readSettings()
{
  EG_VTKDCC(vtkIntArray, cell_code, m_Grid, "cell_code");
  QString buffer = GuiMainWindow::pointer()->getXmlSection("engrid/blayer/settings");
  if(!buffer.isEmpty()) {
    QTextStream in(&buffer, QIODevice::ReadOnly);
    int num_bcs;
    in >> num_bcs;
    QVector<bool> use_bc(num_bcs);
    int N = 0;
    for (int i = 0; i < num_bcs; ++i) {
      int state;
      in >> state;
      use_bc[i] = bool(state);
      if (use_bc[i]) {
        ++N;
      }
    }
    m_BoundaryLayerCodes.resize(N);
    QVector<int> bcs;
    mainWindow()->getAllBoundaryCodes(bcs);
    int j = 0;
    for (int i = 0; i < bcs.size(); ++i) {
      if (use_bc[i]) {
        m_BoundaryLayerCodes[j++] = bcs[i];
      }
    }
    m_LayerAdjacentBoundaryCodes.clear();
    for (vtkIdType id_cell = 0; id_cell < m_Grid->GetNumberOfCells(); ++id_cell) {
      if (isSurface(id_cell, m_Grid)) {
        if (m_BoundaryLayerCodes.contains(cell_code->GetValue(id_cell))) {
          for (int i = 0; i < m_Part.c2cGSize(id_cell); ++i) {
            vtkIdType id_neigh = m_Part.c2cGG(id_cell, i);
            if (!m_BoundaryLayerCodes.contains(cell_code->GetValue(id_neigh))) {
              m_LayerAdjacentBoundaryCodes.insert(cell_code->GetValue(id_neigh));
            }
          }
        }
      }
    }
    in >> m_FeatureAngle;
    m_FeatureAngle = deg2rad(m_FeatureAngle);
    in >> m_StretchingRatio;
    in >> m_FarfieldRatio;
    in >> m_NumBoundaryLayerVectorRelaxations;
    in >> m_NumBoundaryLayerHeightRelaxations;
    in >> m_ShellPassBand;
    m_ShellPassBand = 2*pow(10.0, -3*m_ShellPassBand);
    in >> m_FaceSizeLowerLimit;
    in >> m_FaceSizeUpperLimit;
    in >> m_FaceAngleLimit;
    m_FaceAngleLimit = deg2rad(m_FaceAngleLimit);
    in >> m_MaxHeightInGaps;
    in >> m_RadarAngle;
    int use_grouping;
    in >> use_grouping;
    m_UseGrouping = use_grouping;
    in >> m_GroupingAngle;
    m_GroupingAngle = deg2rad(m_GroupingAngle);
  }
  m_ELSManagerBLayer.clear();
  m_ELSManagerBLayer.readBoundaryLayerRules(m_Grid);
  m_ELSManagerSurface.clear();
  m_ELSManagerSurface.read();
  m_ELSManagerSurface.readRules(m_Grid);
}

void BoundaryLayerOperation::computeBoundaryLayerVectors()
{
  EG_VTKDCC(vtkIntArray, cell_code, m_Grid, "cell_code");
  m_BoundaryLayerVectors.fill(vec3_t(0,0,0), m_Grid->GetNumberOfPoints());
  QVector<int> num_bcs(m_Grid->GetNumberOfPoints());
  QVector<OptimiseNormalVector> n_opt(m_Grid->GetNumberOfPoints(), OptimiseNormalVector(m_UseGrouping, m_GroupingAngle));
  for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
    QSet<int> bcs;
    for (int i = 0; i < m_Part.n2cGSize(id_node); ++i) {
      vtkIdType id_cell = m_Part.n2cGG(id_node, i);
      if (isSurface(id_cell, m_Grid)) {
        int bc = cell_code->GetValue(id_cell);
        if (m_BoundaryLayerCodes.contains(bc)) {
          bcs.insert(bc);
        }
      }
    }
    num_bcs[id_node] = bcs.size();
    QVector<vec3_t> normal(num_bcs[id_node], vec3_t(0,0,0));
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
        EG_GET_CELL(id_cell, m_Grid);
        vec3_t a, b, c;
        for (int j = 0; j < num_pts; ++j) {
          if (pts[j] == id_node) {
            m_Grid->GetPoint(pts[j], a.data());
            if (j > 0) {
              m_Grid->GetPoint(pts[j-1], b.data());
            } else {
              m_Grid->GetPoint(pts[num_pts-1], b.data());
            }
            if (j < num_pts - 1) {
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
        if (m_BoundaryLayerCodes.contains(bc)) {
          normal[bcmap[bc]] += alpha*n;
          n_opt[id_node].addFace(n);
        } else {
          n_opt[id_node].addConstraint(n);
        }
      }
    }
    for (int i = 0; i < num_bcs[id_node]; ++i) {
      normal[i].normalise();
    }
    if (num_bcs[id_node] > 0) {
      if (num_bcs[id_node] > 1) {
        if (num_bcs[id_node] == 3) {
          for (int i = 0; i < num_bcs[id_node]; ++i) {
            for (int j = i + 1; j < num_bcs[id_node]; ++j) {
              vec3_t n = normal[i] + normal[j];
              n.normalise();
              m_BoundaryLayerVectors[id_node] += n;
            }
          }
        } else {
          for (int i = 0; i < num_bcs[id_node]; ++i) {
            m_BoundaryLayerVectors[id_node] += normal[i];
          }
        }
      } else {
        m_BoundaryLayerVectors[id_node] = normal[0];
      }
      m_BoundaryLayerVectors[id_node].normalise();
      m_BoundaryLayerVectors[id_node] = n_opt[id_node](m_BoundaryLayerVectors[id_node]);
      m_BoundaryLayerVectors[id_node].normalise();
    }
    if (!checkVector(m_BoundaryLayerVectors[id_node])) {
      EG_ERR_RETURN("invalid layer vector computed");
    }
  }

  m_BoundaryLayerNode.fill(false, m_Grid->GetNumberOfPoints());
  for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
    for (int i = 0; i < m_Part.n2bcGSize(id_node); ++i) {
      if (m_BoundaryLayerCodes.contains(m_Part.n2bcG(id_node, i))) {
        m_BoundaryLayerNode[id_node] = true;
      }
    }
  }

  computeNodeTypes();
  computeLimitPlaneNormals();
  writeBoundaryLayerVectors("blayer", 1);
  smoothBoundaryLayerVectors(m_NumBoundaryLayerVectorRelaxations);
  writeBoundaryLayerVectors("blayer", 2);

  for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
    if (m_BoundaryLayerVectors[id_node].abs() < 0.1) {
      vec3_t n(0,0,0);
      int N = 0;
      for (int i = 0; i < m_Part.n2cGSize(id_node); ++i) {
        vtkIdType id_cell = m_Part.n2cGG(id_node, i);
        if (isSurface(id_cell, m_Grid)) {
          n += GeometryTools::cellNormal(m_Grid, id_cell);
          ++N;
        }
      }
      if (N) {
        n.normalise();
        m_BoundaryLayerVectors[id_node] = n;
      }
      if (num_bcs[id_node] > 1) {
        m_BoundaryLayerVectors[id_node] = n_opt[id_node](m_BoundaryLayerVectors[id_node]);
      }
      if (!checkVector(m_BoundaryLayerVectors[id_node])) {
        EG_ERR_RETURN("invalid layer vector computed");
      }
    }
  }

  computeHeights();
  for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
    if (m_BoundaryLayerNode[id_node]) {
      m_BoundaryLayerVectors[id_node] *= m_Height[id_node];
      if (!checkVector(m_BoundaryLayerVectors[id_node])) {
        EG_ERR_RETURN("invalid layer vector computed");
      }
    }
  }
  writeBoundaryLayerVectors("blayer", 3);

}

void BoundaryLayerOperation::computeLimitPlaneNormals()
{
  m_LimitPlaneNormals.fill(vec3_t(0,0,0), m_BoundaryLayerVectors.size());
  for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
    if (m_NodeTypes[id_node] == EdgeNode) {
      if (m_SnapPoints[id_node].size() == 2) {
        QList<vec3_t> x;
        foreach (vtkIdType id_snap, m_SnapPoints[id_node]) {
          vec3_t x_snap;
          m_Grid->GetPoint(id_snap, x_snap.data());
          x << x_snap;
        }
        vec3_t v = x[1] - x[0];
        m_LimitPlaneNormals[id_node] = m_BoundaryLayerVectors[id_node].cross(v);
        m_LimitPlaneNormals[id_node].normalise();
      }
    }
  }
}

void BoundaryLayerOperation::addToSnapPoints(vtkIdType id_node, vtkIdType id_snap)
{
  if (m_NodeTypes[id_node] == NormalNode) {
    m_SnapPoints[id_node].insert(id_snap);
  } else if (m_NodeTypes[id_node] == EdgeNode && m_NodeTypes[id_snap] != NormalNode) {
    m_SnapPoints[id_node].insert(id_snap);
  }
}

void BoundaryLayerOperation::computeNodeTypes()
{
  m_NodeTypes.fill(NormalNode, m_Grid->GetNumberOfPoints());
  m_SnapPoints.resize(m_Grid->GetNumberOfPoints());
  QVector<int> num_feature_edges(m_Grid->GetNumberOfPoints(), 0);

  for (vtkIdType id_cell1 = 0; id_cell1 < m_Grid->GetNumberOfCells(); ++id_cell1) {
    if (isSurface(id_cell1, m_Grid)) {
      vec3_t n1 = cellNormal(m_Grid, id_cell1);
      for (int i = 0; i < m_Part.c2cGSize(id_cell1); ++i) {
        vtkIdType id_cell2 = m_Part.c2cGG(id_cell1, i);
        if (id_cell2 > id_cell1) {
          if (!isSurface(id_cell2, m_Grid)) {
            EG_BUG;
          }
          vec3_t n2 = cellNormal(m_Grid, id_cell2);
          if (angle(n1, n2) >= m_FeatureAngle) {
            QVector<vtkIdType> nodes;
            m_Part.getCommonNodes(id_cell1, id_cell2, nodes);
            foreach (vtkIdType id_node, nodes) {
              ++num_feature_edges[id_node];
            }
          }
        }
      }
    }
  }
  for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
    if      (num_feature_edges[id_node] > 2) m_NodeTypes[id_node] = CornerNode;
    else if (num_feature_edges[id_node] > 1) m_NodeTypes[id_node] = EdgeNode;
    else                                     m_NodeTypes[id_node] = NormalNode;
  }
  for (vtkIdType id_cell1 = 0; id_cell1 < m_Grid->GetNumberOfCells(); ++id_cell1) {
    if (isSurface(id_cell1, m_Grid)) {
      vec3_t n1 = cellNormal(m_Grid, id_cell1);
      for (int i = 0; i < m_Part.c2cGSize(id_cell1); ++i) {
        vtkIdType id_cell2 = m_Part.c2cGG(id_cell1, i);
        if (id_cell2 > id_cell1) {
          if (!isSurface(id_cell2, m_Grid)) {
            EG_BUG;
          }
          vec3_t n2 = cellNormal(m_Grid, id_cell2);
          QVector<vtkIdType> nodes;
          m_Part.getCommonNodes(id_cell1, id_cell2, nodes);
          if (nodes.size() != 2) {
            EG_BUG;
          }
          if (angle(n1, n2) >= m_FeatureAngle) {
            addToSnapPoints(nodes[0], nodes[1]);
            addToSnapPoints(nodes[1], nodes[0]);
          } else {
            if (m_NodeTypes[nodes[0]] == NormalNode) {
              m_SnapPoints[nodes[0]].insert(nodes[1]);
            }
            if (m_NodeTypes[nodes[1]] == NormalNode) {
              m_SnapPoints[nodes[1]].insert(nodes[0]);
            }
          }
        }
      }
    }
  }
}

void BoundaryLayerOperation::correctBoundaryLayerVectors()
{
  EG_VTKDCC(vtkIntArray, bc, m_Grid, "cell_code");
  for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
    if (m_BoundaryLayerVectors[id_node].abs() > 0.1) {
      for (int iter = 0; iter < 20; ++iter) {
        for (int i = 0; i < m_Part.n2cGSize(id_node); ++i) {
          vtkIdType id_cell = m_Part.n2cGG(id_node,i);
          if (isSurface(id_cell, m_Grid)) {
            if (!m_BoundaryLayerCodes.contains(bc->GetValue(id_cell))) {
              vec3_t v = m_BoundaryLayerVectors[id_node];
              vec3_t n = GeometryTools::cellNormal(m_Grid, id_cell);
              n.normalise();
              v -= (n*m_BoundaryLayerVectors[id_node])*n;
              v.normalise();
              m_BoundaryLayerVectors[id_node] = v;
            }
          }
        }
      }
    }
  }
}

void BoundaryLayerOperation::smoothBoundaryLayerVectors(int n_iter, double w_iso, double w_dir, QVector<bool> *node_fixed)
{
  EG_VTKDCC(vtkIntArray, cell_code, m_Grid, "cell_code");

  for (int iter = 0; iter < n_iter; ++iter) {
    QVector<vec3_t> v_new = m_BoundaryLayerVectors;
    for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
      bool fixed = false;
      if (node_fixed) {
        fixed = (*node_fixed)[id_node];
      }
      if (m_BoundaryLayerNode[id_node] && !fixed) {
        v_new[id_node] = vec3_t(0,0,0);
        if (m_SnapPoints[id_node].size() > 0) {

          // check for edge between corners
          bool edge_between_corners = false;
          if (m_NodeTypes[id_node] == EdgeNode) {
            edge_between_corners = true;
            foreach (vtkIdType id_snap, m_SnapPoints[id_node]) {
              if (m_NodeTypes[id_snap] != CornerNode) {
                edge_between_corners = false;
                break;
              }
            }
          }

          bool smooth_node = !edge_between_corners;
          if (!smooth_node && m_SnapPoints[id_node].size() > 0) {

            // compute smoothed normal
            vec3_t n_smooth(0,0,0);
            foreach (vtkIdType id_snap, m_SnapPoints[id_node]) {
              n_smooth += m_BoundaryLayerVectors[id_snap];
            }
            n_smooth.normalise();

            // check if it has not been projected onto the plane of an adjacent face
            double h_min = 1.0;
            for (int i = 0; i < m_Part.n2cGSize(id_node); ++i) {
              vtkIdType id_face = m_Part.n2cGG(id_node, i);
              if (m_BoundaryLayerCodes.contains(cell_code->GetValue(id_face))) {
                vec3_t n_face = cellNormal(m_Grid, id_face);
                n_face.normalise();
                n_face *= -1;
                h_min = min(h_min, n_face*n_smooth);
              }
            }
            if (h_min > 0.2) {
              smooth_node = true;
            }

          }

          if (smooth_node) {
            v_new[id_node] = vec3_t(0,0,0);
            vec3_t x_node;
            m_Grid->GetPoint(id_node, x_node.data());
            double w_total = 0.0;
            foreach (vtkIdType id_snap, m_SnapPoints[id_node]) {
              vec3_t x_snap;
              m_Grid->GetPoint(id_snap, x_snap.data());
              QSet<vtkIdType> faces;
              getEdgeCells(id_node, id_snap, faces);
              vec3_t n_edge(0,0,0);
              foreach (vtkIdType id_face, faces) {
                n_edge += cellNormal(m_Grid, id_face);
              }
              n_edge.normalise();
              vec3_t u = m_BoundaryLayerVectors[id_snap] - (m_BoundaryLayerVectors[id_snap]*n_edge)*m_BoundaryLayerVectors[id_snap];
              vec3_t v = x_node - x_snap;
              v.normalise();
              double w = w_iso + w_dir*(u*v);
              w_total += w;
              v_new[id_node] += w*m_BoundaryLayerVectors[id_snap];
            }
            //v_new[id_node].normalise();
            v_new[id_node] *= 1.0/w_total;

            // apply limit plane if required
            /*
            if (m_LimitPlaneNormals[id_node].abs() > 0.5) {
              vec3_t dn = v_new[id_node] - m_BoundaryLayerVectors[id_node];
              dn -= (m_LimitPlaneNormals[id_node]*dn)*m_LimitPlaneNormals[id_node];
              v_new[id_node] = m_BoundaryLayerVectors[id_node] + dn;
            }
            */

          } else {
            v_new[id_node] = m_BoundaryLayerVectors[id_node];
          }
        } else {
          v_new[id_node] = m_BoundaryLayerVectors[id_node];
        }
        if (v_new[id_node].abs() < 0.1) {
          EG_BUG;
        }
      }
    }
    m_BoundaryLayerVectors = v_new;
    correctBoundaryLayerVectors();
  }
}

void BoundaryLayerOperation::writeBoundaryLayerVectors(QString file_name, int counter)
{
  if (counter >= 0) {
    QString counter_txt;
    counter_txt.setNum(counter);
    counter_txt = counter_txt.rightJustified(3, '0');
    file_name += "_" + counter_txt;
  }
  MeshPartition wall_part(m_Grid);
  wall_part.setBCs(m_BoundaryLayerCodes);
  file_name = GuiMainWindow::pointer()->getCwd() + "/" + file_name + ".vtk";
  QFile file(file_name);
  file.open(QIODevice::WriteOnly);
  QTextStream f(&file);
  f << "# vtk DataFile Version 2.0\n";
  f << "m_BoundaryLayerVectors\n";
  f << "ASCII\n";
  f << "DATASET UNSTRUCTURED_GRID\n";
  f << "POINTS " << wall_part.getNumberOfNodes() << " float\n";
  for (int i = 0; i < wall_part.getNumberOfNodes(); ++i) {
    vec3_t x;
    vtkIdType id_node = wall_part.globalNode(i);
    m_Grid->GetPoint(id_node, x.data());
    f << x[0] << " " << x[1] << " " << x[2] << "\n";
  }
  f << "CELLS " << wall_part.getNumberOfCells() << " ";
  int N = 0;
  for (int i = 0; i < wall_part.getNumberOfCells(); ++ i) {
    vtkIdType id_cell = wall_part.globalCell(i);
    EG_GET_CELL(id_cell, m_Grid);
    N += 1 + num_pts;
  }
  f << N << "\n";
  for (int i = 0; i < wall_part.getNumberOfCells(); ++ i) {
    vtkIdType id_cell = wall_part.globalCell(i);
    EG_GET_CELL(id_cell, m_Grid);
    f << num_pts;
    for (int j = 0; j < num_pts; ++j) {
      f << " " << wall_part.localNode(pts[j]);
    }
    f << "\n";
  }

  f << "CELL_TYPES " << wall_part.getNumberOfCells() << "\n";
  for (int i = 0; i < wall_part.getNumberOfCells(); ++ i) {
    vtkIdType id_cell = wall_part.globalCell(i);
    vtkIdType type_cell = m_Grid->GetCellType(id_cell);
    f << type_cell << "\n";
  }
  f << "POINT_DATA " << wall_part.getNumberOfNodes() << "\n";

  f << "VECTORS N float\n";
  for (int i = 0; i < wall_part.getNumberOfNodes(); ++i) {
    vtkIdType id_node = wall_part.globalNode(i);
    f << m_BoundaryLayerVectors[id_node][0] << " ";
    f << m_BoundaryLayerVectors[id_node][1] << " ";
    f << m_BoundaryLayerVectors[id_node][2] << "\n";
  }

  f << "SCALARS node_type int\n";
  f << "LOOKUP_TABLE default\n";
  for (int i = 0; i < wall_part.getNumberOfNodes(); ++i) {
    vtkIdType id_node = wall_part.globalNode(i);
    f << m_NodeTypes[id_node] << "\n";
  }
}

void BoundaryLayerOperation::computeDesiredHeights()
{
  EG_VTKDCN(vtkDoubleArray, cl, m_Grid, "node_meshdensity_desired");

  // first pass (intial height)
  m_Height.fill(0, m_Grid->GetNumberOfPoints());
  int k = 1;
  for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
    if (m_BoundaryLayerNode[id_node]) {
      vec3_t x;
      m_Grid->GetPoint(id_node, x.data());
      double h0 = m_ELSManagerBLayer.minEdgeLength(x);
      double h1 = cl->GetValue(id_node)*m_FarfieldRatio;
      k = max(k, int(logarithm(m_StretchingRatio, h1/h0)));
    }
  }
  for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
    if (m_BoundaryLayerNode[id_node]) {
      vec3_t x;
      m_Grid->GetPoint(id_node, x.data());
      double h0 = m_ELSManagerBLayer.minEdgeLength(x);
      double h1 = cl->GetValue(id_node)*m_FarfieldRatio;
      double s  = pow(h1/h0, 1.0/k);
      double H  = h0*(1 - pow(s, k))/(1 - s);
      if (!checkReal(H)) {
        EG_ERR_RETURN("floating point error while computing heights");
      }
      if (H < 0) {
        EG_ERR_RETURN("negative height computed");
      }
      if (H > 1000*h1) {
        cout << H << ", " << h1 << endl;
        EG_ERR_RETURN("unrealistically large height computed");
      }
      if (H < 1e-3*h0) {
        EG_ERR_RETURN("unrealistically small height computed");
      }
      if (h1 < h0) {
        QString h0_txt, h1_txt, id_txt;
        h0_txt.setNum(h0);
        h1_txt.setNum(h1);
        id_txt.setNum(id_node);
        EG_ERR_RETURN("h1 < h0 (" + h1_txt + " < " + h0_txt + ", for node " + id_txt + ")");
      }
      m_Height[id_node] = H;
    }
  }

  m_NumLayers = k;

  // correct with angle between face normal and propagation direction (node normals)
  for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
    if (m_BoundaryLayerNode[id_node]) {
      int N = 0;
      double scale = 0;
      for (int j = 0; j < m_Part.n2cGSize(id_node); ++j) {
        vtkIdType id_cell = m_Part.n2cGG(id_node, j);
        if (isSurface(id_cell, m_Grid)) {
          scale += m_BoundaryLayerVectors[id_node]*cellNormal(m_Grid, id_cell).normalise();
        }
      }
      if (N > 0) {
        scale /= N;
        m_Height[id_node] /= scale;
      }
    }
  }
}

bool BoundaryLayerOperation::faceFine(vtkIdType id_face, double scale)
{
  EG_GET_CELL(id_face, m_Grid);
  if (type_cell != VTK_TRIANGLE) {
    EG_BUG;
  }
  QVector<vec3_t> x1(num_pts);
  for (vtkIdType i = 0; i < num_pts; ++i) {
    m_Grid->GetPoint(pts[i], x1[i].data());
  }
  vec3_t n1 = triNormal(x1[0], x1[1], x1[2]);
  QVector<vec3_t> x2(num_pts);
  for (vtkIdType i = 0; i < num_pts; ++i) {
    x2[i] = x1[i] + scale*m_Height[pts[i]]*m_BoundaryLayerVectors[pts[i]];
  }
  vec3_t n2 = triNormal(x2[0], x2[1], x2[2]);
  double A1 = n1.abs();
  double A2 = n2.abs();
  if (A2/A1 >= m_FaceSizeLowerLimit && A2/A1 <= m_FaceSizeUpperLimit){
    if (angle(n1, n2) < m_FaceAngleLimit) {
      return true;
    }
  }
  return false;
}

void BoundaryLayerOperation::computeHeights()
{
  EG_VTKDCC(vtkIntArray, bc, m_Grid, "cell_code");

  computeDesiredHeights();
  cout << "initial boundary layer heights computed" << endl;
  // avoid collisions
  limitHeights(1.0);

  // limit face size and angle difference
  //limitSizeAndAngleErrors();


  // The mesh smoothing methods start here
  // General variables kind of useful for push out, smoothing, etc
  // Move to somewhere else later?
  QVector<bool> on_boundary(m_Grid->GetNumberOfPoints(), false);
  QVector<bool> is_convex(m_Grid->GetNumberOfPoints(), false);
  QVector<vec3_t> grid_pnts(m_Grid->GetNumberOfPoints(), vec3_t(0,0,0));
  {
    EG_VTKDCC(vtkIntArray, cell_code, m_Grid, "cell_code");
    for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
      m_Grid->GetPoint(id_node, grid_pnts[id_node].data());
      if (m_BoundaryLayerNode[id_node]) {
        int n_faces = m_Part.n2cGSize(id_node);
        for (int i = 0; i < n_faces; ++i) {
          vtkIdType id_cell = m_Part.n2cGG(id_node, i);
          if (!m_BoundaryLayerCodes.contains(cell_code->GetValue(id_cell))) {
            on_boundary[id_node] = true;
            break;
          }
        }
        if ( (m_NodeTypes[id_node] == EdgeNode || m_NodeTypes[id_node] == CornerNode)
             &&  m_Part.isConvexNode(id_node, m_BoundaryLayerCodes) )
        {
          is_convex[id_node] = true;
        }
      }
    }
  }

  limitHeights(1.0);

  //laplacianIntersectSmoother(on_boundary);
  //angleSmoother(on_boundary, is_convex, grid_pnts);
  smoothUsingBLVectors();
  limitHeights(1.0);

  // laplacian smoothing
  {
    QVector<double> h_safe = m_Height;
    for (int iter = 0; iter < m_NumBoundaryLayerHeightRelaxations; ++iter) {
      QVector<double> h_new = m_Height;
      for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
        if (m_BoundaryLayerNode[id_node]) {
          int count = 0;
          h_new[id_node] = 0;
          for (int i = 0; i < m_Part.n2nGSize(id_node); ++i) {
            vtkIdType id_neigh = m_Part.n2nGG(id_node, i);
            if (m_BoundaryLayerNode[id_neigh]) {
              ++count;
              h_new[id_node] += min(h_safe[id_node], m_Height[id_neigh]);
            }
          }
          if (count == 0) {
            EG_BUG;
          }
          h_new[id_node] /= count;
        }
      }
      m_Height = h_new;
    }
  }

  //limitHeights(1.0);

  cout << "heights computed" << endl;
}

void BoundaryLayerOperation::createSmoothShell()
{
  // Set grid to normal*height
  QVector<vec3_t> x_org(m_Grid->GetNumberOfPoints(), vec3_t(0,0,0));
  QVector<vec3_t> x_new(m_Grid->GetNumberOfPoints(), vec3_t(0,0,0));
  for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
    vec3_t x(0,0,0);
    m_Grid->GetPoint(id_node, x_org[id_node].data());
    x_new[id_node] = x_org[id_node] + m_Height[id_node]*m_BoundaryLayerVectors[id_node];
    m_Grid->GetPoints()->SetPoint(id_node, x_new[id_node].data());
  }

  // extract wall boundaries to a separate grid
  EG_VTKSP(vtkEgBoundaryCodesFilter, extract_wall);
  EG_VTKSP(vtkUnstructuredGrid, wall_grid);
  extract_wall->SetInputData(m_Grid);
  extract_wall->SetBoundaryCodes(m_BoundaryLayerCodes);
  extract_wall->Update();

  EG_VTKSP(vtkDataSetSurfaceFilter, grid_to_pdata);
  grid_to_pdata->SetInputConnection(extract_wall->GetOutputPort());
  EG_VTKSP(vtkLinearSubdivisionFilter, subdiv);
  subdiv->SetInputConnection(grid_to_pdata->GetOutputPort());
  subdiv->SetNumberOfSubdivisions(1);

  EG_VTKSP(vtkWindowedSincPolyDataFilter, smooth);
  smooth->SetInputConnection(subdiv->GetOutputPort());
  smooth->BoundarySmoothingOn();
  smooth->FeatureEdgeSmoothingOn();
  smooth->SetFeatureAngle(180);
  smooth->SetEdgeAngle(180);
  smooth->SetNumberOfIterations(100);
  smooth->NormalizeCoordinatesOn();
  double pb = m_ShellPassBand;
  cout << "pass-band = " << pb << endl;
  smooth->SetPassBand(pb);
  smooth->Update();

  EG_VTKSP(vtkEgPolyDataToUnstructuredGridFilter, pdata_to_grid);
  pdata_to_grid->SetInputConnection(smooth->GetOutputPort());
  pdata_to_grid->Update();
  makeCopy(pdata_to_grid->GetOutput(), m_ShellGrid);

  writeGrid(m_ShellGrid, "shell");

  // reset grid to original points
  for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
    if (m_BoundaryLayerNode[id_node]) {
      m_Grid->GetPoints()->SetPoint(id_node, x_org[id_node].data());
    }
  }
}

double BoundaryLayerOperation::largestAngle(vtkIdType id_node1, vtkIdType id_node2)
{
  double alpha = 0;
  QSet<vtkIdType> faces;
  if (id_node1 == id_node2) {
    for (int i = 0; i < m_Part.n2cGSize(id_node1); ++i) {
      faces.insert(m_Part.n2cGG(id_node1, i));
    }
  } else {
    getEdgeCells(id_node1, id_node2, faces);
  }
  foreach (vtkIdType id_face, faces) {
    vec3_t n = (-1)*cellNormal(m_Grid, id_face);
    alpha = max(alpha, GeometryTools::angle(n, m_BoundaryLayerVectors[id_node1]));
  }
  return alpha;
}

void BoundaryLayerOperation::fixBoundaryLayerVectors(const QList<vtkIdType> &bad_cells, int num_smooth_iter)
{
  QVector<vtkIdType> bad_nodes;
  getNodesFromCells(bad_cells, bad_nodes, m_Grid);
  int num_changes;
  do {
    num_changes = 0;
    foreach (vtkIdType id_node, bad_nodes) {
      foreach (vtkIdType id_snap, m_SnapPoints[id_node]) {
        double node_alpha = largestAngle(id_node, id_snap);
        double snap_alpha = largestAngle(id_snap, id_node);
        if (snap_alpha > node_alpha) {
          vec3_t dv = m_BoundaryLayerVectors[id_snap] - m_BoundaryLayerVectors[id_node];
          if (dv.abs() > 1e-4) {
            m_BoundaryLayerVectors[id_node] = m_BoundaryLayerVectors[id_snap];
            ++num_changes;
          }
        }
      }
    }
  } while (num_changes);
  QVector<bool> node_fixed(m_Grid->GetNumberOfPoints(), false);
  foreach (vtkIdType id_node, bad_nodes) {
    node_fixed[id_node] = true;
  }
  correctBoundaryLayerVectors();
  smoothBoundaryLayerVectors(num_smooth_iter, 1.0, 0.0, &node_fixed);
}

void BoundaryLayerOperation::writeWallGrid(QString file_name, int counter)
{
  if (counter >= 0) {
    QString counter_txt;
    counter_txt.setNum(counter);
    counter_txt = counter_txt.rightJustified(3, '0');
    file_name += "_" + counter_txt;
  }
  MeshPartition wall_part(m_Grid);
  wall_part.setBCs(m_BoundaryLayerCodes);
  EG_VTKSP(vtkUnstructuredGrid, wall_grid);
  wall_part.extractToVtkGrid(wall_grid);
  writeGrid(wall_grid, file_name);
}

void BoundaryLayerOperation::smoothUsingBLVectors()
{
  // create shell
  createSmoothShell();

  newHeightFromShellIntersect(1.0);

  EG_VTKDCC(vtkIntArray, cell_code, m_Grid, "cell_code");

  QList<vtkIdType> bad_cells;

  double w_iso = 1.0;
  double w_dir = 0.0;
  double dw    = 0.01;

  int iter = 0;

  int last_num_bad = m_Grid->GetNumberOfCells();

  //snapAllVectorsToShell(m_ShellGrid);
  //return;

  smoothBoundaryLayerVectors(m_NumBoundaryLayerVectorRelaxations, w_iso, w_dir);

  QVector<vec3_t> old_vecs = m_BoundaryLayerVectors;
  QVector<double> old_height = m_Height;

  do {
    bad_cells.clear();
    for (vtkIdType id_cell = 0; id_cell < m_Grid->GetNumberOfCells(); ++id_cell) {
      if (m_BoundaryLayerCodes.contains(cell_code->GetValue(id_cell))) {
        if (!faceFine(id_cell, 1.0)) {
          bad_cells.append(id_cell);
        }
      }
    }
    cout << "found " << bad_cells.size() << " distorted faces" << endl;

    if (bad_cells.size() > 0) {
      if (bad_cells.size() < last_num_bad) {
        last_num_bad = bad_cells.size();
        old_vecs = m_BoundaryLayerVectors;
        old_height = m_Height;
        smoothBoundaryLayerVectors(m_NumBoundaryLayerVectorRelaxations, w_iso, w_dir);
        newHeightFromShellIntersect(1.0);
      } else {
        cout << "cannot fix completely -- terminating the loop!" << endl;
        //cout << "moving to global under relaxation now ..." << endl;
        m_BoundaryLayerVectors = old_vecs;
        m_Height = old_height;
        break;
      }
    }
    w_iso -= dw;
    w_dir += dw;
    ++iter;
  } while (bad_cells.size() > 0 && w_iso >= 0.0);

  /*
  double relax = 0.95;

  do {

    cout << "relaxation factor: " << relax << endl;
    newHeightFromShellIntersect(shell_grid, relax);

    bad_cells.clear();
    for (vtkIdType id_cell = 0; id_cell < m_Grid->GetNumberOfCells(); ++id_cell) {
      if (m_BoundaryLayerCodes.contains(cell_code->GetValue(id_cell))) {
        if (!faceFine(id_cell, 1.0)) {
          bad_cells.append(id_cell);
        }
      }
    }
    cout << "found " << bad_cells.size() << " distorted faces" << endl;
    relax -= 0.05;
  } while (bad_cells.size() > 0 && relax >= 0.25);
  */

  //swapEdgesToMatchShell(shell_grid, deg2rad(5.0));
}

bool BoundaryLayerOperation::checkVectorForNode(vec3_t v, vtkIdType id_node)
{
  for (int i = 0; i < m_Part.n2cGSize(id_node); ++i) {
    vtkIdType id_face = m_Part.n2cGG(id_node, i);
    EG_VTKDCC(vtkIntArray, cell_code, m_Grid, "cell_code");
    if (m_BoundaryLayerCodes.contains(cell_code->GetValue(id_face))) {
      vec3_t n = cellNormal(m_Grid, id_face);
      if (n*v > 0) {
        return false;
      }
    }
  }
  return true;
}

vec3_t BoundaryLayerOperation::snapToShell(CadInterface* cad, vtkIdType id_node)
{
  bool dbg = false;

  if (id_node == 0) {
    cout << "break" << endl;
    dbg = true;
  }

  vec3_t x_node;
  m_Grid->GetPoint(id_node, x_node.data());
  vec3_t x_snap = cad->snap(x_node);
  if (dbg) cout << x_snap << endl;
  if (checkVectorForNode(x_snap - x_node, id_node)) {
    return x_snap;
  }

  // initial guess
  {
    x_snap = x_node + m_Height[id_node]*m_BoundaryLayerVectors[id_node];
    QVector<QPair<vec3_t, vtkIdType> > intersections;
    cad->computeIntersections(x_node, m_BoundaryLayerVectors[id_node], intersections);
    double h_min = EG_LARGE_REAL;
    bool found = false;
    vec3_t shell_vector = m_BoundaryLayerVectors[id_node];
    for (int i = 0; i < intersections.size(); ++i) {
      vec3_t xi = intersections[i].first;
      if (dbg) cout << "xi=" << xi << endl;
      vec3_t vi = xi - x_node;
      if (vi*m_BoundaryLayerVectors[id_node] > 0) {
        double h = vi.abs();
        if (h < h_min) {
          if (dbg) cout << "h=" << h << endl;
          found = true;
          shell_vector = vi;
          h_min = h;
        }
      }
    }
    if (found) {
      x_snap = x_node + shell_vector;
    }
  }
  if (dbg) cout << x_snap << endl;


  double dist_old = 0;
  double dist_new = m_Height[id_node];
  while (fabs(dist_new - dist_old) > 1e-3*m_Height[id_node]) {
    double w = 1.0;
    vec3_t x_new = x_snap;
    do {
      w *= 0.5;
      if (w < 1e-3) {
        QString msg;
        msg.setNum(id_node);
        msg = "unable to snap node " + msg + " to shell";
        EG_ERR_RETURN(msg);
      }
      x_new = cad->snap(w*x_node + (1-w)*x_snap);
      if (dbg) cout << "w=" << w << ", x_new=" << x_new << endl;
    } while (!checkVectorForNode(x_new - x_node, id_node));
    dist_old = dist_new;
    dist_new = (x_new - x_snap).abs();
    x_snap = x_new;
  }

  return x_snap;
}

void BoundaryLayerOperation::snapAllVectorsToShell(vtkUnstructuredGrid *shell_grid)
{
  CgalTriCadInterface cad(shell_grid);
  for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
    if (m_BoundaryLayerNode[id_node]) {
      //if (m_NodeTypes[id_node] == CornerNode || m_NodeTypes[id_node] == EdgeNode) {
      vec3_t x1;
      m_Grid->GetPoint(id_node, x1.data());
      vec3_t x2 = snapToShell(&cad, id_node);
      m_BoundaryLayerVectors[id_node] = x2 - x1;
      m_Height[id_node] = m_BoundaryLayerVectors[id_node].abs();
      m_BoundaryLayerVectors[id_node].normalise();
    }
  }
  correctBoundaryLayerVectors();
}

// Compute intersection points with a shell following m_BoundaryLayerVector
// Updates m_Height
void BoundaryLayerOperation::newHeightFromShellIntersect(double relax)
{
  CgalTriCadInterface cad(m_ShellGrid);
  for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
    if (m_BoundaryLayerNode[id_node]) {
      vec3_t x;
      m_Grid->GetPoint(id_node, x.data());
      QVector<QPair<vec3_t, vtkIdType> > intersections;
      cad.computeIntersections(x, m_BoundaryLayerVectors[id_node], intersections);
      double h = 2*m_Height[id_node];
      bool found = false;
      vec3_t layer_vector = m_BoundaryLayerVectors[id_node];
      for (int i = 0; i < intersections.size(); ++i) {
        vec3_t xi = intersections[i].first;
        vec3_t vi = xi - x;
        if (vi*m_BoundaryLayerVectors[id_node] > 0) {
          double hi = vi.abs();
          if (hi < h) {
            h = hi;
            layer_vector = vi.normalise();
            found = true;
          }
        }
      }
      if (found) {
        m_Height[id_node] = relax*h;
        m_BoundaryLayerVectors[id_node] = layer_vector;
      }
    }
  }
}

void BoundaryLayerOperation::limitSizeAndAngleErrors()
{
  bool done;
  int iter = 0;
  do {
    done = true;
    QVector<double> new_scale(m_Grid->GetNumberOfPoints(), 1.0);
    QVector<int> count(m_Grid->GetNumberOfPoints(), 1);
    for (vtkIdType id_cell = 0; id_cell < m_Grid->GetNumberOfCells(); ++id_cell) {
      if (isSurface(id_cell, m_Grid)) {
        EG_GET_CELL(id_cell, m_Grid);
        bool check_face = true;
        for (vtkIdType i = 0; i < num_pts; ++i) {
          if (!m_BoundaryLayerNode[pts[i]]) {
            check_face = false;
            break;
          }
        }
        if (check_face) {
          if (!faceFine(id_cell, 1)) {
            done = false;
            double scale1 = 0.1;
            double scale2 = 1;

            bool found;
            do {
              found = false;
              int num_steps = 10;
              for (int i = 1; i <= num_steps; ++i) {
                double s = scale2 - i*(scale2 - scale1)/num_steps;
                if (faceFine(id_cell, s)) {
                  found = true;
                  scale1 = s;
                  scale2 -= (i-1)*(scale2 - scale1)/num_steps;
                  break;
                }
              }
              if (!found) {
                scale1 = scale2;
              }

            } while ((scale2 - scale1) > 1e-4 && found);

            double scale = 0.5*(scale1 + scale2);
            for (vtkIdType i = 0; i < num_pts; ++i) {
              new_scale[pts[i]] += scale;
              ++count[pts[i]];
            }
          }
        }
      }
    }
    double relax = 0.2;
    for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
      double h_new = m_Height[id_node]*new_scale[id_node]/count[id_node];
      m_Height[id_node] = relax*h_new + (1 - relax)*m_Height[id_node];
    }
    //done = true;
    ++iter;
    if (iter >= 10) {
      done = true;
    }
  } while (!done);
}


void BoundaryLayerOperation::angleSmoother(const QVector<bool>& on_boundary, const QVector<bool>& is_convex, QVector<vec3_t>& grid_pnts)
{
  int n_iter = 20;
  double weight_const = 1.;
  const double PI = 3.14159265359;

  EG_VTKDCC(vtkIntArray, cell_code, m_Grid, "cell_code");
  for (int iter = 0; iter < n_iter; ++iter) {
    // Set points to bl_normal*height
    for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
      if (m_BoundaryLayerNode[id_node]) {
        vec3_t x(0,0,0);
        m_Grid->GetPoint(id_node, x.data());
        x += m_Height[id_node]*m_BoundaryLayerVectors[id_node];
        m_Grid->GetPoints()->SetPoint(id_node, x.data());
      }
    }

    QVector<int> move_count(grid_pnts.size());
    QVector<vec3_t> grid_smoothed(grid_pnts.size(), vec3_t(0,0,0));
    for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
      if (m_BoundaryLayerNode[id_node]) {
        for (vtkIdType i = 0; i < m_Part.n2nGSize(id_node); ++i) {
          vtkIdType id_neigh = m_Part.n2nGG(id_node, i);
          if (!m_BoundaryLayerNode[id_neigh]) continue;
          if (id_neigh < id_node) continue;
          if (on_boundary[id_node] && on_boundary[id_neigh]) continue;

          QList<vtkIdType> edge_faces;
          m_Part.getEdgeFaces(id_node, id_neigh, edge_faces);
          // Test for 2 faces
          int count = 0;
          for (int j = 0; j < edge_faces.size(); ++j) {
            if (m_BoundaryLayerCodes.contains(cell_code->GetValue(edge_faces[j]))) {
              ++count;
            }
          }
          if (count < 2 ) continue;
          if (edge_faces.size() > 2) EG_BUG;

          // Prepare cell info
          vtkIdType id_cell_1 = edge_faces[0];
          vtkIdType id_cell_2 = edge_faces[1];

          vtkSmartPointer<vtkIdList> pts_c1 = vtkSmartPointer<vtkIdList>::New();
          vtkSmartPointer<vtkIdList> pts_c2 = vtkSmartPointer<vtkIdList>::New();
          m_Grid->GetCellPoints(id_cell_1, pts_c1);
          m_Grid->GetCellPoints(id_cell_2, pts_c2);
          vtkIdType npts_c1 = pts_c1->GetNumberOfIds();
          vtkIdType npts_c2 = pts_c2->GetNumberOfIds();

          if (npts_c1 != 3 || npts_c2 !=3) EG_BUG;

          vtkIdType id_n3_c1;
          vtkIdType id_n3_c2;
          for (int j = 0; j < npts_c1; j++) {
            if ( pts_c1->GetId(j) != id_node && pts_c1->GetId(j) != id_neigh) {
              id_n3_c1 = pts_c1->GetId(j);
            }
          }
          for (int j = 0; j < npts_c2; j++) {
            if ( pts_c2->GetId(j) != id_node && pts_c2->GetId(j) != id_neigh) {
              id_n3_c2 = pts_c2->GetId(j);
            }
          }

          vec3_t normal_c1 = cellNormal(m_Grid, id_cell_1);
          vec3_t normal_c2 = cellNormal(m_Grid, id_cell_2);

          double angle = GeometryTools::angle(normal_c1, normal_c2);
          if (rad2deg(angle) < 1) continue;

          double spring_angle = weight_const*angle*angle/(PI*PI);

          vec3_t x_node(0,0,0);
          vec3_t x_neigh(0,0,0);
          vec3_t x_n3_c1(0,0,0);
          vec3_t x_n3_c2(0,0,0);
          m_Grid->GetPoint(id_node,  x_node.data());
          m_Grid->GetPoint(id_neigh, x_neigh.data());
          m_Grid->GetPoint(id_n3_c1, x_n3_c1.data());
          m_Grid->GetPoint(id_n3_c2, x_n3_c2.data());

          vec3_t axis = x_node.cross(x_neigh);
          vec3_t v1 = x_n3_c1 - x_node;
          vec3_t v2 = x_n3_c2 - x_node;

          vec3_t cross_vector = axis.cross(v1);
          vec3_t v1_rot(0,0,0);
          vec3_t v2_rot(0,0,0);
          if (cross_vector*normal_c1 > 0) {
            v1_rot = GeometryTools::rotate(v1, axis,  spring_angle);
            v2_rot = GeometryTools::rotate(v2, axis, -spring_angle);
          }
          else {
            v1_rot = GeometryTools::rotate(v1, axis, -spring_angle);
            v2_rot = GeometryTools::rotate(v2, axis,  spring_angle);
          }

          grid_smoothed[id_n3_c1] += x_node + v1_rot;
          move_count[id_n3_c1] += 1;
          grid_smoothed[id_n3_c2] += x_node + v2_rot;
          move_count[id_n3_c2] += 1;
        }
      }
    }

    QVector<double> h_new = m_Height;
    QVector<vec3_t> new_BoundaryLayerVectors = m_BoundaryLayerVectors;
    for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
      m_Grid->GetPoints()->SetPoint(id_node, grid_pnts[id_node].data());
      if (m_BoundaryLayerNode[id_node]) {
        if (move_count[id_node] > 0) {
          grid_smoothed[id_node] *= 1.0/move_count[id_node];

          vec3_t new_norm = grid_smoothed[id_node] - grid_pnts[id_node];
          h_new[id_node] = new_norm.abs();
          new_norm.normalise();
          new_BoundaryLayerVectors[id_node] = new_norm;
        }
      }
    }
    double max_diff = 0;
    for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
      double diff = m_Height[id_node] - h_new[id_node];
      diff = std::sqrt(diff*diff);
      max_diff = std::max(max_diff, diff);
    }
    cout << "==========================  max diff->" << max_diff << endl;
    m_Height = h_new;
    m_BoundaryLayerVectors = new_BoundaryLayerVectors;
  }
}

void BoundaryLayerOperation::limitHeights(double safety_factor)
{
  EG_VTKDCC(vtkIntArray, bc, m_Grid, "cell_code");

  // save original node positions to x_old
  QVector<vec3_t> x_old(m_Grid->GetNumberOfPoints());
  for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
    m_Grid->GetPoint(id_node, x_old[id_node].data());
  }

  double beta = m_MaxHeightInGaps/(1.0 - m_MaxHeightInGaps);
  QVector<double> h_save = m_Height;

  int max_pass = 5;
  for (int pass = 0; pass <= max_pass; ++pass) {

    // move nodes (small steps per pass)
    {
      double w = double(pass)/max_pass;
      for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
        if (m_BoundaryLayerNode[id_node]) {
          vec3_t x = x_old[id_node] + w*m_Height[id_node]*m_BoundaryLayerVectors[id_node];
          m_Grid->GetPoints()->SetPoint(id_node, x.data());
        }
      }
    }

    CgalTriCadInterface cad(m_Grid);
    m_Height = h_save;

    for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
      if (m_BoundaryLayerNode[id_node]) {
        QList<vtkIdType> cells_of_node;
        for (int i = 0; i < m_Part.n2cGSize(id_node); ++i) {
          cells_of_node.append(m_Part.n2cGG(id_node, i));
        }
        QVector<QPair<vec3_t, vtkIdType> > intersections;

        vec3_t x_start = x_old[id_node];
        foreach (int adj_bc, m_LayerAdjacentBoundaryCodes) {
          if (m_Part.hasBC(id_node, adj_bc)) {
            int count = 0;
            vec3_t xs(0, 0, 0);
            for (int i = 0; i < m_Part.n2cGSize(id_node); ++i) {
              vtkIdType id_face = m_Part.n2cGG(id_node, i);
              if (!m_LayerAdjacentBoundaryCodes.contains(bc->GetValue(id_face))) {
                vec3_t x = cellCentre(m_Grid, id_face);
                xs += x;
                ++count;
              }
            }
            if (count > 0) {
              double w = 0.1;
              xs *= 1.0/count;
              x_start = w*xs + (1-w)*x_start;
            }
            break;
          }
        }

        cad.computeIntersections(x_start, m_BoundaryLayerVectors[id_node], intersections);
        for (int i = 0; i < intersections.size(); ++i) {
          QPair<vec3_t,vtkIdType> inters = intersections[i];
          vec3_t xi = inters.first;
          vtkIdType id_tri = inters.second;
          if (!cells_of_node.contains(id_tri)) {

            double crit_angle = deg2rad(200.0); // consider all intersections

            if (m_LayerAdjacentBoundaryCodes.contains(bc->GetValue(id_tri))) {
              crit_angle = deg2rad(85.0); // different angle for adjacent boundaries
            }

            vec3_t dx = xi - x_old[id_node];
            double alpha = angle(dx, cellNormal(m_Grid, id_tri));
            if (dx*m_BoundaryLayerVectors[id_node] > 0 && alpha < crit_angle) {
              double h_max = safety_factor*beta*dx.abs();
              if (h_max < m_Height[id_node]) {
                m_Height[id_node] = h_max;
              }
            }
          }
        }
      }
    }
  }

  // reset node positions
  for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
    m_Grid->GetPoints()->SetPoint(id_node, x_old[id_node].data());
  }

}

void BoundaryLayerOperation::laplacianIntersectSmoother(const QVector<bool>& on_boundary)
{
  int n_loops = 1;
  int n_iter = 4;

  // Set grid to normal*height
  // And set weight factor of edges and corners to 1.
  QVector<vec3_t> org_grid(m_Grid->GetNumberOfPoints(), vec3_t(0,0,0));
  QVector<vec3_t> new_grid(m_Grid->GetNumberOfPoints(), vec3_t(0,0,0));
  for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
    vec3_t x(0,0,0);
    m_Grid->GetPoint(id_node, x.data());
    org_grid[id_node] = x;
  }

  // Laplacian on points
  for (int loops = 0; loops < n_loops; ++loops) {
    for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
      if (m_BoundaryLayerNode[id_node]) {
        vec3_t x(0,0,0);
        m_Grid->GetPoint(id_node, x.data());
        x += m_Height[id_node]*m_BoundaryLayerVectors[id_node];
        m_Grid->GetPoints()->SetPoint(id_node, x.data());
        new_grid[id_node] =  x;
      }
    }
    for (int iter = 0; iter < n_iter; iter++) {
      for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
        if (m_BoundaryLayerNode[id_node]) {

          // !!!!!!!   Warning, this might not be robust enough for complex geometries !!!!!!!
          bool shared_boundary = false;
          {
            QVector<int> bc_ids;
            for (int i = 0; i < m_Part.n2bcGSize(id_node); ++i) {
              int bc = m_Part.n2bcG(id_node, i);
              if (!m_BoundaryLayerCodes.contains(bc)
                  && !bc_ids.contains(bc) )
              {
                bc_ids.append(bc);
              }
            }

            if ( bc_ids.size() > 1 ) continue;

            if (bc_ids.size() == 1) {
              shared_boundary = true;
            }
          }

          vec3_t avg_pnt(0,0,0);
          int count = 0;
          if (!shared_boundary) {
            double area = 0;
            for (int i = 0; i < m_Part.n2cGSize(id_node); ++i) {
              vtkIdType id_cell = m_Part.n2cGG(id_node, i);
              vec3_t cell_ctr = cellCentre(m_Grid, id_cell);
              double cell_area = cellVA(m_Grid, id_cell);
              avg_pnt += cell_area*cell_ctr;
              area += cell_area;
              ++count;
            }
            if (count == 0)
              continue;
            avg_pnt *= 1.0/area;
          }
          if (shared_boundary) {
            for (int i = 0; i < m_Part.n2nGSize(id_node); ++i) {
              vtkIdType id_neigh = m_Part.n2nGG(id_node, i);
              if (on_boundary[id_neigh]) {
                vec3_t x(0,0,0);
                m_Grid->GetPoint(id_neigh, x.data());
                avg_pnt += x;
                ++count;
              }
            }
            if (count == 0)
              continue;
            avg_pnt *= 1.0/count;
          }

          if (on_boundary[id_node]) {
            vec3_t node_pnt(0,0,0);
            m_Grid->GetPoint(id_node, node_pnt.data());
            double n = 1./3.;
            avg_pnt = (1-n)*node_pnt + n*(avg_pnt);
          }
          //vec3_t new_vec = avg_pnt - org_grid[id_node];
          //m_Height[id_node] = new_vec.abs();
          //new_vec.normalise();
          //m_BoundaryLayerVectors[id_node] = new_vec;
          new_grid[id_node] = avg_pnt;
        }
      }
      // Set grid to new points
      for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
        if (m_BoundaryLayerNode[id_node]) {
          m_Grid->GetPoints()->SetPoint(id_node, new_grid[id_node].data());
        }
      }
    }
    return;

    EG_VTKSP(vtkUnstructuredGrid, bl_grid);
    EG_VTKDCC(vtkIntArray, cell_code, m_Grid, "cell_code");
    MeshPartition bl_part(m_Grid);
    QList<vtkIdType> bl_faces;
    for (vtkIdType id_cell = 0; id_cell < m_Grid->GetNumberOfCells(); ++id_cell) {
      if (isSurface(id_cell, m_Grid)) {
        if (m_BoundaryLayerCodes.contains(cell_code->GetValue(id_cell))) {
          bl_faces.append(id_cell);
        }
      }
    }
    bl_part.setCells(bl_faces);
    bl_part.extractToVtkGrid(bl_grid);
    CgalTriCadInterface cad(bl_grid);
    for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
      if (m_BoundaryLayerNode[id_node]) {
        vec3_t org_pnt = org_grid[id_node];
        QVector<QPair<vec3_t, vtkIdType> > intersections;
        cad.computeIntersections(org_pnt, m_BoundaryLayerVectors[id_node], intersections);

        QVector<QPair<double, vec3_t> > intersect_vec;
        for (int i = 0; i < intersections.size(); ++i) {
          vec3_t xi = intersections[i].first;
          vec3_t new_vec = xi - org_pnt;
          if (new_vec*m_BoundaryLayerVectors[id_node] < 0) {
            continue;
          }
          double length = new_vec.abs();
          intersect_vec.append(QPair<double, vec3_t>(length, xi));
        }
        vec3_t new_pnt(0,0,0);
        double length = 1e99;
        for (int i = 0; i < intersect_vec.size(); ++i) {
          if (intersect_vec[i].first < length) {
            length = intersect_vec[i].first;
            new_pnt = intersect_vec[i].second;
          }
        }
        //vec3_t new_vec = new_pnt - org_pnt;
        //m_Height[id_node] = new_vec.abs();
        //new_vec.normalise();
        //m_BoundaryLayerVectors[id_node] = new_vec;
        m_Height[id_node] = length;
      }
    }

    // Set grid to original points before return
    for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
      if (m_BoundaryLayerNode[id_node]) {
        vec3_t x = org_grid[id_node];
        m_Grid->GetPoints()->SetPoint(id_node, x.data());
      }
    }


    QVector<double> projected_height = m_Height;
    limitSizeAndAngleErrors();
    return;

    QVector<double> height_diff(m_Height.size(), 0.0);
    for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
      if (m_BoundaryLayerNode[id_node]) {
        height_diff[id_node] = projected_height[id_node] - m_Height[id_node];
      }
    }

    int n_iter_smooth = 0;
    if (loops == 1) {
      n_iter_smooth = 1000;
    }
    else {
      n_iter_smooth = 1000;
    }
    for (int iter = 0; iter < n_iter_smooth; iter++) {
      QVector<double> new_height = height_diff;
      //QVector<double> new_height(m_Height.size(), 0.0);
      for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
        if (m_BoundaryLayerNode[id_node]) {
          int count = 1;
          double org_height = height_diff[id_node];
          double avg_height = 0;
          for(int i = 0; i < m_Part.n2nGSize(id_node); ++i) {
            vtkIdType id_neigh = m_Part.n2nGG(id_node, i);
            avg_height += height_diff[id_neigh];
            ++count;
          }
          avg_height /= count;
          if (org_height < avg_height) {
            new_height[id_node] = avg_height;
          }
        }
      }
      height_diff = new_height;
    }
    for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
      if (m_BoundaryLayerNode[id_node]) {
        m_Height[id_node] = projected_height[id_node] - height_diff[id_node];
      }
    }
  // End of global pass loop:
  }

  for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
    if (m_BoundaryLayerNode[id_node]) {
      vec3_t x = org_grid[id_node];
      m_Grid->GetPoints()->SetPoint(id_node, x.data());
    }
  }
}

bool BoundaryLayerOperation::swapRequired(stencil_t stencil, CadInterface *cad, double threshold_angle)
{

  QVector<vec3_t> x(4);
  m_Grid->GetPoint(stencil.id_node[0], x[0].data());
  m_Grid->GetPoint(stencil.p1,         x[1].data());
  m_Grid->GetPoint(stencil.id_node[1], x[2].data());
  m_Grid->GetPoint(stencil.p2,         x[3].data());

  // centres of triangles
  vec3_t xc_tri_11 = 0.333333*(x[0] + x[1] + x[3]);
  vec3_t xc_tri_12 = 0.333333*(x[1] + x[2] + x[3]);
  vec3_t xc_tri_21 = 0.333333*(x[0] + x[1] + x[2]);
  vec3_t xc_tri_22 = 0.333333*(x[2] + x[3] + x[0]);

  cad->snap(xc_tri_11);
  vec3_t n_snap_11 = cad->getLastNormal();
  cad->snap(xc_tri_12);
  vec3_t n_snap_12 = cad->getLastNormal();
  cad->snap(xc_tri_21);
  vec3_t n_snap_21 = cad->getLastNormal();
  cad->snap(xc_tri_22);
  vec3_t n_snap_22 = cad->getLastNormal();

  vec3_t n_tri_11 = GeometryTools::triNormal(x[0], x[1], x[3]).normalise();
  vec3_t n_tri_12 = GeometryTools::triNormal(x[1], x[2], x[3]).normalise();
  vec3_t n_tri_21 = GeometryTools::triNormal(x[0], x[1], x[2]).normalise();
  vec3_t n_tri_22 = GeometryTools::triNormal(x[2], x[3], x[0]).normalise();

  double alpha_11 = GeometryTools::angle(n_tri_11, n_snap_11);
  double alpha_12 = GeometryTools::angle(n_tri_12, n_snap_12);
  double alpha_21 = GeometryTools::angle(n_tri_21, n_snap_21);
  double alpha_22 = GeometryTools::angle(n_tri_22, n_snap_22);

  if (alpha_11 > threshold_angle || alpha_12 > threshold_angle) {
    if (alpha_11 > alpha_21 && alpha_12 > alpha_22) {
      return true;
    }
  }

  return false;
}

void BoundaryLayerOperation::swapEdgesToMatchShell(vtkUnstructuredGrid *shell_grid, double threshold_angle)
{
  // Set grid to normal*height
  QVector<vec3_t> x_org(m_Grid->GetNumberOfPoints());
  for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
    if (m_BoundaryLayerNode[id_node]) {
      m_Grid->GetPoint(id_node, x_org[id_node].data());
      vec3_t x = x_org[id_node] + m_Height[id_node]*m_BoundaryLayerVectors[id_node];
      m_Grid->GetPoints()->SetPoint(id_node, x.data());
    }
  }

  CgalTriCadInterface cad(shell_grid);
  int num_swaps = 0;
  int count = 0;

  do {
    num_swaps = 0;
    m_Part.setAllCells();
    EG_VTKDCC(vtkIntArray, cell_code, m_Grid, "cell_code");
    EG_VTKDCN(vtkCharArray, node_type, m_Grid, "node_type");
    QVector<bool> swapped(m_Grid->GetNumberOfCells(), false);
    QVector<bool> marked(m_Grid->GetNumberOfCells(), false);
    for (vtkIdType id_cell = 0; id_cell < m_Grid->GetNumberOfCells(); ++id_cell) {
      if (m_BoundaryLayerCodes.contains(cell_code->GetValue(id_cell)) && m_Grid->GetCellType(id_cell) == VTK_TRIANGLE) { //if it is a selected triangle
        if (!marked[id_cell] && !swapped[id_cell]) {
          for (int j = 0; j < 3; ++j) {
            stencil_t stencil = getStencil(id_cell, j);
            if (stencil.id_cell.size() == 2 && stencil.sameBC) {
              if (swapRequired(stencil, &cad, threshold_angle)) {
                marked[stencil.id_cell[0]] = true;
                marked[stencil.id_cell[1]] = true;
                for (int k = 0; k < m_Part.n2cGSize(stencil.id_node[0]); ++k) {
                  vtkIdType id_neigh = m_Part.n2cGG(stencil.id_node[0], k);
                  marked[id_neigh] = true;
                }
                for (int k = 0; k < m_Part.n2cGSize(stencil.id_node[1]); ++k) {
                  vtkIdType id_neigh = m_Part.n2cGG(stencil.id_node[1], k);
                  marked[id_neigh] = true;
                }
                for (int k = 0; k < m_Part.n2cGSize(stencil.p1); ++k) {
                  vtkIdType id_neigh = m_Part.n2cGG(stencil.p1, k);
                  marked[id_neigh] = true;
                }
                for (int k = 0; k < m_Part.n2cGSize(stencil.p2); ++k) {
                  vtkIdType id_neigh = m_Part.n2cGG(stencil.p2, k);
                  marked[id_neigh] = true;
                }

                vtkIdType new_pts1[3], new_pts2[3];
                new_pts1[0] = stencil.p1;
                new_pts1[1] = stencil.id_node[1];
                new_pts1[2] = stencil.id_node[0];
                new_pts2[0] = stencil.id_node[1];
                new_pts2[1] = stencil.p2;
                new_pts2[2] = stencil.id_node[0];
                /*
                  vtkIdType old_pts1[3], old_pts2[3];
                  old_pts1[0] = stencil.id_node[0];
                  old_pts1[1] = stencil.p1;
                  old_pts1[2] = stencil.p2;
                  old_pts2[0] = stencil.id_node[1];
                  old_pts2[1] = stencil.p2;
                  old_pts2[2] = stencil.p1;
                  */
                m_Grid->ReplaceCell(stencil.id_cell[0], 3, new_pts1);
                m_Grid->ReplaceCell(stencil.id_cell[1], 3, new_pts2);

                swapped[stencil.id_cell[0]] = true;
                swapped[stencil.id_cell[1]] = true;
                ++num_swaps;
                break;
              }
            }
          }
        }
      }
    }
    ++count;
    cout << count << ": " << num_swaps << endl;
  } while (num_swaps > 0 && count < 3);


  // reset grid to original points
  for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
    if (m_BoundaryLayerNode[id_node]) {
      m_Grid->GetPoints()->SetPoint(id_node, x_org[id_node].data());
    }
  }
}


