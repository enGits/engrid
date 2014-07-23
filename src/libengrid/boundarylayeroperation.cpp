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

  computeNodeTypes();
  relaxBoundaryLayerVectors();

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
    m_BoundaryLayerVectors[id_node] *= m_Height[id_node];
    if (!checkVector(m_BoundaryLayerVectors[id_node])) {
      EG_ERR_RETURN("invalid layer vector computed");
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

void BoundaryLayerOperation::relaxBoundaryLayerVectors()
{
  m_BoundaryLayerNode.fill(false, m_Grid->GetNumberOfPoints());
  for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
    for (int i = 0; i < m_Part.n2bcGSize(id_node); ++i) {
      if (m_BoundaryLayerCodes.contains(m_Part.n2bcG(id_node, i))) {
        m_BoundaryLayerNode[id_node] = true;
      }
    }
  }
  for (int iter = 0; iter < m_NumBoundaryLayerVectorRelaxations; ++iter) {
    QVector<vec3_t> v_new(m_BoundaryLayerVectors.size(), vec3_t(0,0,0));
    for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
      if (m_BoundaryLayerNode[id_node]) {
        if (m_SnapPoints[id_node].size() > 0) {
          v_new[id_node] = vec3_t(0,0,0);
          foreach (vtkIdType id_snap, m_SnapPoints[id_node]) {
            v_new[id_node] += m_BoundaryLayerVectors[id_snap];
          }
          v_new[id_node].normalise();
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

void BoundaryLayerOperation::writeBoundaryLayerVectors(QString file_name)
{
  QVector<vtkIdType> bcells;
  QSet<int> bcs = getAllBoundaryCodes(m_Grid);
  getSurfaceCells(bcs, bcells, m_Grid);
  MeshPartition bpart(m_Grid);
  bpart.setCells(bcells);
  file_name = GuiMainWindow::pointer()->getCwd() + "/" + file_name + ".vtk";
  QFile file(file_name);
  file.open(QIODevice::WriteOnly);
  QTextStream f(&file);
  f << "# vtk DataFile Version 2.0\n";
  f << "m_BoundaryLayerVectors\n";
  f << "ASCII\n";
  f << "DATASET UNSTRUCTURED_GRID\n";
  f << "POINTS " << bpart.getNumberOfNodes() << " float\n";
  for (int i = 0; i < bpart.getNumberOfNodes(); ++i) {
    vec3_t x;
    vtkIdType id_node = bpart.globalNode(i);
    m_Grid->GetPoint(id_node, x.data());
    f << x[0] << " " << x[1] << " " << x[2] << "\n";
  }
  f << "CELLS " << bpart.getNumberOfCells() << " ";
  int N = 0;
  for (int i = 0; i < bpart.getNumberOfCells(); ++ i) {
    vtkIdType id_cell = bpart.globalCell(i);
    vtkIdType N_pts, *pts;
    m_Grid->GetCellPoints(id_cell, N_pts, pts);
    N += 1 + N_pts;
  }
  f << N << "\n";
  for (int i = 0; i < bpart.getNumberOfCells(); ++ i) {
    vtkIdType id_cell = bpart.globalCell(i);
    vtkIdType N_pts, *pts;
    m_Grid->GetCellPoints(id_cell, N_pts, pts);
    f << N_pts;
    for (int j = 0; j < N_pts; ++j) {
      f << " " << bpart.localNode(pts[j]);
    }
    f << "\n";
  }

  f << "CELL_TYPES " << bpart.getNumberOfCells() << "\n";
  for (int i = 0; i < bpart.getNumberOfCells(); ++ i) {
    vtkIdType id_cell = bpart.globalCell(i);
    vtkIdType type_cell = m_Grid->GetCellType(id_cell);
    f << type_cell << "\n";
  }
  f << "POINT_DATA " << bpart.getNumberOfNodes() << "\n";

  f << "VECTORS N float\n";
  for (int i = 0; i < bpart.getNumberOfNodes(); ++i) {
    vtkIdType id_node = bpart.globalNode(i);
    f << m_BoundaryLayerVectors[id_node][0] << " ";
    f << m_BoundaryLayerVectors[id_node][1] << " ";
    f << m_BoundaryLayerVectors[id_node][2] << "\n";
  }

  f << "SCALARS node_type int\n";
  f << "LOOKUP_TABLE default\n";
  for (int i = 0; i < bpart.getNumberOfNodes(); ++i) {
    vtkIdType id_node = bpart.globalNode(i);
    f << m_NodeTypes[id_node] << "\n";
  }
}

void BoundaryLayerOperation:: computeDesiredHeights()
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
      double h1 = m_ELSManagerSurface.minEdgeLength(x)*m_FarfieldRatio;
      k = max(k, int(logarithm(m_StretchingRatio, h1/h0)));
    }
  }
  for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
    if (m_BoundaryLayerNode[id_node]) {
      vec3_t x;
      m_Grid->GetPoint(id_node, x.data());
      double h0 = m_ELSManagerBLayer.minEdgeLength(x);
      double h1 = cl->GetValue(id_node);
      double s  = pow(h1/h0, 1.0/k);
      double H  = h0*(1 - pow(s, k))/(1 - s);
      if (!checkReal(H)) {
        EG_ERR_RETURN("floating point error while computing heights");
      }
      if (H < 0) {
        EG_ERR_RETURN("negative height computed");
      }
      if (H > 1000*h1) {
        EG_ERR_RETURN("unrealistically large height computed");
      }
      if (H < 1e-3*h0) {
        EG_ERR_RETURN("unrealistically small height computed");
      }
      if (h1 < h0) {
        EG_ERR_RETURN("h1 < h0");
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

  // gaps
  if (true) {

    // check for potential collisions
    CgalTriCadInterface cad(m_Grid);

    for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
      if (m_BoundaryLayerNode[id_node]) {
        vec3_t x;
        m_Grid->GetPoint(id_node, x.data());
        QVector<QPair<vec3_t, vtkIdType> > intersections;
        cad.computeIntersections(x, m_BoundaryLayerVectors[id_node], intersections);
        vec3_t n = m_Part.globalNormal(id_node);
        for (int i = 0; i < intersections.size(); ++i) {
          vec3_t xi = intersections[i].first;
          vtkIdType id_tri = intersections[i].second;
          vec3_t ni = GeometryTools::cellNormal(m_Grid, id_tri);
          ni.normalise();
          vec3_t dx = xi - x;
          if (dx*n < 0 && ni*n < 0) {
            m_Height[id_node] = min(m_Height[id_node], m_MaxHeightInGaps*dx.abs());
          }
        }
      }
    }

  }

  // limit face size difference (neighbour collisions)
  if (true) {
    bool done;
    int iter = 0;
    do {
      done = true;
      QVector<double> new_scale(m_Grid->GetNumberOfPoints(), 1.0);
      QVector<int> count(m_Grid->GetNumberOfPoints(), 1);
      for (vtkIdType id_cell = 0; id_cell < m_Grid->GetNumberOfCells(); ++id_cell) {
        if (isSurface(id_cell, m_Grid)) {
          vtkIdType num_pts, *pts;
          m_Grid->GetCellPoints(id_cell, num_pts, pts);
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

              /*
              while (scale2 - scale1 > 1e-3) {
                if (faceFine(id_cell, 0.5*(scale1 + scale2))) {
                  scale1 = 0.5*(scale1 + scale2);
                } else {
                  scale2 = 0.5*(scale1 + scale2);
                }
              }
              */

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

  // smoothing
  for (int iter = 0; iter < m_NumBoundaryLayerHeightRelaxations; ++iter) {
    QVector<double> h_new = m_Height;
    for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
      if (m_BoundaryLayerNode[id_node]) {
        if (m_SnapPoints[id_node].size() > 0) {
          h_new[id_node] = 0;
          if (m_SnapPoints[id_node].size() == 2) {
            h_new[id_node] += m_Height[id_node];
          }
          foreach (vtkIdType id_snap, m_SnapPoints[id_node]) {
            h_new[id_node] += m_Height[id_snap];
          }
          if (m_SnapPoints[id_node].size() == 2) {
            h_new[id_node] /= 3;
          } else {
            h_new[id_node] /= m_SnapPoints[id_node].size();
          }
        } else {
          h_new[id_node] = m_Height[id_node];
        }
      }
    }
    m_Height = h_new;
  }

  cout << "heights computed" << endl;

}


