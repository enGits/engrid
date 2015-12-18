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

#include "vtkEgNormalExtrusion.h"
#include "geometrytools.h"

vtkStandardNewMacro(vtkEgNormalExtrusion)

vtkEgNormalExtrusion::vtkEgNormalExtrusion()
{
  m_Mode = normal;
  m_ScaleX = 1;
  m_ScaleY = 1;
  m_ScaleZ = 1;
  m_RemoveInternalFaces = true;
  m_Curve1 = NULL;
  m_Curve2 = NULL;
  m_OrthoExtrusion = false;
}

int vtkEgNormalExtrusion::getNewBc(MeshPartition *part, vtkIdType id_node1, vtkIdType id_node2)
{
  QVector<vtkIdType> nodes(2);
  nodes[0] = id_node1;
  nodes[1] = id_node2;
  QSet<int> common_bcs = part->getCommonBcs(nodes);
  QSet<int> potential_bcs = common_bcs - m_BoundaryCodes;
  if (potential_bcs.size() > 0) {
    return potential_bcs.values().first();
  }
  return m_UnknownBc;
}

void vtkEgNormalExtrusion::SetCurves(Curve *curve1, Curve *curve2)
{
  m_Curve1 = curve1;
  m_Curve2 = curve2;
  if (m_Curve1 && m_Curve2) {
    m_Mode = curve;
  } else {
    m_Mode = normal;
  }
}

bool vtkEgNormalExtrusion::PrepareCurveExtrusion(QVector<vtkIdType> nodes)
{
  if (!m_Curve1 || !m_Curve2) {
    return false;
  }
  mat3_t A;
  if (m_OrthoExtrusion) {
    A = m_Curve1->computeOrthoBase(0, m_Curve2);
  } else {
    A = m_Curve2->computeBase(0, m_Curve2);
  }
  A = A.inverse();
  m_InitialBaseCoords.resize(nodes.size());
  vec3_t x0 = m_Curve1->position(0);
  for (int i = 0; i < nodes.size(); ++i) {
    vec3_t x;
    m_Input->GetPoint(nodes[i], x.data());
    x -= x0;
    m_InitialBaseCoords[i] = A*x;
  }
  return true;
}

mat3_t vtkEgNormalExtrusion::ComputeBase(double l)
{
  mat3_t A;
  if (m_OrthoExtrusion) {
    A = m_Curve1->computeOrthoBase(l, m_Curve2);
  } else {
    A = m_Curve2->computeBase(l, m_Curve2);
  }
  return A;
}

void vtkEgNormalExtrusion::ExecuteEg()
{
  QVector<vtkIdType> cells, nodes, n1, n2;
  QVector<vec3_t> cell_normals, node_normals;
  ExtractBoundary(cells, nodes, m_BoundaryCodes, m_Input);
  if (m_Mode == normal) {
    computeNormals(cell_normals, node_normals, cells, nodes,m_Input);
  } else if (m_Mode == cylinder) {
    axis.normalise();
    node_normals.resize(nodes.size());
    for (int i = 0; i < node_normals.size(); ++i) {
      vec3_t x;
      m_Input->GetPoint(nodes[i],x.data());
      vec3_t x0 = x - ((x-origin)*axis)*axis;
      node_normals[i] = x0 - origin;
      node_normals[i].normalise();
    }
  } else if ((m_Mode == fixed) || (m_Mode == planar)) {
    fixed_normal.normalise();
    node_normals.resize(nodes.size());
    for (int i = 0; i < node_normals.size(); ++i) {
      node_normals[i] = fixed_normal;
    }
  }
  for (int i = 0; i < node_normals.size(); ++ i) {
    node_normals[i][0] *= m_ScaleX;
    node_normals[i][1] *= m_ScaleY;
    node_normals[i][2] *= m_ScaleZ;
  }
  n1.resize(nodes.size());
  n2.resize(nodes.size());
  if (m_Mode == curve) {
    if (!PrepareCurveExtrusion(nodes)) {
      return;
    }
  }

  // mapping
  QVector<int> _cells, _nodes;
  createNodeMapping(nodes, _nodes, m_Input);
  createCellMapping(cells, _cells, m_Input);
  QVector<QSet<int> > n2c;
  createNodeToCell(cells, nodes, _nodes, n2c, m_Input);

  vtkIdType NnewNodes = m_Input->GetNumberOfPoints() + (layer_y.size()-1)*nodes.size();
  vtkIdType NnewCells = m_Input->GetNumberOfCells() + (layer_y.size()-1)*cells.size();

  // count the number of new surface elements (side walls)
  for (int i_cell = 0; i_cell < cells.size(); ++i_cell) {
    vtkIdType *pts;
    vtkIdType  Npts;
    m_Input->GetCellPoints(cells[i_cell], Npts, pts);
    QVector<vtkIdType> surf_pts(Npts);
    for (int i_pts = 0; i_pts < Npts; ++i_pts) {
      surf_pts[i_pts] = _nodes[pts[i_pts]];
    }
    for (int i_surf_pts = 0; i_surf_pts < Npts; ++i_surf_pts) {
      vtkIdType p1 = surf_pts[i_surf_pts];
      vtkIdType p2;
      if (i_surf_pts < Npts - 1) {
        p2 = surf_pts[i_surf_pts + 1];
      } else {
        p2 = surf_pts[0];
      }
      bool add_bd = false;
      {
        int N = 0;
        int cell;
        foreach(cell, n2c[p1]) {
          if (n2c[p2].contains(cell)) ++N;
        }
        if (N == 1) add_bd = true;
      }
      if (add_bd) {
        NnewCells += layer_y.size();
      }
    }
  }

  // count the number of new surface elements (base)
  QVector<bool> is_boundary;
  is_boundary.fill(false, cells.size());
  {
    QVector<int> nvol;
    nvol.fill(0, nodes.size());
    for (vtkIdType id_cell = 0; id_cell < m_Input->GetNumberOfCells(); ++id_cell) {
      if (isVolume(id_cell, m_Input)) {
        vtkIdType Npts, *pts;
        m_Input->GetCellPoints(id_cell, Npts, pts);
        for (int i = 0; i < Npts; ++i) {
          if (_nodes[pts[i]] >= 0) {
            ++nvol[_nodes[pts[i]]];
          }
        }
      }
    }
    for (int i_cell = 0; i_cell < cells.size(); ++i_cell) {
      vtkIdType id_cell = cells[i_cell];
      vtkIdType Npts, *pts;
      m_Input->GetCellPoints(id_cell, Npts, pts);
      for (int i = 0; i < Npts; ++i) {
        if (nvol[_nodes[pts[i]]] == 0) {
          is_boundary[i_cell] = true;
        }
      }
      if (is_boundary[i_cell] || !m_RemoveInternalFaces) {
        ++NnewCells;
      }
    }
  }

  // allocate memory for the new grid
  allocateGrid(m_Output, NnewCells, NnewNodes);

  // boundary conditions
  EG_VTKDCC(vtkIntArray, cell_code1, m_Input, "cell_code");
  EG_VTKDCC(vtkIntArray, cell_code2, m_Output, "cell_code");
  EG_VTKDCC(vtkIntArray, orgdir, m_Output, "cell_orgdir");
  EG_VTKDCC(vtkIntArray, voldir, m_Output, "cell_voldir");
  EG_VTKDCC(vtkIntArray, curdir, m_Output, "cell_curdir");

  m_UnknownBc = 1;
  for (vtkIdType id_cell = 0; id_cell < m_Input->GetNumberOfCells(); ++id_cell) {
    if (isSurface(id_cell, m_Input)) {
      if (cell_code1->GetValue(id_cell) >= m_UnknownBc) {
        m_UnknownBc = cell_code1->GetValue(id_cell) + 1;
      }
    }
  }

  QMap<QPair<int,int>,int> new_bcs;

  MeshPartition part(m_Input, true);

  for (int i = 0; i < nodes.size(); ++i) {
    n2[i] = nodes[i];
    vec3_t x;
    m_Input->GetPoint(nodes[i],x.data());
    m_Output->GetPoints()->SetPoint(n2[i],x.data());
  }
  double total_layers = layer_y[layer_y.size()-1] - layer_y[0];
  QVector<double> total_dist(nodes.size());
  double L_max = -1e99;
  int i_max = 0;
  if (m_Mode == planar) {
    for (int i = 0; i < nodes.size(); ++i) {
      total_dist[i] = total_layers;
      vec3_t x_origin, x_target;
      m_Input->GetPoint(nodes[i],x_origin.data());
      x_target = x_origin + total_layers*node_normals[i];
      double L = x_target*fixed_normal;
      if (L > L_max) {
        i_max = i;
        L_max = L;
      }
    }
    vec3_t x_far;
    m_Input->GetPoint(nodes[i_max],x_far.data());
    x_far += min_dist*fixed_normal;
    for (int i = 0; i < nodes.size(); ++i) {
      total_dist[i] = total_layers;
      if (m_Mode == planar) {
        vec3_t x_origin;
        m_Input->GetPoint(nodes[i],x_origin.data());
        total_dist[i] = (x_far-x_origin)*fixed_normal;
      }
    }
  }
  for (int i_layer = 0; i_layer < layer_y.size() - 1; ++i_layer) {
    mat3_t current_base;
    vec3_t current_origin;
    if (m_Mode == curve) {
      current_origin = m_Curve1->position(layer_y[i_layer + 1]);
      current_base = ComputeBase(layer_y[i_layer + 1]);
    }
    for (int i = 0; i < n1.size(); ++i) {
      n1[i] = n2[i];
      n2[i] = i_layer*nodes.size() + i + m_Input->GetNumberOfPoints();
      vec3_t x1, x2;
      m_Output->GetPoint(n1[i],x1.data());
      if (m_Mode == rotation) {

        x2 = x1 - origin;
        double alpha = (layer_y[i_layer + 1] - layer_y[i_layer])*M_PI/180.0;
        x2 = GeometryTools::rotate(x2,axis,alpha);
        x2 += origin;

      } else if (m_Mode == curve) {

        x2 = current_origin + current_base*m_InitialBaseCoords[i];

      } else {

        double dist = (layer_y[i_layer + 1] - layer_y[i_layer]);
        if (m_Mode == planar) {
          dist *= total_dist[i]/total_layers;
        }
        x2 = x1 + dist*node_normals[i];

      }
      m_Output->GetPoints()->SetPoint(n2[i],x2.data());
    }

    for (int i_cell = 0; i_cell < cells.size(); ++i_cell) {
      vtkIdType *pts;
      vtkIdType  Npts;
      m_Input->GetCellPoints(cells[i_cell], Npts, pts);
      QVector<vtkIdType> surf_pts(Npts);
      for (int i_pts = 0; i_pts < Npts; ++i_pts) {
        surf_pts[i_pts] = _nodes[pts[i_pts]];
      }
      for (int i_surf_pts = 0; i_surf_pts < Npts; ++i_surf_pts) {
        vtkIdType p1 = surf_pts[i_surf_pts];
        vtkIdType p2;
        if (i_surf_pts < Npts - 1) {
          p2 = surf_pts[i_surf_pts + 1];
        } else {
          p2 = surf_pts[0];
        }
        bool add_bd = false;
        {
          int N = 0;
          int cell;
          foreach(cell, n2c[p1]) {
            if (n2c[p2].contains(cell)) ++N;
          }
          if (N == 1) add_bd = true;
        }
        if (add_bd) {
          int new_bc = new_bcs[QPair<int,int>(p1,p2)];
          if (new_bc == 0) {
            new_bc = getNewBc(&part, nodes[p1], nodes[p2]);
            new_bcs[QPair<int,int>(p1,p2)] = new_bc;
          }
          vtkIdType quad_pts[4];
          quad_pts[0] = n1[p1];
          quad_pts[1] = n1[p2];
          quad_pts[2] = n2[p2];
          quad_pts[3] = n2[p1];
          vtkIdType id_new_cell = m_Output->InsertNextCell(VTK_QUAD,4,quad_pts);
          cell_code2->SetValue(id_new_cell, new_bc);
          orgdir->SetValue(id_new_cell, 0);
          curdir->SetValue(id_new_cell, 0);
          voldir->SetValue(id_new_cell, 0);
        }
      }
      if (Npts == 3) {
        {
          vtkIdType pri_pts[6];
          pri_pts[0] = n1[surf_pts[0]];
          pri_pts[1] = n1[surf_pts[2]];
          pri_pts[2] = n1[surf_pts[1]];
          pri_pts[3] = n2[surf_pts[0]];
          pri_pts[4] = n2[surf_pts[2]];
          pri_pts[5] = n2[surf_pts[1]];
          vtkIdType id_new_cell = m_Output->InsertNextCell(VTK_WEDGE,6,pri_pts);
          cell_code2->SetValue(id_new_cell, 0);
          orgdir->SetValue(id_new_cell, 0);
          curdir->SetValue(id_new_cell, 0);
          voldir->SetValue(id_new_cell, 0);
        }
        if (i_layer == layer_y.size() - 2) {
          vtkIdType tri_pts[3];
          tri_pts[0] = n2[surf_pts[0]];
          tri_pts[1] = n2[surf_pts[1]];
          tri_pts[2] = n2[surf_pts[2]];
          vtkIdType id_new_cell = m_Output->InsertNextCell(VTK_TRIANGLE,3,tri_pts);
          cell_code2->SetValue(id_new_cell, cell_code1->GetValue(cells[i_cell]));
          orgdir->SetValue(id_new_cell, 0);
          curdir->SetValue(id_new_cell, 0);
          voldir->SetValue(id_new_cell, 0);
        }
      }
      if (Npts == 4) {
        {
          vtkIdType pri_pts[8];
          pri_pts[0] = n1[surf_pts[0]];
          pri_pts[1] = n1[surf_pts[1]];
          pri_pts[2] = n1[surf_pts[2]];
          pri_pts[3] = n1[surf_pts[3]];
          pri_pts[4] = n2[surf_pts[0]];
          pri_pts[5] = n2[surf_pts[1]];
          pri_pts[6] = n2[surf_pts[2]];
          pri_pts[7] = n2[surf_pts[3]];
          vtkIdType id_new_cell = m_Output->InsertNextCell(VTK_HEXAHEDRON,8,pri_pts);
          cell_code2->SetValue(id_new_cell, 0);
          orgdir->SetValue(id_new_cell, 0);
          curdir->SetValue(id_new_cell, 0);
          voldir->SetValue(id_new_cell, 0);
        }
        if (i_layer == layer_y.size() - 2) {
          vtkIdType quad_pts[4];
          quad_pts[0] = n2[surf_pts[0]];
          quad_pts[1] = n2[surf_pts[1]];
          quad_pts[2] = n2[surf_pts[2]];
          quad_pts[3] = n2[surf_pts[3]];
          vtkIdType id_new_cell = m_Output->InsertNextCell(VTK_QUAD,4,quad_pts);
          cell_code2->SetValue(id_new_cell, cell_code1->GetValue(cells[i_cell]));
          orgdir->SetValue(id_new_cell, 0);
          curdir->SetValue(id_new_cell, 0);
          voldir->SetValue(id_new_cell, 0);
        }
      }
      if (Npts > 4) {
        {
          EG_VTKSP(vtkIdList, stream);
          stream->SetNumberOfIds(1 + 2*(Npts + 1) + Npts*5);
          vtkIdType id = 0;
          stream->SetId(id++, 2 + Npts);

          // bottom face
          stream->SetId(id++, Npts);
          for (int i = Npts - 1; i >= 0; --i) {
            stream->SetId(id++, n1[surf_pts[i]]);
          }

          // top face
          stream->SetId(id++, Npts);
          for (int i = 0; i < Npts; ++i) {
            stream->SetId(id++, n2[surf_pts[i]]);
          }

          // side faces
          for (int i = 0; i < Npts; ++i) {
            stream->SetId(id++, 4);
            vtkIdType p1 = n1[surf_pts[i]];
            vtkIdType p2 = n1[surf_pts[0]];
            vtkIdType p3 = n2[surf_pts[0]];
            vtkIdType p4 = n2[surf_pts[i]];
            if (i < Npts - 1) {
              p2 = n1[surf_pts[i+1]];
              p3 = n2[surf_pts[i+1]];
            }
            stream->SetId(id++, p1);
            stream->SetId(id++, p2);
            stream->SetId(id++, p3);
            stream->SetId(id++, p4);
          }

          vtkIdType id_new_cell = m_Output->InsertNextCell(VTK_POLYHEDRON, stream);
          cell_code2->SetValue(id_new_cell, 0);
          orgdir->SetValue(id_new_cell, 0);
          curdir->SetValue(id_new_cell, 0);
          voldir->SetValue(id_new_cell, 0);
        }
        if (i_layer == layer_y.size() - 2) {
          EG_VTKSP(vtkIdList, poly_pts);
          poly_pts->SetNumberOfIds(Npts);
          for (vtkIdType id = 0; id < Npts; ++id) {
            poly_pts->SetId(id, n2[surf_pts[id]]);
          }
          vtkIdType id_new_cell = m_Output->InsertNextCell(VTK_POLYGON, poly_pts);
          cell_code2->SetValue(id_new_cell, cell_code1->GetValue(cells[i_cell]));
          orgdir->SetValue(id_new_cell, 0);
          curdir->SetValue(id_new_cell, 0);
          voldir->SetValue(id_new_cell, 0);
        }
      }
    }
  }


  for (vtkIdType nodeId = 0; nodeId < m_Input->GetNumberOfPoints(); ++nodeId) {
    vec3_t x;
    m_Input->GetPoints()->GetPoint(nodeId, x.data());
    m_Output->GetPoints()->SetPoint(nodeId, x.data());
  }

  // copy all original cells that were not part of the extrusion
  for (vtkIdType id_cell = 0; id_cell < m_Input->GetNumberOfCells(); ++id_cell) {
    if (_cells[id_cell] == -1 || !m_RemoveInternalFaces) {
      vtkIdType id_new_cell = copyCell(m_Input, id_cell, m_Output);
      copyCellData(m_Input, id_cell, m_Output, id_new_cell);
    }
  }

  // close boundary where no volume cells were present in the original grid

  for (int i_cell = 0; i_cell < cells.size(); ++i_cell) {
    if (is_boundary[i_cell]) {
      vtkIdType id_cell = cells[i_cell];
      vtkIdType *pts;
      vtkIdType  Npts;
      m_Input->GetCellPoints(id_cell, Npts, pts);
      vtkIdType id_new_cell = m_Output->InsertNextCell(m_Input->GetCellType(id_cell), Npts, pts);
      m_Output->GetCellPoints(id_new_cell, Npts, pts);
      QVector<vtkIdType> nodes(Npts);
      for (vtkIdType j = 0; j < Npts; ++j) nodes[j]          = pts[j];
      for (vtkIdType j = 0; j < Npts; ++j) pts[Npts - j - 1] = nodes[j];
      copyCellData(m_Input, id_cell, m_Output, id_new_cell);
    }
  }

  UpdateCellIndex(m_Output);
}

void vtkEgNormalExtrusion::SetLayers(const QVector<double> &y)
{
  layer_y.resize(y.size());
  qCopy(y.begin(), y.end(), layer_y.begin());
}
