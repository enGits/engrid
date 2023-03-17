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

#include "stitchholes.h"
#include "guimainwindow.h"
#include "deletestraynodes.h"

#include <vtkCellArray.h>
#include <vtkPolyDataWriter.h>

StitchHoles::StitchHoles(int bc)
{
  m_Bc = bc;
  m_Cad = GuiMainWindow::pointer()->getCadInterface(bc);
}

QList<vtkIdType> StitchHoles::getNextHole()
{
  QList<vtkIdType> loop_nodes;
  EG_VTKDCC(vtkIntArray, cell_code, m_Grid, "cell_code");

  // find the first two nodes of the hole
  //
  vtkIdType id_node1 = -1;
  vtkIdType id_node2 = -1;
  EG_FORALL_CELLS (id_cell, m_Grid) {
    if (isSurface(id_cell, m_Grid) && cell_code->GetValue(id_cell) == m_Bc) {
      for (int i = 0; i < m_Part.c2cGSize(id_cell); ++i) {
        if (m_Part.c2cGG(id_cell, i) == -1) {
          QList<vtkIdType> pts;
          getPointsOfCell(m_Grid, id_cell, pts);
          pts << pts.first();
          id_node1 = pts[i+1];
          id_node2 = pts[i];
          break;
        }
      }
      if (id_node1 != -1) {
        break;
      }
    }
  }
  if (id_node1 == -1) {
    return loop_nodes;
  }

  // create node loop around hole
  //
  vtkIdType id_start = id_node1;
  loop_nodes << id_node1;
  vec3_t x0;
  m_Grid->GetPoint(id_node1, x0.data());
  while (id_node2 != id_start && loop_nodes.size() < m_Grid->GetNumberOfPoints()) {
    loop_nodes << id_node2;
    vec3_t x;
    m_Grid->GetPoint(id_node2, x.data());
    x0 += x;
    bool found = false;
    for (int i = 0; i < m_Part.n2nGSize(id_node2); ++i) {
      vtkIdType id_node3 = m_Part.n2nGG(id_node2, i);
      if (id_node3 != id_node1) {
        QList<vtkIdType> edge_faces;
        m_Part.getEdgeFaces(id_node2, id_node3, edge_faces);
        if (edge_faces.size() == 1) {
          found = true;
          id_node1 = id_node2;
          id_node2 = id_node3;
          break;
        }
      }
    }
    if (!found) {
      loop_nodes.clear();
      return loop_nodes;
    }
  }
  x0 *= 1.0/loop_nodes.size();
  m_Cad->snap(x0);
  vec3_t n = m_Cad->getLastNormal();
  setOrigin(x0);
  setNormal(n);
  setupTransformation();
  return loop_nodes;
}

vec3_t StitchHoles::transformFromPlane(vec3_t x)
{
  vec3_t x_best = x;
  double dist_min = EG_LARGE_REAL;
  for (int i = 0; i < m_X2.size(); ++i) {
    double dist = (x - m_X2[i]).abs();
    if (dist < dist_min) {
      dist_min = dist;
      x_best = m_X3[i];
    }
  }
  return x_best;
}

void StitchHoles::stitchHole(QList<vtkIdType> loop_nodes)
{
  EG_VTKSP(vtkPolyData, edge_pdata);
  EG_VTKSP(vtkPoints, points);
  EG_VTKSP(vtkCellArray, polys);

  m_X2.clear();
  m_X3.clear();

  for (int i = 0; i < loop_nodes.size(); ++i) {
    vec3_t x;
    m_Grid->GetPoint(loop_nodes[i], x.data());
    m_X3 << x;
    x = toPlane(x);
    m_X2 << x;
    points->InsertNextPoint(x.data());
  }
  EG_VTKSP(vtkIdList, pts);
  pts->SetNumberOfIds(loop_nodes.size());
  for (vtkIdType i = 0; i < pts->GetNumberOfIds(); ++i) {
    pts->SetId(i, i);
  }
  polys->InsertNextCell(pts);
  edge_pdata->SetPoints(points);
  edge_pdata->SetPolys(polys);
  EG_VTKSP(vtkUnstructuredGrid, tri_grid);
  triangulate(edge_pdata, tri_grid, m_Bc);
  for (vtkIdType id_node = 0; id_node < tri_grid->GetNumberOfPoints(); ++id_node) {
    vec3_t x;
    tri_grid->GetPoint(id_node, x.data());
    x = transformFromPlane(x);
    tri_grid->GetPoints()->SetPoint(id_node, x.data());
  }


  EG_FORALL_CELLS(id_cell, tri_grid) {
    if (cellNormal(tri_grid, id_cell)*m_N < 0) {
      EG_GET_CELL(id_cell, tri_grid);
      QVector<vtkIdType> nodes(num_pts);
      for (vtkIdType j = 0; j < num_pts; ++j) {
        nodes[j] = pts[j];
      }
      for (vtkIdType j = 0; j < num_pts; ++j) {
        pts[num_pts - j - 1] = nodes[j];
      }
    }
  }

  MeshPartition tri_part(tri_grid, true);
  m_Part.addPartition(tri_part);

  DeleteStrayNodes del_stray;
  del_stray.setGrid(m_Grid);
  del_stray.setAllCells();
  del_stray();

  ++m_Counter;
  QString counter_txt;
  counter_txt.setNum(m_Counter);
  counter_txt = counter_txt.rightJustified(3, '0');
}

void StitchHoles::operate()
{
  m_Counter = 0;
  bool hole_found = false;
  cout << "stitching holes of \"" << qPrintable(GuiMainWindow::pointer()->getBC(m_Bc).getName()) << "\"" << endl;
  int count = 0;
  m_Part.trackGrid(m_Grid);
  do {
    QList<vtkIdType> loop_nodes = getNextHole();
    hole_found = loop_nodes.size() > 0;
    if (hole_found) {
      stitchHole(loop_nodes);
      ++count;
      cout << "  " << count << ". hole stiched" << endl;
    }
  } while (hole_found && count < 50);
  cout << "  finished" << endl;
}
