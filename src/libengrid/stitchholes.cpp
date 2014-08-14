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
  while (id_node2 != id_start) {
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
      EG_BUG;
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

void StitchHoles::stitchHole(QList<vtkIdType> loop_nodes)
{
  EG_VTKSP(vtkPolyData, edge_pdata);
  EG_VTKSP(vtkPoints, points);
  EG_VTKSP(vtkCellArray, polys);
  for (int i = 0; i < loop_nodes.size(); ++i) {
    vec3_t x;
    m_Grid->GetPoint(loop_nodes[i], x.data());
    x = toPlane(x);
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
  /*
  {
    QString name = GuiMainWindow::pointer()->getCwd() + "/input.vtk";
    EG_VTKSP(vtkPolyDataWriter, vtk);
    vtk->SetFileName(qPrintable(name));
    vtk->SetInputData(edge_pdata);
    vtk->Write();
  }
  */
  triangulate(edge_pdata, tri_grid, m_Bc);
  //writeGrid(tri_grid, "tri");
  gridFromPlane(tri_grid);

  EG_FORALL_CELLS(id_cell, tri_grid) {
    if (cellNormal(tri_grid, id_cell)*m_N < 0) {
      vtkIdType num_pts, *pts;
      tri_grid->GetCellPoints(id_cell, num_pts, pts);
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
  /*
  int n1 = tri_grid->GetNumberOfPoints();
  int n2 = tri_grid->GetNumberOfCells();
  writeGrid(m_Grid, "before");
  */
  m_Part.addPartition(tri_part);
  /*
  int N1 = m_Grid->GetNumberOfPoints();
  int N2 = m_Grid->GetNumberOfCells();
  writeGrid(m_Grid, "after");
  */
}

void StitchHoles::operate()
{
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
  } while (hole_found);
  cout << "  finished" << endl;
}
