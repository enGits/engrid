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
#include "tetgenoperation.h"
#include "guimainwindow.h"

void TetGenOperation::copyToTetGen(tetgenio &tgio)
{
  EG_VTKDCC(vtkIntArray, cell_code, m_Grid, "cell_code");
  EG_VTKDCC(vtkIntArray, cell_orgdir, m_Grid, "cell_orgdir");
  EG_VTKDCC(vtkIntArray, cell_curdir, m_Grid, "cell_curdir");
  EG_VTKDCC(vtkIntArray, cell_voldir, m_Grid, "cell_voldir");
  tgio.initialize();
  tgio.numberofpoints = m_Grid->GetNumberOfPoints();
  tgio.pointlist = new REAL [3*tgio.numberofpoints];
  for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
    vec3_t x;
    m_Grid->GetPoint(id_node, x.data());
    for (int i = 0; i < 3; ++i) {
      tgio.pointlist[3*id_node + i] = x[i];
    }
  }
  QList<vtkIdType> triangles, tetrahedra;
  for (vtkIdType id_cell = 0; id_cell < m_Grid->GetNumberOfCells(); ++id_cell) {
    if      (m_Grid->GetCellType(id_cell) == VTK_TRIANGLE) triangles.append(id_cell);
    else if (m_Grid->GetCellType(id_cell) == VTK_TETRA)    tetrahedra.append(id_cell);
    else {
      EG_BUG;
    }
  }
  tgio.numberoffacets = triangles.size();
  tgio.facetmarkerlist = new int [triangles.size()];
  tgio.facetlist = new tetgenio::facet [triangles.size()];
  tgio.numberoftetrahedra = tetrahedra.size();
  tgio.tetrahedronlist = new int [tetrahedra.size()*4];
  {
    int i = 0;
    foreach (vtkIdType id_cell, triangles) {
      m_OrgDir = cell_orgdir->GetValue(id_cell);
      m_CurDir = cell_curdir->GetValue(id_cell);
      m_VolDir = cell_voldir->GetValue(id_cell);
      vtkIdType num_pts, *pts;
      m_Grid->GetCellPoints(id_cell, num_pts, pts);
      tgio.facetlist[i].numberofpolygons = 1;
      tgio.facetlist[i].polygonlist = new tetgenio::polygon[1];
      tgio.facetlist[i].numberofholes = 0;
      tgio.facetlist[i].holelist = NULL;
      tetgenio::polygon &poly = tgio.facetlist[i].polygonlist[0];
      poly.numberofvertices = 3;
      poly.vertexlist = new int[3];
      for (int j = 0; j < num_pts; ++j) {
        poly.vertexlist[j] = pts[j];
      }
      tgio.facetmarkerlist[i] = cell_code->GetValue(id_cell);
      ++i;
    }
  }
  {
    int i = 0;
    foreach (vtkIdType id_cell, tetrahedra) {
      vtkIdType num_pts, *pts;
      m_Grid->GetCellPoints(id_cell, num_pts, pts);
      for (int j = 0; j < num_pts; ++j) {
        tgio.tetrahedronlist[i*num_pts + j] = pts[j];
      }
      ++i;
    }
  }
}

void TetGenOperation::copyFromTetGen(tetgenio &tgio)
{
  EG_VTKSP(vtkUnstructuredGrid, new_grid);
  allocateGrid(new_grid, tgio.numberoftetrahedra + tgio.numberoftrifaces, tgio.numberofpoints);
  EG_VTKDCC(vtkIntArray, cell_code, new_grid, "cell_code");
  EG_VTKDCC(vtkIntArray, cell_orgdir, new_grid, "cell_orgdir");
  EG_VTKDCC(vtkIntArray, cell_curdir, new_grid, "cell_curdir");
  EG_VTKDCC(vtkIntArray, cell_voldir, new_grid, "cell_voldir");
  for (vtkIdType id_node = 0; id_node < tgio.numberofpoints; ++id_node) {
    vec3_t x;
    for (int i = 0; i < 3; ++i) {
      x[i] = tgio.pointlist[3*id_node + i];
    }
    new_grid->GetPoints()->SetPoint(id_node, x.data());
  }
  for (int i = 0; i < tgio.numberoftetrahedra; ++i) {
    vtkIdType num_pts = 4, pts[4];
    for (int j = 0; j < num_pts; ++j) {
      pts[j] = tgio.tetrahedronlist[i*num_pts + j];
    }
    vtkIdType id_cell = new_grid->InsertNextCell(VTK_TETRA, num_pts, pts);
    cell_code->SetValue(id_cell, 0);
  }
  for (int i = 0; i < tgio.numberoftrifaces; ++i) {
    vtkIdType num_pts = 3, pts[3];
    for (int j = 0; j < num_pts; ++j) {
      pts[2-j] = tgio.trifacelist[i*num_pts + j];
    }
    vtkIdType id_cell = new_grid->InsertNextCell(VTK_TRIANGLE, num_pts, pts);
    cell_code->SetValue(id_cell, tgio.trifacemarkerlist[i]);
    cell_orgdir->SetValue(id_cell, m_OrgDir);
    cell_curdir->SetValue(id_cell, m_CurDir);
    cell_voldir->SetValue(id_cell, m_VolDir);
  }
  makeCopy(new_grid, m_Grid);
  //correctSurfaceOrientation();
}


void TetGenOperation::readSettings()
{
  QString buffer = GuiMainWindow::pointer()->getXmlSection("engrid/surface/settings");
  if(!buffer.isEmpty()) {
    QTextStream in(&buffer, QIODevice::ReadOnly);
    in >> m_MaximalEdgeLength;
    in >> m_MinimalEdgeLength;
    in >> m_GrowthFactor;
    in >> m_NodesPerQuarterCircle;
    int num_bcs;
    in >> num_bcs;
    int check_state;
    for (int i = 0; i < num_bcs; ++i) {
      in >> check_state;
    }
    if (!in.atEnd()) {
      in >> m_2dFeatureResolution;
      in >> m_3dFeatureResolution;
    }
  }
  m_ELSManager.read();
}

void TetGenOperation::correctSurfaceOrientation()
{
  m_Part.setGrid(m_Grid);
  m_Part.setAllCells();
  for (vtkIdType id_cell = 0; id_cell < m_Grid->GetNumberOfCells(); ++id_cell) {
    if (isVolume(id_cell, m_Grid)) {
      for (int i = 0; i < m_Part.c2cGSize(id_cell); ++i) {
        vtkIdType id_face = m_Part.c2cGG(id_cell, i);
        if (isSurface(id_face, m_Grid)) {
          vec3_t n = cellNormal(m_Grid, id_face);
          vec3_t v = cellCentre(m_Grid, id_cell) - cellCentre(m_Grid, id_face);
          if (n*v < 0) {
            vtkIdType num_pts, *pts;
            m_Grid->GetCellPoints(id_face, num_pts, pts);
            swap(pts[0], pts[1]);
          }
        }
      }
    }
  }
}

void TetGenOperation::tetgen(QString flags)
{
  tetgenio in, out;
  copyToTetGen(in);
  char *char_flags = new char [flags.size() + 1];
  strcpy(char_flags, qPrintable(flags));
  tetrahedralize(char_flags, &in, &out);
  copyFromTetGen(out);
  delete [] char_flags;
}

