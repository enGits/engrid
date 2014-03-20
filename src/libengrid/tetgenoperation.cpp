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
#include "optimisenormalvector.h"

#include <QDir>

TetGenOperation::TetGenOperation()
{
  // The following code will be used once the option to start TetGen in a
  // separate process has been implemented

  /*
  QDir cwd(GuiMainWindow::pointer()->getCwd());
  m_TetGenPath = GuiMainWindow::pointer()->getCwd() + "tetgen";
  QDir tetgen_dir(m_TetGenPath);
  if (!tetgen_dir.exists()) {
    cwd.mkdir("tetgen");
  }
  */
}

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
    tgio.numberoftetrahedra = tetrahedra.size();
    tgio.tetrahedronlist = new int [tetrahedra.size()*4];
    int i = 0;
    foreach (vtkIdType id_cell, tetrahedra) {
      vtkIdType num_pts, *pts;
      m_Grid->GetCellPoints(id_cell, num_pts, pts);
      for (int j = 0; j < num_pts; ++j) {
        tgio.tetrahedronlist[i*num_pts + j] = pts[j];
      }
      ++i;
    }
    tgio.numberofpointmtrs = 1;
    tgio.pointmtrlist = new REAL [tgio.numberofpoints];
    for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
      vec3_t x;
      m_Grid->GetPoint(id_node, x.data());
      double l = min(m_MaximalEdgeLength, m_ELSManager.minEdgeLength(x));
      tgio.pointmtrlist[id_node] = l;
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
    in >> m_2dFeatureResolution;
    in >> m_3dFeatureResolution;
  }

  m_ELSManager.read();
  m_ELSManager.readRules(m_Grid);
}

QString TetGenOperation::qualityText()
{
  using namespace GeometryTools;
  vec3_t x1(1.0/sqrt(3.0), 0, 0);
  vec3_t x2 = rotate(x1, vec3_t(0, 0, 1), deg2rad(60.0));
  vec3_t x3 = rotate(x2, vec3_t(0, 0, 1), deg2rad(60.0));
  vec3_t x4(0, 0, 1.0);
  vec3_t n1 = triNormal(x1, x2, x3);
  vec3_t n2 = triNormal(x1, x4, x2);
  vec3_t n3 = triNormal(x2, x4, x3);
  double alpha = angle(n1, n2);
  alpha = min(alpha, angle(n1, n3));
  alpha = min(alpha, angle(n2, n3));
  QString alpha_txt;
  alpha_txt.setNum(rad2deg(alpha));
  return QString("1.4") + "/" + alpha_txt;
}

void TetGenOperation::tetgen(QString flags)
{
  try {
    tetgenio in, out, background;
    readSettings();
    copyToTetGen(in);
    copyToTetGen(background);
    char *char_flags = new char [flags.size() + 1];
    strcpy(char_flags, qPrintable(flags));
    tetrahedralize(char_flags, &in, &out, NULL, &background);
    copyFromTetGen(out);
    delete [] char_flags;
  } catch (int x) {
    switch (x) {
    case 1: // Out of memory.
      EG_ERR_RETURN("TetGen error:  Out of memory");
      break;
    case 2: // Encounter an internal error.
      EG_ERR_RETURN("TetGen error:  internal error");
      break;
    case 3:
      EG_ERR_RETURN("TetGen error:  A self-intersection was detected.");
      break;
    case 4:
      EG_ERR_RETURN("TetGen error:  A very small input feature size was detected.");
      break;
    case 5:
      EG_ERR_RETURN("TetGen error:  Two very close input facets were detected.");
      break;
    case 10:
      EG_ERR_RETURN("TetGen error:  An input error was detected.");
      EG_BUG;
      break;
    default:
      EG_ERR_RETURN("TetGen error:  unknown error");
      break;
    }
  }
}

