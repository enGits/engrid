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

#include "polymolecule.h"
#include "guimainwindow.h"

#include <vtkUnstructuredGridWriter.h>
#include <vtkGeometryFilter.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkDelaunay3D.h>
#include <vtkHull.h>
#include <vtkTriangleFilter.h>

PolyMolecule::PolyMolecule()
{
  m_AllPositive = false;
  m_PolyMesh = NULL;
  m_SplitCell1 = NULL;
  m_SplitCell2 = NULL;
}

PolyMolecule::PolyMolecule(PolyMesh *poly_mesh, int i_pcell)
{
  m_PCell = i_pcell;
  QVector<QList<int> > faces(poly_mesh->numFacesOfPCell(i_pcell));
  for (int i = 0; i < poly_mesh->numFacesOfPCell(i_pcell); ++i) {
    int i_face = poly_mesh->pcell2Face(i_pcell, i);
    int num_nodes = poly_mesh->numNodes(i_face);
    for (int j = 0; j < num_nodes; ++j) {
      int node = poly_mesh->nodeIndex(i_face, j);
      if (poly_mesh->owner(i_face) == i_pcell) {
        faces[i].append(node);
      } else {
        faces[i].prepend(node);
      }
    }
  }
  m_AllPositive = false;
  m_SplitCell1 = NULL;
  m_SplitCell2 = NULL;
  init(poly_mesh, faces);
}

void PolyMolecule::writeVtkFile(QString file_name)
{
  EG_VTKSP(vtkUnstructuredGrid, grid);
  allocateGrid(grid, m_Faces.size(), m_Nodes.size());

  EG_VTKSP(vtkDoubleArray, normals);
  normals->SetNumberOfComponents(3);
  normals->SetNumberOfTuples(m_Nodes.size());
  normals->SetName("N");

  EG_VTKSP(vtkDoubleArray, badness);
  badness->SetNumberOfComponents(1);
  badness->SetNumberOfTuples(m_Nodes.size());
  badness->SetName("badness");

  for (int i = 0; i < m_Nodes.size(); ++i) {
    vec3_t x = getXNode(i);
    grid->GetPoints()->SetPoint(i, x.data());
    normals->SetTuple(i, m_NodeNormals[i].data());
    badness->SetValue(i, double(m_N2BadFaces[i].size()));
  }
  grid->GetPointData()->AddArray(normals);
  grid->GetPointData()->AddArray(badness);

  for (int i = 0; i < m_Faces.size(); ++i) {
    QVector<vtkIdType> pts(m_Faces[i].size());
    for (int j = 0; j < m_Faces[i].size(); ++j) {
      pts[j] = m_Faces[i][j];
    }
    grid->InsertNextCell(VTK_POLYGON, m_Faces[i].size(), pts.data());
  }
  EG_VTKSP(vtkXMLUnstructuredGridWriter, vtu);
  file_name = GuiMainWindow::pointer()->getCwd() + "/" + file_name + ".vtu";
  vtu->SetFileName(qPrintable(file_name));
  vtu->SetDataModeToBinary();
  vtu->SetInput(grid);
  vtu->Write();
}

void PolyMolecule::createPolyData(vtkPolyData *poly_data)
{
  int N = m_Nodes.size();
  EG_VTKSP(vtkPoints, points);
  points->SetNumberOfPoints(N);
  for (int i = 0; i < N; ++i) {
    vec3_t x = getXNode(i);
    points->SetPoint(i, x.data());
  }
  poly_data->Allocate(m_Faces.size());
  poly_data->SetPoints(points);
  for (int i = 0; i < m_Faces.size(); ++i) {
    QVector<vtkIdType> pts(m_Faces[i].size());
    for (int j = 0; j < m_Faces[i].size(); ++j) {
      pts[j] = m_Faces[i][j];
    }
    if (pts.size() == 3) {
      poly_data->InsertNextCell(VTK_TRIANGLE, pts.size(), pts.data());
    } else if (pts.size() == 4) {
      poly_data->InsertNextCell(VTK_QUAD, pts.size(), pts.data());
    } else {
      poly_data->InsertNextCell(VTK_POLYGON, pts.size(), pts.data());
    }
  }
}

void PolyMolecule::computeCentreOfGravity()
{
  m_AllPositive = true;
  // first guess of centre of gravity
  vec3_t x_centre(0,0,0);
  for (int i = 0; i < m_Faces.size(); ++i) {
    x_centre += getXFace(i);
  }
  m_N2BadFaces.resize(m_Nodes.size());
  for (int i = 0; i < m_Nodes.size(); ++i) {
    m_N2BadFaces[i].clear();
  }
  x_centre *= 1.0/m_Faces.size();
  for (int iter = 0; iter < 2; ++iter) {
    double V_total = 0;
    m_CentreOfGravity = vec3_t(0,0,0);
    m_MaxPyramidVolume = -1e99;
    m_MinPyramidVolume =  1e99;
    for (int i = 0; i < m_Faces.size(); ++i) {
      double V_pyramid = 0;
      vec3_t x_face = getXFace(i);
      QVector<vec3_t> x(m_Faces[i].size() + 1);
      for (int j = 0; j < m_Faces[i].size(); ++j) {
        x[j] = getXNode(m_Faces[i][j]);
      }
      x.last() = x.first();      
      for (int j = 0; j < m_Faces[i].size(); ++j) {
        double V = GeometryTools::tetraVol(x_face, x[j+1], x[j], x_centre, true);
        if (V <= 0) {
          m_AllPositive = false;
        }
        V_total += V;
        V_pyramid += V;
        m_CentreOfGravity += 0.25*V*(x_face + x[j+1] + x[j] + x_centre);
      }
      if (V_pyramid <= 0) {
        for (int j = 0; j < m_Faces[i].size(); ++j) {
          m_N2BadFaces[m_Faces[i][j]].insert(i);
        }
      }
      m_MaxPyramidVolume = max(V_pyramid, m_MaxPyramidVolume);
      m_MinPyramidVolume = min(V_pyramid, m_MinPyramidVolume);
    }
    m_CentreOfGravity *= 1.0/V_total;
    x_centre = m_CentreOfGravity;
  }
}

void PolyMolecule::buildNode2Node()
{
  m_N2N.clear();
  m_N2N.resize(m_Nodes.size());
  for (int i = 0; i < m_Faces.size(); ++i) {
    for (int j = 0; j < m_Faces[i].size() - 1; ++j) {
      m_N2N[m_Faces[i][j]].insert(m_Faces[i][j+1]);
      m_N2N[m_Faces[i][j+1]].insert(m_Faces[i][j]);
    }
  }
}

void PolyMolecule::computeNormals()
{
  m_FaceNormals.fill(vec3_t(0,0,0), m_Faces.size());
  m_NodeNormals.fill(vec3_t(0,0,0), m_Nodes.size());
  for (int i = 0; i < m_Faces.size(); ++i) {
    QVector<vec3_t> x(m_Faces[i].size() + 1);
    vec3_t x_face(0,0,0);
    for (int j = 0; j < m_Faces[i].size(); ++j) {
      x[j] = getXNode(m_Faces[i][j]);
      x_face += x[j];
    }
    x.last() = x.first();
    x_face *= 1.0/m_Faces[i].size();
    m_FaceNormals[i] = vec3_t(0,0,0);
    for (int j = 0; j < m_Faces[i].size(); ++j) {
      m_FaceNormals[i] += triNormal(x_face, x[j], x[j+1]);
    }
    for (int j = 0; j < m_Faces[i].size(); ++j) {
      vec3_t n = m_FaceNormals[i];
      //n.normalise();
      m_NodeNormals[m_Faces[i][j]] += n;
    }
  }
  for (int i = 0; i < m_Nodes.size(); ++i) {
    m_NodeNormals[i].normalise();
  }
}

vec3_t PolyMolecule::getXFace(int face)
{
  vec3_t x(0,0,0);
  foreach (int i, m_Faces[face]) {
    x += getXNode(i);
  }
  if (m_Nodes.size() == 0) {
    EG_BUG;
  }
  return (1.0/m_Faces[face].size())*x;
}

void PolyMolecule::smooth(bool delaunay, bool write)
{
  computeNormals();
  if (write) writeVtkFile("bad_cell");

  // collect all "healthy" cells in the neighbourhood
  // they will be checked for volume inversion caused by the smoothing of this cell
  //
  QSet<int> check_cells;
  check_cells.insert(m_PCell);
  for (int i = 0; i < m_PolyMesh->numFacesOfPCell(m_PCell); ++i) {
    int i_face = m_PolyMesh->pcell2Face(m_PCell, i);
    int i_cell = m_PolyMesh->owner(i_face);
    PolyMolecule pm(m_PolyMesh, i_cell);
    if (pm.minPyramidVolume() > 0) {
      check_cells.insert(i_cell);
    }
    i_cell = m_PolyMesh->neighbour(i_face);
    if (i_cell != -1) {
      PolyMolecule pm(m_PolyMesh, i_cell);
      if (pm.minPyramidVolume() > 0) {
        check_cells.insert(i_cell);
      }
    }
  }

  // VTK pipeline to get convex hull
  // The cell will be "inflated" into this convex hull.
  //
  QVector<QVector<vec3_t> > convex_hull;
  {
    EG_VTKSP(vtkPolyData, poly_data);
    EG_VTKSP(vtkPoints, points);
    points->SetNumberOfPoints(m_Nodes.size());
    for (vtkIdType id_node = 0; id_node < m_Nodes.size(); ++id_node) {
      vec3_t x = getXNode(id_node);
      points->SetPoint(id_node, x.data());
    }
    poly_data->Allocate(m_Faces.size());
    poly_data->SetPoints(points);
    for (int i = 0; i < m_Faces.size(); ++i) {
      QVector<vtkIdType> pts(m_Faces[i].size());
      for (int j = 0; j < m_Faces[i].size(); ++j) {
        pts[j] = m_Faces[i][j];
      }
      poly_data->InsertNextCell(VTK_POLYGON, pts.size(), pts.data());
    }
    EG_VTKSP(vtkGeometryFilter, surface);
    if (delaunay) {
      EG_VTKSP(vtkDelaunay3D, delaunay);
      delaunay->SetOffset(100);
      delaunay->SetInput(poly_data);
      surface->SetInput(delaunay->GetOutput());
    } else {
      EG_VTKSP(vtkHull, hull);
      hull->SetInput(poly_data);
      hull->AddRecursiveSpherePlanes(5);
      EG_VTKSP(vtkTriangleFilter, tri);
      tri->SetInput(hull->GetOutput());
      surface->SetInput(tri->GetOutput());
    }
    surface->Update();
    if (surface->GetOutput()->GetNumberOfCells() < 4) {
      return;
    }
    EG_VTKSP(vtkXMLPolyDataWriter, vtp);
    if (write) {
      QString file_name = GuiMainWindow::pointer()->getCwd() + "/" + "convex_hull.vtp";
      vtp->SetFileName(qPrintable(file_name));
      vtp->SetDataModeToAscii();
      vtp->SetInput(surface->GetOutput());
      vtp->Write();
    }
    convex_hull.fill(QVector<vec3_t>(3), surface->GetOutput()->GetNumberOfCells());
    for (vtkIdType id_tri = 0; id_tri < surface->GetOutput()->GetNumberOfCells(); ++id_tri) {
      EG_GET_CELL(id_tri, surface->GetOutput());
      if (type_cell != VTK_TRIANGLE) {
        EG_BUG;
      }
      for (int i = 0; i < 3; ++i) {
        vec3_t x;
        surface->GetOutput()->GetPoint(pts[i], x.data());
        convex_hull[id_tri][i] = x;
      }
    }
  }

  // Save old nodes for under-relaxations
  //
  QVector<vec3_t> x_old(m_Nodes.size());
  for (int i = 0; i < m_Nodes.size(); ++i) {
    x_old[i] = getXNode(i);
  }

  // "inflate" cell
  //
  for (int iter = 1; iter <= 100; ++iter) {
    for (int i = 0; i < m_Nodes.size(); ++i) {
      vec3_t x = getXNode(i);
      double L = (x - m_CentreOfGravity).abs();

      vec3_t Dx_n = 0.1*L*m_NodeNormals[i];
      vec3_t Dx_lp(0,0,0);
      foreach (int neigh, m_N2N[i]) {
        Dx_lp += 0.1*(getXNode(neigh) - x);
      }
      //Dx_lp -= (Dx_lp * m_NodeNormals[i])*m_NodeNormals[i];
      vec3_t Dx = Dx_n + Dx_lp;

      for (int j = 0; j < convex_hull.size(); ++j) {
        vec3_t np = triNormal(convex_hull[j][0], convex_hull[j][1], convex_hull[j][2]);
        np.normalise();
        double d = GeometryTools::intersection(x, Dx, convex_hull[j][0], np);
        if (d > -1e-3 && d < 1) {
          Dx -= (Dx*np)*np;
        }
      }

      for (int j = 0; j < convex_hull.size(); ++j) {
        vec3_t np = triNormal(convex_hull[j][0], convex_hull[j][1], convex_hull[j][2]);
        np.normalise();
        double d = GeometryTools::intersection(x, Dx, convex_hull[j][0], np);
        if (d > -1e-3 && d < 1) {
          Dx *= 0;
        }
      }

      x += Dx;
      m_PolyMesh->setNodeVector(m_Nodes[i], x);
    }
    QString file_name;
    file_name.setNum(iter);
    file_name = file_name.rightJustified(4, '0');
    file_name = "improved_cell." + file_name;
    if (write) writeVtkFile(file_name);
    computeCentreOfGravity();
    computeNormals();
  }

  // get new state of the cell for under relaxation
  //
  QVector<vec3_t> x_new(m_Nodes.size());
  for (int i = 0; i < m_Nodes.size(); ++i) {
    x_new[i] = getXNode(i);
  }

  // check all formerly "healthy" cells and under relax as rquired
  //
  bool done = false;
  double relax = 1;
  while (!done) {
    done = true;
    foreach (int i_pcell, check_cells) {
      PolyMolecule pm(m_PolyMesh, i_pcell);
      if (!pm.allPositive()) {
        done = false;
        break;
      }
    }
    if (!done) {
      relax -= 0.01;
      if (relax < 1e-3) {
        //relax = 1.0;
        done = true;
      }
      for (int i = 0; i < m_Nodes.size(); ++i) {
        m_PolyMesh->setNodeVector(m_Nodes[i], x_old[i] + relax*(x_new[i] - x_old[i]));
      }
    }
  }


  if (write) writeVtkFile("improved_cell");
}

void PolyMolecule::split(bool write)
{
  if (write) writeVtkFile("concave_cell");
}

void PolyMolecule::fix(bool write)
{
  /*
  if (m_MinPyramidVolume/m_MaxPyramidVolume < -0.1) {
    split(true);
    EG_BUG;
  }
  */

  smooth(true, write);
  if (m_MinPyramidVolume <= 0) {
    smooth(false, write);
  }
  //EG_BUG;
}



