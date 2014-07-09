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
    m_FacesInPMesh.append(i_face);
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
  vtu->SetInputData(grid);
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

void PolyMolecule::buildFace2Face()
{
  m_F2F.clear();
  m_F2F.resize(m_Faces.size());
  for (int face = 0; face < m_Faces.size(); ++face) {
    for (int i = 0; i < m_Faces[face].size(); ++i) {
      int j = i + 1;
      if (j >= m_Faces[face].size()) {
        j = 0;
      }
      edge_t e = getEdge(m_Faces[face][i], m_Faces[face][j]);
      if (e.face1 != face) {
        m_F2F[face].insert(e.face1);
      }
      if (e.face2 != face) {
        m_F2F[face].insert(e.face2);
      }
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
  if (write) {
    writeVtkFile("bad_cell");
  }

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
      delaunay->SetInputData(poly_data);
      surface->SetInputConnection(delaunay->GetOutputPort());
    } else {
      EG_VTKSP(vtkHull, hull);
      hull->SetInputData(poly_data);
      hull->AddRecursiveSpherePlanes(5);
      EG_VTKSP(vtkTriangleFilter, tri);
      tri->SetInputConnection(hull->GetOutputPort());
      surface->SetInputConnection(tri->GetOutputPort());
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
      vtp->SetInputConnection(surface->GetOutputPort());
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

PolyMolecule::edge_t PolyMolecule::getEdge(int node1, int node2)
{
  edge_t e;
  e.node1 = node1;
  e.node2 = node2;
  e.face1 = -1;
  e.face2 = -1;
  for (int i = 0; i < m_Faces.size(); ++i) {
    for (int j = 0; j < m_Faces[i].size(); ++j) {
      int n1 = m_Faces[i][j];
      int n2 = m_Faces[i][0];
      if (j < m_Faces[i].size() - 1) {
        n2 = m_Faces[i][j+1];
      }
      if (n1 == node1 && n2 == node2) {
        e.face1 = i;
      }
      if (n1 == node2 && n2 == node1) {
        e.face2 = i;
      }
    }
  }
  if (e.face1 == -1 || e.face2 == -1) {
    EG_BUG;
  }
  return e;
}

QList<PolyMolecule::edge_t> PolyMolecule::findConcaveEdges()
{
  if (m_PCell == 32) {
    cout << "break" << endl;
  }
  QList<edge_t> concave_edges;
  for (int node1 = 0; node1 < m_Nodes.size(); ++node1) {
    foreach (int node2, m_N2N[node1]) {
      if (node2 > node1) {
        edge_t e = getEdge(node1, node2);
        vec3_t x1 = getXNode(node1);
        vec3_t x2 = getXNode(node2);
        vec3_t xe  = 0.5*(getXNode(node1) + getXNode(node2));
        vec3_t xf1 = getXFace(e.face1);
        vec3_t xf2 = getXFace(e.face2);
        vec3_t xc  = m_CentreOfGravity;
        vec3_t xn  = 0.5*(xf1 + xf2);
        double Lnorm = (xf2 - xf1).abs();
        double L = ((xn - xc).abs() - (xe - xc).abs())/Lnorm;
        if (L > 0.01) {
          concave_edges.append(e);
        }
      }
    }
  }
  return concave_edges;
}

void PolyMolecule::split(bool write)
{
  if (write) {
    writeVtkFile("concave_cell");
  }
  buildNode2Node();
  buildFace2Face();
  QList<edge_t> concave_edges = findConcaveEdges();
  if (concave_edges.size() == 0) {
    return;
  }
  QList<int> seed_faces;
  foreach (edge_t e, concave_edges) {
    if (!seed_faces.contains(e.face1)) {
      seed_faces.append(e.face1);
    }
    if (!seed_faces.contains(e.face2)) {
      seed_faces.append(e.face2);
    }
  }
  if (seed_faces.size() != concave_edges.size() + 1) {
    //EG_BUG;
  }
  QVector<int> old2new(m_Faces.size());
  int num_faces = m_Faces.size();
  for (int face = 0; face < m_Faces.size(); ++face) {
    old2new[face] = seed_faces.indexOf(face);
    if (old2new[face] != -1) {
      --num_faces;
    }
  }
  while (num_faces > 0) {
    for (int face1 = 0; face1 < m_Faces.size(); ++face1) {
      if (old2new[face1] != -1) {
        foreach (int face2, m_F2F[face1]) {
          if (old2new[face2] == -1) {
            old2new[face2] = old2new[face1];
            --num_faces;
          }
        }
      }
    }
  }
  QVector<int> new_pcell_index(seed_faces.size());
  new_pcell_index[0] = m_PCell;
  for (int i = 1; i < seed_faces.size(); ++i) {
    new_pcell_index[i] = m_PolyMesh->numCells() + i - 1;
  }
  QVector<PolyMesh::face_t> new_faces(concave_edges.size());
  QVector<vec3_t> new_cell_centre(seed_faces.size(), vec3_t(0,0,0));
  {
    QVector<int> count(seed_faces.size(), 0);
    for (int face = 0; face < m_Faces.size(); ++face) {
      new_cell_centre[old2new[face]] += getXFace(face);
      ++count[old2new[face]];
    }
    for (int i = 0; i < new_cell_centre.size(); ++i) {
      if (count[i] == 0) {
        EG_BUG;
      }
      new_cell_centre[i] *= 1.0/count[i];
    }
  }
  for (int i = 0; i < concave_edges.size(); ++i) {
    edge_t e = concave_edges[i];
    int start = e.node1;
    int first = e.node2;
    new_faces[i].node.append(m_Nodes[e.node1]);
    int o2n1 = old2new[e.face1];
    int o2n2 = old2new[e.face2];
    new_faces[i].ref_vec = new_cell_centre[o2n2] - new_cell_centre[o2n1];
    new_faces[i].owner = new_pcell_index[o2n1];
    new_faces[i].neighbour = new_pcell_index[o2n2];
    new_faces[i].bc = 0;
    bool reversed = false;
    while (e.node2 != start) {
      if (reversed) {
        new_faces[i].node.prepend(m_Nodes[e.node2]);
      } else {
        new_faces[i].node.append(m_Nodes[e.node2]);
      }
      bool found = false;
      foreach (int n, m_N2N[e.node2]) {
        if (n != e.node1) {
          edge_t ne = getEdge(e.node2, n);
          if (old2new[ne.face1] == o2n1 && old2new[ne.face2] == o2n2) {
            e = ne;
            found = true;
            break;
          }
        }
      }

      // reached a dead-end here ...
      if (!found && !reversed) {
        reversed = true;
        swap(o2n1, o2n2);
        swap(start, first);
        foreach (int n, m_N2N[first]) {
          if (n != start) {
            edge_t ne = getEdge(first, n);
            if (old2new[ne.face1] == o2n1 && old2new[ne.face2] == o2n2) {
              e = ne;
              found = true;
              break;
            }
          }
        }
      }
      if (!found) {
        break;
      }
    }
    if (new_faces[i].owner > new_faces[i].neighbour && new_faces[i].neighbour != -1) {
      swap(new_faces[i].owner, new_faces[i].neighbour);
      new_faces[i].ref_vec *= -1;
      invertQContainer(new_faces[i].node);
    }
  }
  for (int face = 0; face < m_Faces.size(); ++face) {
    int pface = m_FacesInPMesh[face];
    //m_PolyMesh->m_Faces[pface].bc = 1;
    if (m_PolyMesh->m_Faces[pface].owner == m_PCell) {
      m_PolyMesh->m_Faces[pface].owner = new_pcell_index[old2new[face]];
    }
    if (m_PolyMesh->m_Faces[pface].neighbour == m_PCell) {
      m_PolyMesh->m_Faces[pface].neighbour = new_pcell_index[old2new[face]];
    }
    if (m_PolyMesh->m_Faces[pface].owner > m_PolyMesh->m_Faces[pface].neighbour && m_PolyMesh->m_Faces[pface].neighbour != -1) {
      swap(m_PolyMesh->m_Faces[pface].owner, m_PolyMesh->m_Faces[pface].neighbour);
      m_PolyMesh->m_Faces[pface].ref_vec *= -1;
      invertQContainer(m_PolyMesh->m_Faces[pface].node);
    }
  }
  foreach (PolyMesh::face_t f, new_faces) {
    if (f.owner == -1) {
      EG_BUG;
    }
    if (f.neighbour == -1) {
      EG_BUG;
    }
    if (f.owner > f.neighbour) {
      swap(f.owner, f.neighbour);
      invertQContainer(f.node);
    }
    //f.bc = 2;
    m_PolyMesh->m_Faces.append(f);
  }
  m_PolyMesh->m_NumPolyCells += seed_faces.size() - 1;
}

void PolyMolecule::updatePMesh()
{
  int offset = 0;
  foreach (int face, m_FacesInPMesh) {
    m_PolyMesh->m_Faces.removeAt(face - offset);
    ++offset;
  }
}

void PolyMolecule::optimise(bool write)
{
  smooth(true, write);
  if (m_MinPyramidVolume <= 0) {
    smooth(false, write);
  }
}

void PolyMolecule::updateFace(int face, int new_cell_index)
{
  PolyMesh::face_t f = m_PolyMesh->m_Faces[m_FacesInPMesh[face]];
  if (f.owner == m_PCell) {
    f.owner = new_cell_index;
  }
  if (f.neighbour == m_PCell) {
    f.neighbour = new_cell_index;
  }
  if (f.owner > f.neighbour && f.neighbour != -1) {
    swap(f.owner, f.neighbour);
    invertQContainer(f.node);
    f.ref_vec *= -1;
  }
  m_PolyMesh->m_Faces[m_FacesInPMesh[face]] = f;
}

void PolyMolecule::centreSplit()
{
  buildNode2Node();
  buildFace2Face();
  QList<edge_t> concave_edges = findConcaveEdges();
  if (concave_edges.size() == 0) {
    return;
  }

  // find common nodes in all edges and return if none exists
  QSet<int> common_nodes;
  common_nodes.insert(concave_edges.first().node1);
  common_nodes.insert(concave_edges.first().node2);
  foreach (edge_t e, concave_edges) {
    QSet<int> edge_nodes;
    edge_nodes.insert(e.node1);
    edge_nodes.insert(e.node2);
    common_nodes.intersect(edge_nodes);
  }
  if (common_nodes.size() == 0) {
    return;
  }

  // start point of straight for divide point computation
  vec3_t x1;
  {
    int count = 0;
    foreach(int node, common_nodes) {
      x1 += getXNode(node);
      ++count;
    }
    x1 *= 1.0/count;
  }

  // end point of straight for divide point computation
  computeNormals();
  vec3_t x2;
  {
    vec3_t n(0,0,0);
    QList<int> adjacent_faces;
    foreach (edge_t e, concave_edges) {
      if (!adjacent_faces.contains(e.face1)) {
        adjacent_faces.append(e.face1);
      }
      if (!adjacent_faces.contains(e.face2)) {
        adjacent_faces.append(e.face2);
      }
    }
    foreach (int face, adjacent_faces) {
      n += m_FaceNormals[face];
    }
    n.normalise();
    bool x2_set = false;
    double max_dist = 0;
    for (int node = 0; node < m_Nodes.size(); ++node) {
      vec3_t x = getXNode(node);
      double dist = (x1 - x)*n;
      if (dist > max_dist) {
        x2 = x1 - dist*n;
        max_dist = dist;
        x2_set = true;
      }
    }
    if (!x2_set) {
      EG_BUG;
    }
  }

  // tip point for all polygonial pyramids
  vec3_t x_tip = 0.5*(x1 + x2);

  // vector with new cell indices
  QVector<int> cell_index(m_Faces.size());
  cell_index[0] = m_PCell;
  for (int i = 1; i < cell_index.size(); ++i) {
    cell_index[i] = m_PolyMesh->m_NumPolyCells++;
  }

  // new node index
  int new_node = m_PolyMesh->m_Points.size();
  m_PolyMesh->m_Points.append(x_tip);

  // re-map existing faces
  for (int face = 0; face < m_FacesInPMesh.size(); ++face) {
    updateFace(face, cell_index[face]);
  }

  // create new faces
  for (int face = 0; face < m_Faces.size(); ++face) {
    QList<int> nodes = m_Faces[face];
    nodes.append(nodes.first());
    for (int i_node = 0; i_node < nodes.size() - 1; ++i_node) {
      int node1 = nodes[i_node];
      int node2 = nodes[i_node + 1];
      edge_t e = getEdge(node1, node2);
      PolyMesh::face_t f;
      if (e.face1 == face) {
        f.neighbour = cell_index[e.face2];
      } else {
        f.neighbour = cell_index[e.face1];
      }
      if (f.owner < f.neighbour) {
        f.bc = 0;
        f.node.append(m_Nodes[node2]);
        f.node.append(m_Nodes[node1]);
        f.node.append(new_node);
        f.owner = cell_index[face];
        vec3_t x1 = m_PolyMesh->m_Points[f.node[0]];
        vec3_t x2 = m_PolyMesh->m_Points[f.node[1]];
        vec3_t x3 = m_PolyMesh->m_Points[f.node[2]];
        f.ref_vec = GeometryTools::triNormal(x1, x2, x3);
        /*
        if (f.owner > f.neighbour) {
          swap(f.owner, f.neighbour);
          f.ref_vec *= -1;
          invertQContainer(f.node);
        }
        */
        m_PolyMesh->m_Faces.append(f);
      }
    }
  }
}


