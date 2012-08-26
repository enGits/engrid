#include "orthogonalityoptimiser.h"
#include "surfaceprojection.h"
#include "guimainwindow.h"

OrthogonalityOptimiser::OrthogonalityOptimiser()
{
  m_StrictPrismChecking = true;
}

double OrthogonalityOptimiser::faceError(const face_t &face, vec3_t x)
{
  vec3_t n = face.c1 - face.c2;
  n.normalise();
  vec3_t xm = x;
  foreach (vtkIdType id_node, face.nodes) {
    if (id_node != m_IdNode) {
      vec3_t xn;
      m_Grid->GetPoint(id_node, xn.data());
      xm += xn;
    }
  }
  xm *= 1.0/face.nodes.size();
  vec3_t v = x - xm;
  if (!checkVector(v)) {
    EG_BUG;
  }
  return sqr(v*n);
}

double OrthogonalityOptimiser::func(vec3_t x)
{
  if (m_IdNode == 427) {
    //cout << "stop" << endl;
  }
  double err = 0;
  foreach (int i_face, m_Node2Face[m_IdNode]) {
    double e = faceError(m_Faces[i_face], x);
    if (isnan(e)) {
      EG_BUG;
    }
    if (isinf(e)) {
      EG_BUG;
    }
    err += e;
  }
  return err;
}

void OrthogonalityOptimiser::getFace(vtkIdType id_cell, int subidx, QSet<vtkIdType> &nodes)
{
  nodes.clear();
  vtkIdType N_pts, *pts;
  m_Grid->GetCellPoints(id_cell, N_pts, pts);
  vtkIdType type_cell = m_Grid->GetCellType(id_cell);
  if (type_cell == VTK_TETRA) {
    if (subidx == 0) {
      nodes.insert(pts[2]);
      nodes.insert(pts[1]);
      nodes.insert(pts[0]);
    }
    if (subidx == 1) {
      nodes.insert(pts[1]);
      nodes.insert(pts[3]);
      nodes.insert(pts[0]);
    }
    if (subidx == 2) {
      nodes.insert(pts[3]);
      nodes.insert(pts[2]);
      nodes.insert(pts[0]);
    }
    if (subidx == 3) {
      nodes.insert(pts[2]);
      nodes.insert(pts[3]);
      nodes.insert(pts[1]);
    }
  }
  if (type_cell == VTK_WEDGE) {
    if (subidx == 0) {
      nodes.insert(pts[0]);
      nodes.insert(pts[1]);
      nodes.insert(pts[2]);
    }
    if (subidx == 1) {
      nodes.insert(pts[3]);
      nodes.insert(pts[5]);
      nodes.insert(pts[4]);
    }
    if (subidx == 2) {
      nodes.insert(pts[3]);
      nodes.insert(pts[4]);
      nodes.insert(pts[1]);
      nodes.insert(pts[0]);
    }
    if (subidx == 3) {
      nodes.insert(pts[1]);
      nodes.insert(pts[4]);
      nodes.insert(pts[2]);
      nodes.insert(pts[5]);
    }
    if (subidx == 4) {
      nodes.insert(pts[0]);
      nodes.insert(pts[2]);
      nodes.insert(pts[5]);
      nodes.insert(pts[3]);
    }
  }
  if (type_cell == VTK_HEXAHEDRON) {
    if (subidx == 0) {
      nodes.insert(pts[0]);
      nodes.insert(pts[3]);
      nodes.insert(pts[2]);
      nodes.insert(pts[1]);
    }
    if (subidx == 1) {
      nodes.insert(pts[4]);
      nodes.insert(pts[5]);
      nodes.insert(pts[6]);
      nodes.insert(pts[7]);
    }
    if (subidx == 2) {
      nodes.insert(pts[0]);
      nodes.insert(pts[1]);
      nodes.insert(pts[5]);
      nodes.insert(pts[4]);
    }
    if (subidx == 3) {
      nodes.insert(pts[3]);
      nodes.insert(pts[7]);
      nodes.insert(pts[6]);
      nodes.insert(pts[2]);
    }
    if (subidx == 4) {
      nodes.insert(pts[0]);
      nodes.insert(pts[4]);
      nodes.insert(pts[7]);
      nodes.insert(pts[3]);
    }
    if (subidx == 5) {
      nodes.insert(pts[1]);
      nodes.insert(pts[2]);
      nodes.insert(pts[6]);
      nodes.insert(pts[5]);
    }
  }
}

int OrthogonalityOptimiser::getNumberOfFaces(vtkIdType id_cell)
{
  vtkIdType cell_type = m_Grid->GetCellType(id_cell);
  int N = 0;
  if      (cell_type == VTK_TETRA)      N = 4;
  else if (cell_type == VTK_PYRAMID)    N = 5;
  else if (cell_type == VTK_WEDGE)      N = 5;
  else if (cell_type == VTK_HEXAHEDRON) N = 6;
  return N;
}

void OrthogonalityOptimiser::buildNode2Face()
{
  l2g_t cells = m_Part.getCells();
  m_Node2Face.resize(m_Grid->GetNumberOfPoints());
  QList<face_t> faces;
  foreach (vtkIdType id_cell1, cells) {
    if (!isSurface(id_cell1, m_Grid)) {
      for (int i = 0; i < m_Part.c2cGSize(id_cell1); ++i) {
        vtkIdType id_cell2 = m_Part.c2cGG(id_cell1, i);
        if (id_cell1 < id_cell2 && !isSurface(id_cell2, m_Grid)) {
          for (int j = 0; j < getNumberOfFaces(id_cell1); ++j) {
            QSet<vtkIdType> nodes1;
            getFace(id_cell1, j, nodes1);
            for (int k = 0; k < getNumberOfFaces(id_cell2); ++k) {
              QSet<vtkIdType> nodes2;
              getFace(id_cell2, j, nodes2);
              if (nodes1 == nodes2) {
                face_t F;
                F.c1 = cellCentre(m_Grid, id_cell1);
                F.c2 = cellCentre(m_Grid, id_cell2);
                F.nodes.resize(nodes1.size());
                qCopy(nodes1.begin(), nodes1.end(), F.nodes.begin());
                faces.append(F);
              }
            }
          }
        }
      }
    }
  }
  m_Faces.resize(faces.size());
  qCopy(faces.begin(), faces.end(), m_Faces.begin());
  for (int i = 0; i < m_Faces.size(); ++i) {
    foreach (vtkIdType id_node, m_Faces[i].nodes) {
      m_Node2Face[id_node].append(i);
    }
  }
}

void OrthogonalityOptimiser::fixBoundaryNodes()
{
  l2g_t cells = m_Part.getCells();
  m_Fixed.fill(false, m_Grid->GetNumberOfPoints());
  EG_VTKDCC(vtkIntArray, cell_code, m_Grid, "cell_code");
  QVector<QSet<int> > n2bc(m_Grid->GetNumberOfPoints());
  foreach (vtkIdType id_cell, cells) {
    if (isSurface(id_cell, m_Grid)) {
      vtkIdType N_pts, *pts;
      m_Grid->GetCellPoints(id_cell, N_pts, pts);
      for (int i = 0; i < N_pts; ++i) {
        n2bc[pts[i]].insert(cell_code->GetValue(id_cell));
      }
    }
  }
  m_Projection.fill(NULL, m_Grid->GetNumberOfPoints());
  for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
    if (n2bc[id_node].size() > 1) {
      m_Fixed[id_node] = true;
    } else if (n2bc[id_node].size() == 1) {
      GuiMainWindow::pointer()->getSurfProj(n2bc[id_node].toList().first());
    }
    if (m_Node2Face[id_node].size() == 0) {
      m_Fixed[id_node] = true;
    }
  }
}

vec3_t OrthogonalityOptimiser::newPosition(vtkIdType id_node)
{
  vec3_t x_old;
  m_Grid->GetPoint(id_node, x_old.data());
  if (m_Node2Face[id_node].size() == 0) {
    return x_old;
  }
  vec3_t x_new(0,0,0);
  foreach (int i_face, m_Node2Face[id_node]) {
    const face_t &face = m_Faces[i_face];
    vec3_t n = face.c1 - face.c2;
    n.normalise();
    vec3_t xm(0,0,0);
    foreach (vtkIdType id_node, face.nodes) {
      vec3_t xn;
      m_Grid->GetPoint(id_node, xn.data());
      xm += xn;
    }
    xm *= 1.0/face.nodes.size();
    vec3_t v = x_old - xm;
    if (!checkVector(v)) {
      EG_BUG;
    }
    vec3_t Dx = (v*n)*n;
    x_new += x_old - Dx;
  }
  x_new *= (1.0/m_Node2Face[id_node].size());
  return x_new;
}

void OrthogonalityOptimiser::operate()
{
  computeNormals();
  buildNode2Face();
  fixBoundaryNodes();
  l2g_t nodes = m_Part.getNodes();
  int N = 0;
  foreach (vtkIdType id_node, nodes) {
    if (!m_Fixed[id_node]) {
      m_IdNode = id_node;
      vec3_t x_old;
      m_Grid->GetPoint(id_node, x_old.data());
      //vec3_t x_new = optimise(x_old);
      vec3_t x_new = newPosition(id_node);
      if (checkVector(x_new)) {
        if (m_Projection[id_node]) {
          x_new = m_Projection[id_node]->project(x_new, id_node);
        }
        if (checkVector(x_new)) {
          vec3_t Dx = x_new - x_old;
          if (moveNode(id_node, Dx)) ++N;
        }
      }
    }
  }
  cout << N << " nodes optimised" << endl;
}
