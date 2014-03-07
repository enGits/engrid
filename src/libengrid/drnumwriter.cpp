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

#include "drnumwriter.h"
#include "guimainwindow.h"

DrNumWriter::DrNumWriter()
{
  m_OverlapLayers = 2;
}

QString DrNumWriter::boundaryCode(vtkIdType id_cell, int i)
{
  EG_VTKDCC(vtkIntArray, cell_code, m_Grid, "cell_code");
  QString bcode = "0";
  vtkIdType id_neigh = m_Part.c2cGG(id_cell, i);
  if (id_neigh != -1) {
    if (m_CellToCartPatch[id_neigh] == -1) {
      if (isSurface(id_neigh, m_Grid)) {
        BoundaryCondition bc = GuiMainWindow::pointer()->getBC(cell_code->GetValue(id_neigh));
        if (bc.getName() != "unknown") {
          bcode = bc.getType();
        }
      }
    }
  }
  return bcode;
}

void DrNumWriter::computeMeshDensity()
{
  QString buffer = GuiMainWindow::pointer()->getXmlSection("engrid/surface/settings").replace("\n", " ");
  if (!buffer.isEmpty()) {
    QTextStream in(&buffer, QIODevice::ReadOnly);
    in >> m_MaxEdgeLength;
    in >> m_MinEdgeLength;
    in >> m_GrowthFactor;
  } else {
    m_MaxEdgeLength = 1000.0;
    m_MinEdgeLength = 0.0;
    m_GrowthFactor = 1.5;
  }
  m_ELSManager.read();
  m_H.fill(m_MaxEdgeLength, m_Grid->GetNumberOfCells());

  for (vtkIdType id_cell = 0; id_cell < m_Grid->GetNumberOfCells(); ++id_cell) {
    vec3_t x = cellCentre(m_Grid, id_cell);
    double cl_src = m_ELSManager.minEdgeLength(x);
    if (cl_src > 0) {
      if (cl_src < m_H[id_cell]) {
        m_H[id_cell] = cl_src;
      }
    }
  }

  double H_min = 1e99;
  vtkIdType id_min = -1;
  for (vtkIdType id_cell = 0; id_cell < m_Grid->GetNumberOfCells(); ++id_cell) {
    if (m_H[id_cell] < H_min && isVolume(id_cell, m_Grid)) {
      id_min = id_cell;
      H_min = m_H[id_cell];
    }
  }
  if (id_min < 0) {
    EG_BUG;
  }

  QVector<bool> marked(m_Grid->GetNumberOfCells(), false);
  marked[id_min] = true;
  bool done = false;
  while (!done) {
    done = true;
    for (vtkIdType id_cell = 0; id_cell < m_Grid->GetNumberOfCells(); ++id_cell) {
      if (marked[id_cell] && m_H[id_cell] <= H_min && isVolume(id_cell, m_Grid)) {
        vtkIdType num_pts, *pts;
        m_Grid->GetCellPoints(id_cell, num_pts, pts);
        for (int i = 0; i < num_pts; ++i) {
          for (int j = 0; j < m_Part.n2cGSize(pts[i]); ++j) {
            vtkIdType id_neigh = m_Part.n2cGG(pts[i], j);
            if (!marked[id_neigh] && isVolume(id_neigh, m_Grid)) {
              double h = m_GrowthFactor*m_H[id_cell];
              if (h < 0) {
                EG_BUG;
              }
              m_H[id_neigh] = min(m_H[id_neigh], h);
              marked[id_neigh] = true;
              done = false;
            }
          }
        }
      }
    }
    H_min *= m_GrowthFactor;
  }

  for (vtkIdType id_cell = 0; id_cell < m_Grid->GetNumberOfCells(); ++id_cell) {
    vec3_t x = cellCentre(m_Grid, id_cell);
    double cl_src = m_ELSManager.minEdgeLength(x);
    if (cl_src > 0) {
      if (cl_src < m_H[id_cell]) {
        m_H[id_cell] = cl_src;
      }
    }
  }

}

void DrNumWriter::createAllCartPatches()
{
  m_CellToCartPatch.fill(-1, m_Grid->GetNumberOfCells());
  int num_cart_patches = 0;
  for (vtkIdType id_cell = 0; id_cell < m_Grid->GetNumberOfCells(); ++id_cell) {
    vtkIdType type_cell = m_Grid->GetCellType(id_cell);
    if (type_cell == VTK_HEXAHEDRON) {
      m_CellToCartPatch[id_cell] = num_cart_patches;
      ++num_cart_patches;
    }
  }
  m_CartPatches.resize(num_cart_patches);
  EG_VTKDCN(vtkDoubleArray, cl, m_Grid, "node_meshdensity_desired");
  for (vtkIdType id_cell = 0; id_cell < m_Grid->GetNumberOfCells(); ++id_cell) {
    if (m_CellToCartPatch[id_cell] != -1) {
      cart_patch_t patch = m_CartPatches[m_CellToCartPatch[id_cell]];
      patch.id_cell = id_cell;
      vtkIdType num_pts, *pts;
      vec3_t x[8];
      m_Grid->GetCellPoints(id_cell, num_pts, pts);
      if (num_pts != 8) {
        EG_BUG;
      }
      vec3_t xc = cellCentre(m_Grid, id_cell);
      for (int i = 0; i < 8; ++i) {
        m_Grid->GetPoint(pts[i], x[i].data());
      }
      double h = m_H[id_cell];
      vec3_t gi = 0.25*(x[1] + x[2] + x[6] + x[5] - x[0] - x[3] - x[7] - x[4]);
      vec3_t gj = 0.25*(x[4] + x[5] + x[6] + x[7] - x[0] - x[1] - x[2] - x[3]);
      gi.normalise();
      gj.normalise();
      vec3_t gk = gi.cross(gj);
      double Li = 0;
      double Lj = 0;
      double Lk = 0;
      for (int i = 0; i < 8; ++i) {
        vec3_t Dx = x[i] - xc;
        Li = max(2*fabs(Dx*gi), Li);
        Lj = max(2*fabs(Dx*gj), Lj);
        Lk = max(2*fabs(Dx*gk), Lk);
      }
      int Ni = max(5, int(Li/h));
      int Nj = max(5, int(Lj/h));
      int Nk = max(5, int(Lk/h));
      double hi = Li/Ni;
      double hj = Lj/Nj;
      double hk = Lk/Nk;
      vec3_t x0 = xc - 0.5*(Li*gi + Lj*gj + Lk*gk);

      {
        QString bc = boundaryCode(id_cell, 4);
        patch.bx = bc;
        if (bc == "0") {
          patch.sl_i1 = m_OverlapLayers;
          Ni += m_OverlapLayers;
          Li += m_OverlapLayers*hi;
          x0 -= m_OverlapLayers*hi*gi;
        } else {
          patch.sl_i1 = 0;
        }
      }
      {
        QString bc = boundaryCode(id_cell, 5);
        patch.bX = bc;
        if (bc == "0") {
          patch.sl_i2 = m_OverlapLayers;
          Ni += m_OverlapLayers;
          Li += m_OverlapLayers*hi;
        } else {
          patch.sl_i2 = 0;
        }
      }
      {
        QString bc = boundaryCode(id_cell, 0);
        patch.by = bc;
        if (bc == "0") {
          patch.sl_j1 = m_OverlapLayers;
          Nj += m_OverlapLayers;
          Lj += m_OverlapLayers*hj;
          x0 -= m_OverlapLayers*hj*gj;
        } else {
          patch.sl_j1 = 0;
        }
      }
      {
        QString bc = boundaryCode(id_cell, 1);
        patch.bY = bc;
        if (bc == "0") {
          patch.sl_j2 = m_OverlapLayers;
          Nj += m_OverlapLayers;
          Lj += m_OverlapLayers*hj;
        } else {
          patch.sl_j2 = 0;
        }
      }
      {
        QString bc = boundaryCode(id_cell, 3);
        patch.bz = bc;
        if (bc == "0") {
          patch.sl_k1 = m_OverlapLayers;
          Nk += m_OverlapLayers;
          Lk += m_OverlapLayers*hk;
          x0 -= m_OverlapLayers*hk*gk;
        } else {
          patch.sl_k1 = 0;
        }
      }
      {
        QString bc = boundaryCode(id_cell,2);
        patch.bZ = bc;
        if (bc == "0") {
          patch.sl_k2 = m_OverlapLayers;
          Nk += m_OverlapLayers;
          Lk += m_OverlapLayers*hk;
        } else {
          patch.sl_k2 = 0;
        }
      }

      patch.gi = gi;
      patch.gj = gj;
      patch.Li = Li;
      patch.Lj = Lj;
      patch.Lk = Lk;
      patch.Ni = Ni;
      patch.Nj = Nj;
      patch.Nk = Nk;
      patch.x0 = x0;

      patch.fx = "fx";
      patch.fy = "fy";
      patch.fz = "fz";
      patch.s  = "0";

      m_CartPatches[m_CellToCartPatch[id_cell]] = patch;
    }
  }
}

void DrNumWriter::writeCartPatches(QTextStream &s)
{
  foreach (cart_patch_t patch, m_CartPatches) {
    s << "1001\n";
    s << "{\n";
    s << "  " << patch.x0[0] << " " << patch.x0[1] << " " << patch.x0[2] << "\n";
    s << "  " << patch.gi[0] << " " << patch.gi[1] << " " << patch.gi[2] << "\n";
    s << "  " << patch.gj[0] << " " << patch.gj[1] << " " << patch.gj[2] << "\n";
    s << "  1.0\n";
    s << "  " << patch.Ni << " " << patch.Nj << " " << patch.Nk << "\n";
    s << "  " << patch.sl_i1 << " " << patch.sl_i2 << " " << patch.sl_j1 << " " << patch.sl_j2 << " " << patch.sl_k1 << " " << patch.sl_k2 << "\n";
    s << "  " << patch.Li << " " << patch.Lj << " " << patch.Lk << "\n";
    s << "  " << patch.fx << " " << patch.fy << " " << patch.fz << "\n";
    s << "  " << patch.bx << " " << patch.bX << " " << patch.by << " " << patch.bY << " " << patch.bz << " " << patch.bZ << "\n";
    s << "  " << patch.s << "\n";
    s << "}\n";
    s << "\n";
  }
}

void DrNumWriter::operate()
{
  try {
    readOutputDirectory();
    if (isValid()) {
      QString p1 = getFileName();
      QString p2 = p1 + "/patches";
      QDir d1(p1);
      QDir d2(p2);
      if (!d1.exists()) {
        EG_BUG;
      }
      if (!d2.exists()) {
        d1.mkdir("patches");
        d2 = QDir(p2);
      }
      d1 = d2;
      p1 = p2;
      p2 = p1 + "/complex";
      d2 = QDir(p2);
      if (!d2.exists()) {
        d1.mkdir("complex");
      };
      m_PatchFile   = getFileName() + "/patches/from_enGrid";
      m_ComplexPath = getFileName() + "/patches/complex";
      if (!QDir(m_ComplexPath).exists()) {
        EG_BUG;
      }
      QFile file(m_PatchFile);
      file.open(QIODevice::WriteOnly);
      QTextStream s(&file);
      computeMeshDensity();
      createAllCartPatches();
      writeCartPatches(s);
      s << "0\n";
    }
  } catch (Error err) {
    err.display();
  }
}
