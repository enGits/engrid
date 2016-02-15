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

#include "guicreatehexshell.h"
#include "guimainwindow.h"

GuiCreateHexShell::GuiCreateHexShell()
{
  connect(m_Ui.spinBoxNx, SIGNAL(valueChanged(int)), this, SLOT(updateNumberOfCells(int)));
  connect(m_Ui.spinBoxNy, SIGNAL(valueChanged(int)), this, SLOT(updateNumberOfCells(int)));
  connect(m_Ui.spinBoxNz, SIGNAL(valueChanged(int)), this, SLOT(updateNumberOfCells(int)));
  connect(m_Ui.spinBoxLx1, SIGNAL(valueChanged(int)), this, SLOT(updateNumberOfCells(int)));
  connect(m_Ui.spinBoxLx2, SIGNAL(valueChanged(int)), this, SLOT(updateNumberOfCells(int)));
  connect(m_Ui.spinBoxLy1, SIGNAL(valueChanged(int)), this, SLOT(updateNumberOfCells(int)));
  connect(m_Ui.spinBoxLy2, SIGNAL(valueChanged(int)), this, SLOT(updateNumberOfCells(int)));
  connect(m_Ui.spinBoxLz1, SIGNAL(valueChanged(int)), this, SLOT(updateNumberOfCells(int)));
  connect(m_Ui.spinBoxLz2, SIGNAL(valueChanged(int)), this, SLOT(updateNumberOfCells(int)));
}

void GuiCreateHexShell::ijkCell(int idx, int &i, int &j, int &k)
{
  i = idx/(m_NumJCells*m_NumKCells);
  int rest = idx - i*(m_NumJCells*m_NumKCells);
  j = rest/m_NumKCells;
  k = rest - j*m_NumKCells;
}

void GuiCreateHexShell::ijkNode(int idx, int &i, int &j, int &k)
{
  i = idx/(m_NumJNodes*m_NumKNodes);
  int rest = idx - i*(m_NumJNodes*m_NumKNodes);
  j = rest/m_NumKNodes;
  k = rest - j*m_NumKNodes;
}

void GuiCreateHexShell::updateNumberOfCells(int)
{
  int num_i1 = m_Ui.spinBoxNx->value();
  int num_i2 = num_i1 + m_Ui.spinBoxLx1->value() + m_Ui.spinBoxLx2->value();
  int num_j1 = m_Ui.spinBoxNy->value();
  int num_j2 = num_j1 + m_Ui.spinBoxLy1->value() + m_Ui.spinBoxLy2->value();
  int num_k1 = m_Ui.spinBoxNz->value();
  int num_k2 = num_k1 + m_Ui.spinBoxLz1->value() + m_Ui.spinBoxLz2->value();
  m_TotalNumberOfCells = num_i2*num_j2*num_k2 - num_i1*num_j1*num_k1;
  int num_cells = m_TotalNumberOfCells;
  QList<int> numbers;
  while (num_cells > 0) {
    numbers << num_cells % 1000;
    num_cells /= 1000;
  }
  QString txt;
  txt.setNum(numbers.last());
  numbers.pop_back();
  while (numbers.size() > 0) {
    QString num;
    num.setNum(numbers.last());
    while (num.size() < 3) {
      num = "0" + num;
    }
    txt += "," + num;
    numbers.pop_back();
  }
  m_Ui.labelTotalNumberOfCells->setText(txt);
}

void GuiCreateHexShell::before()
{
  EG_VTKDCN(vtkDoubleArray, cl, m_Grid, "node_meshdensity_desired");
  m_H = 1e99;
  for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
    m_H = min(m_H, cl->GetValue(id_node));
  }
  vec3_t x1( 1e99,  1e99,  1e99);
  vec3_t x2(-1e99, -1e99, -1e99);
  for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
    vec3_t x;
    m_Grid->GetPoint(id_node, x.data());
    for (int i = 0; i < 3; ++i) {
      x1[i] = min(x1[i], x[i]);
      x2[i] = max(x2[i], x[i]);
    }
  }
  x1 -= 5*vec3_t(m_H, m_H, m_H);
  x2 += 5*vec3_t(m_H, m_H, m_H);
  m_Ui.spinBoxNx->setValue(int((x2[0] - x1[0])/m_H));
  m_Ui.spinBoxNy->setValue(int((x2[1] - x1[1])/m_H));
  m_Ui.spinBoxNz->setValue(int((x2[2] - x1[2])/m_H));
  setDouble(x1[0], m_Ui.lineEditX1);
  setDouble(x1[1], m_Ui.lineEditY1);
  setDouble(x1[2], m_Ui.lineEditZ1);
  setDouble(x2[0], m_Ui.lineEditX2);
  setDouble(x2[1], m_Ui.lineEditY2);
  setDouble(x2[2], m_Ui.lineEditZ2);
  updateNumberOfCells();
}

void GuiCreateHexShell::createGridWithNodes(vtkUnstructuredGrid *grid)
{
  int num_nodes = 0;
  m_NodeIDs.fill(-1, m_NumINodes*m_NumJNodes*m_NumKNodes);
  vtkIdType id_node = 0;
  for (int i = 0; i < m_NumINodes; ++i) {
    for (int j = 0; j < m_NumJNodes; ++j) {
      for (int k = 0; k < m_NumKNodes; ++k) {
        bool in_core = true;
        in_core = in_core && i > m_Ui.spinBoxLx1->value();
        in_core = in_core && i < m_NumINodes - m_Ui.spinBoxLx2->value() - 1;
        in_core = in_core && j > m_Ui.spinBoxLy1->value();
        in_core = in_core && j < m_NumJNodes - m_Ui.spinBoxLy2->value() - 1;
        in_core = in_core && k > m_Ui.spinBoxLz1->value();
        in_core = in_core && k < m_NumKNodes - m_Ui.spinBoxLz2->value() - 1;
        if (!in_core) {
          ++num_nodes;
        }
      }
    }
  }
  allocateGrid(grid, m_TotalNumberOfCells, num_nodes);
  for (int i = 0; i < m_NumINodes; ++i) {
    for (int j = 0; j < m_NumJNodes; ++j) {
      for (int k = 0; k < m_NumKNodes; ++k) {
        bool in_core = true;
        in_core = in_core && i > m_Ui.spinBoxLx1->value();
        in_core = in_core && i < m_NumINodes - m_Ui.spinBoxLx2->value() - 1;
        in_core = in_core && j > m_Ui.spinBoxLy1->value();
        in_core = in_core && j < m_NumJNodes - m_Ui.spinBoxLy2->value() - 1;
        in_core = in_core && k > m_Ui.spinBoxLz1->value();
        in_core = in_core && k < m_NumKNodes - m_Ui.spinBoxLz2->value() - 1;
        if (!in_core) {
          m_NodeIDs[indexNode(i,j,k)] = id_node;
          grid->GetPoints()->SetPoint(id_node, x(i), y(j), z(k));
          ++id_node;
        }
      }
    }
  }
}

void GuiCreateHexShell::createHexCells(vtkUnstructuredGrid *grid)
{
  EG_VTKDCC(vtkIntArray, cell_code, grid, "cell_code");
  m_CellIDs.fill(-1, m_NumICells*m_NumJCells*m_NumKCells);
  for (int i = 0; i < m_NumICells; ++i) {
    for (int j = 0; j < m_NumJCells; ++j) {
      for (int k = 0; k < m_NumKCells; ++k) {
        bool in_core = true;
        in_core = in_core && i >= m_Ui.spinBoxLx1->value();
        in_core = in_core && i < m_NumICells - m_Ui.spinBoxLx2->value();
        in_core = in_core && j >= m_Ui.spinBoxLy1->value();
        in_core = in_core && j < m_NumJCells - m_Ui.spinBoxLy2->value();
        in_core = in_core && k >= m_Ui.spinBoxLz1->value();
        in_core = in_core && k < m_NumKCells - m_Ui.spinBoxLz2->value();
        if (!in_core) {
          EG_VTKSP(vtkIdList, pts);
          pts->SetNumberOfIds(8);
          pts->SetId(0, m_NodeIDs[indexNode(i  ,j  ,k  )]);
          pts->SetId(1, m_NodeIDs[indexNode(i+1,j  ,k  )]);
          pts->SetId(2, m_NodeIDs[indexNode(i+1,j+1,k  )]);
          pts->SetId(3, m_NodeIDs[indexNode(i  ,j+1,k  )]);
          pts->SetId(4, m_NodeIDs[indexNode(i  ,j  ,k+1)]);
          pts->SetId(5, m_NodeIDs[indexNode(i+1,j  ,k+1)]);
          pts->SetId(6, m_NodeIDs[indexNode(i+1,j+1,k+1)]);
          pts->SetId(7, m_NodeIDs[indexNode(i  ,j+1,k+1)]);
          for (int i_pts = 0; i_pts < 8; ++i_pts) {
            if (pts->GetId(i_pts) == -1) {
              EG_BUG;
            }
          }
          m_CellIDs[indexCell(i,j,k)] = grid->InsertNextCell(VTK_HEXAHEDRON, pts);
          if (m_CellIDs[indexCell(i,j,k)] >= m_TotalNumberOfCells) {
            EG_BUG;
          }
          cell_code->SetValue(m_CellIDs[indexCell(i,j,k)], 0);
        }
      }
    }
  }
}

void GuiCreateHexShell::defineBoundaryCodes()
{
  QSet<int> bcs = mainWindow()->getAllBoundaryCodes();
  m_InnerBC = 1;
  foreach (int bc, bcs) {
    m_InnerBC = max(bc, m_InnerBC);
  }
  ++m_InnerBC;
  m_OuterBC = m_InnerBC + 1;
  mainWindow()->setBC(m_InnerBC, BoundaryCondition("HexShellInside", "patch", m_InnerBC));
  mainWindow()->setBC(m_OuterBC, BoundaryCondition("HexShellOutside", "patch", m_OuterBC));
}

void GuiCreateHexShell::createOuterBoundary(vtkUnstructuredGrid *grid)
{
  EG_VTKSP(vtkUnstructuredGrid, new_grid);
  allocateGrid(new_grid, grid->GetNumberOfCells() + 2*(m_NumICells*m_NumJCells + m_NumICells*m_NumKCells + m_NumJCells*m_NumKCells), grid->GetNumberOfPoints());
  makeCopyNoAlloc(grid, new_grid);
  EG_VTKDCC(vtkIntArray, cell_code, new_grid, "cell_code");

  // left boundary
  {
    int i = 0;
    for (int j = 0; j < m_NumJCells; ++j) {
      for (int k = 0; k < m_NumKCells; ++k) {
        vtkIdType id_cell = m_CellIDs[indexCell(i,j,k)];
        if (id_cell == -1) {
          EG_BUG;
        }
        QVector<vtkIdType> pts;
        getFaceOfCell(new_grid, id_cell, 4, pts);
        vtkIdType id_new_cell = new_grid->InsertNextCell(VTK_QUAD, 4, pts.data());
        cell_code->SetValue(id_new_cell, m_OuterBC);
      }
    }
  }

  // right boundary
  {
    int i = m_NumICells - 1;
    for (int j = 0; j < m_NumJCells; ++j) {
      for (int k = 0; k < m_NumKCells; ++k) {
        vtkIdType id_cell = m_CellIDs[indexCell(i,j,k)];
        if (id_cell == -1) {
          EG_BUG;
        }
        QVector<vtkIdType> pts;
        getFaceOfCell(new_grid, id_cell, 5, pts);
        vtkIdType id_new_cell = new_grid->InsertNextCell(VTK_QUAD, 4, pts.data());
        cell_code->SetValue(id_new_cell, m_OuterBC);
      }
    }
  }

  // front boundary
  {
    for (int i = 0; i < m_NumICells; ++i) {
      int j = 0;
      for (int k = 0; k < m_NumKCells; ++k) {
        vtkIdType id_cell = m_CellIDs[indexCell(i,j,k)];
        if (id_cell == -1) {
          EG_BUG;
        }
        QVector<vtkIdType> pts;
        getFaceOfCell(new_grid, id_cell, 2, pts);
        vtkIdType id_new_cell = new_grid->InsertNextCell(VTK_QUAD, 4, pts.data());
        cell_code->SetValue(id_new_cell, m_OuterBC);
      }
    }
  }

  // back boundary
  {
    for (int i = 0; i < m_NumICells; ++i) {
      int j = m_NumJCells - 1;
      for (int k = 0; k < m_NumKCells; ++k) {
        vtkIdType id_cell = m_CellIDs[indexCell(i,j,k)];
        if (id_cell == -1) {
          EG_BUG;
        }
        QVector<vtkIdType> pts;
        getFaceOfCell(new_grid, id_cell, 3, pts);
        vtkIdType id_new_cell = new_grid->InsertNextCell(VTK_QUAD, 4, pts.data());
        cell_code->SetValue(id_new_cell, m_OuterBC);
      }
    }
  }

  // bottom boundary
  {
    for (int i = 0; i < m_NumICells; ++i) {
      for (int j = 0; j < m_NumJCells; ++j) {
        int k = 0;
        vtkIdType id_cell = m_CellIDs[indexCell(i,j,k)];
        if (id_cell == -1) {
          EG_BUG;
        }
        QVector<vtkIdType> pts;
        getFaceOfCell(new_grid, id_cell, 0, pts);
        vtkIdType id_new_cell = new_grid->InsertNextCell(VTK_QUAD, 4, pts.data());
        cell_code->SetValue(id_new_cell, m_OuterBC);
      }
    }
  }

  // top boundary
  {
    for (int i = 0; i < m_NumICells; ++i) {
      for (int j = 0; j < m_NumJCells; ++j) {
        int k = m_NumKCells - 1;
        vtkIdType id_cell = m_CellIDs[indexCell(i,j,k)];
        if (id_cell == -1) {
          EG_BUG;
        }
        QVector<vtkIdType> pts;
        getFaceOfCell(new_grid, id_cell, 1, pts);
        vtkIdType id_new_cell = new_grid->InsertNextCell(VTK_QUAD, 4, pts.data());
        cell_code->SetValue(id_new_cell, m_OuterBC);
      }
    }
  }

  makeCopy(new_grid, grid);
}

void GuiCreateHexShell::createInnerBoundary(vtkUnstructuredGrid *grid)
{
  int num_new_nodes = 0;
  int num_new_faces = 0;
  MeshPartition part(grid, true);
  QVector<bool> copy_cell(grid->GetNumberOfCells(), true);
  for (vtkIdType id_cell = 0; id_cell < grid->GetNumberOfCells(); ++id_cell) {
    if (isVolume(id_cell, grid)) {
      bool face_found = false;
      for (int i = 0; i < part.c2cGSize(id_cell); ++i) {
        vtkIdType id_neigh = part.c2cGG(id_cell, i);
        if (id_neigh == -1) {
          if (face_found) {
            EG_BUG;
          } else {
            num_new_faces += 4;
            num_new_nodes += 1;
            face_found = true;
            copy_cell[id_cell] = false;
          }
        }
      }
    }
  }
  EG_VTKSP(vtkUnstructuredGrid, new_grid);
  allocateGrid(new_grid, grid->GetNumberOfCells() + num_new_faces, grid->GetNumberOfPoints() + num_new_nodes);
  for (vtkIdType id_node = 0; id_node < grid->GetNumberOfPoints(); ++id_node) {
    vec3_t x;
    grid->GetPoint(id_node, x.data());
    new_grid->GetPoints()->SetPoint(id_node, x.data());
  }
  for (vtkIdType id_cell = 0; id_cell < grid->GetNumberOfCells(); ++id_cell) {
    if (copy_cell[id_cell]) {
      vtkIdType id_new_cell = copyCell(grid, id_cell, new_grid);
      copyCellData(grid, id_cell, new_grid, id_new_cell);
    }
  }
  EG_VTKDCC(vtkIntArray, cell_code, new_grid, "cell_code");
  vtkIdType id_new_node = grid->GetNumberOfPoints();
  for (vtkIdType id_cell = 0; id_cell < grid->GetNumberOfCells(); ++id_cell) {
    if (isVolume(id_cell, grid)) {
      int i_replace_face = -1;
      QVector<QVector<vtkIdType> > triangles;
      for (int i = 0; i < part.c2cGSize(id_cell); ++i) {
        vtkIdType id_neigh = part.c2cGG(id_cell, i);
        if (id_neigh == -1) {
          i_replace_face = i;
          QVector<vtkIdType> pts;
          getFaceOfCell(grid, id_cell, i, pts);
          if (pts.size() != 4) {
            EG_BUG;
          }
          vec3_t xc(0,0,0);
          foreach (vtkIdType id_node, pts) {
            vec3_t x;
            grid->GetPoint(id_node, x.data());
            xc += x;
          }
          xc *= 0.25;
          new_grid->GetPoints()->SetPoint(id_new_node, xc.data());
          pts.append(pts.first());
          for (int j = 0; j < 4; ++j) {
            QVector<vtkIdType> triangle(3);
            triangle[0] = pts[j];
            triangle[1] = pts[j+1];
            triangle[2] = id_new_node;
            vtkIdType id_new_face = new_grid->InsertNextCell(VTK_TRIANGLE, 3, triangle.data());
            cell_code->SetValue(id_new_face, m_InnerBC);
            triangles.push_back(triangle);
          }
          ++id_new_node;
        }
      }
      if (i_replace_face != -1) {
        EG_VTKSP(vtkIdList, stream);
        stream->SetNumberOfIds(42);
        vtkIdType id = 0;
        stream->SetId(id++, 9);
        for (int i = 0; i < 6; ++i) {
          if (i != i_replace_face) {
            stream->SetId(id++, 4);
            QVector<vtkIdType> pts;
            getFaceOfCell(grid, id_cell, i, pts);
            if (pts.size() != 4) {
              EG_BUG;
            }
            foreach (vtkIdType id_node, pts) {
              stream->SetId(id++, id_node);
            }
          }
        }
        foreach (QVector<vtkIdType> triangle, triangles) {
          stream->SetId(id++, 3);
          foreach (vtkIdType id_node, triangle) {
            stream->SetId(id++, id_node);
          }
        }
        vtkIdType id_new_cell = new_grid->InsertNextCell(VTK_POLYHEDRON, stream);
        copyCellData(grid, id_cell, new_grid, id_new_cell);
      }
    }
  }
  makeCopy(new_grid, grid);
}

void GuiCreateHexShell::createSourceBox()
{
  QString xml_text = mainWindow()->getXmlSection("engrid/sources");
  {
    QTextStream s(&xml_text);
    s << "box: HexShell; " << m_X0 << "; " << m_Y0 << "; " << m_Z0 << "; ";
    s << m_X0 + m_NumICells*m_Dx << "; ";
    s << m_Y0 + m_NumJCells*m_Dy << "; ";
    s << m_Z0 + m_NumKCells*m_Dz << "; ";
    s << 1.5*min(m_Dx, min(m_Dy, m_Dz)) << ";\n";
  }
  mainWindow()->setXmlSection("engrid/sources", xml_text);
}

void GuiCreateHexShell::operate()
{
  m_NumICells =  m_Ui.spinBoxNx->value() + m_Ui.spinBoxLx1->value() + m_Ui.spinBoxLx2->value();
  m_NumJCells =  m_Ui.spinBoxNy->value() + m_Ui.spinBoxLy1->value() + m_Ui.spinBoxLy2->value();
  m_NumKCells =  m_Ui.spinBoxNz->value() + m_Ui.spinBoxLz1->value() + m_Ui.spinBoxLz2->value();
  m_NumINodes = m_NumICells + 1;
  m_NumJNodes = m_NumJCells + 1;
  m_NumKNodes = m_NumKCells + 1;
  double x1 = m_Ui.lineEditX1->text().toDouble();
  double x2 = m_Ui.lineEditX2->text().toDouble();
  double y1 = m_Ui.lineEditY1->text().toDouble();
  double y2 = m_Ui.lineEditY2->text().toDouble();
  double z1 = m_Ui.lineEditZ1->text().toDouble();
  double z2 = m_Ui.lineEditZ2->text().toDouble();
  m_Dx = (x2 - x1)/m_Ui.spinBoxNx->value();
  m_Dy = (y2 - y1)/m_Ui.spinBoxNy->value();
  m_Dz = (z2 - z1)/m_Ui.spinBoxNz->value();
  m_X0 = x1 - m_Ui.spinBoxLx1->value()*m_Dx;
  m_Y0 = y1 - m_Ui.spinBoxLy1->value()*m_Dy;
  m_Z0 = z1 - m_Ui.spinBoxLz1->value()*m_Dz;
  EG_VTKSP(vtkUnstructuredGrid, shell_grid);
  createGridWithNodes(shell_grid);
  createHexCells(shell_grid);
  defineBoundaryCodes();
  createOuterBoundary(shell_grid);
  createInnerBoundary(shell_grid);
  MeshPartition shell_part(shell_grid, true);
  m_Part.addPartition(shell_part);
  createSourceBox();
  mainWindow()->updateBoundaryCodes(true);
}
