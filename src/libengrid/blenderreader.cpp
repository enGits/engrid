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

#include "blenderreader.h"
#include "guimainwindow.h"
#include "pointfinder.h"

BlenderReader::BlenderReader()
{
  setFormat("Blender/Engrid files(*.begc *.BEGC)");
  getSet("General", "tolerance for importing geometries (% of smallest edge length)", 0.1, m_RelativeTolerance);
  m_RelativeTolerance *= 1e-2;
  m_Append = false;
}

void BlenderReader::operate()
{
  try {
    QFileInfo file_info(GuiMainWindow::pointer()->getFilename());
    readInputFileName(file_info.completeBaseName() + ".begc", false);
    if (isValid()) {
      // read raw data from exported file
      QFile file(getFileName());
      file.open(QIODevice::ReadOnly | QIODevice::Text);
      QTextStream f(&file);
      QList<vec3_t> rnodes;
      QList<QVector<int> > rfaces;
      int num_parts;
      f >> num_parts;
      QVector<QString> part_name(num_parts);
      for (int i_part = 0; i_part < num_parts; ++i_part) {
        f >> part_name[i_part];
      }
      QVector<QString> sorted_part_name = part_name;
      qSort(sorted_part_name);
      QVector<int> part_bc(part_name.size());
      for (int i_part = 0; i_part < num_parts; ++i_part) {
        part_bc[i_part] = sorted_part_name.indexOf(part_name[i_part]) + 1;
      }
      for (int i_part = 0; i_part < num_parts; ++i_part) {
        int num_nodes, num_faces;
        f >> num_nodes >> num_faces;
        for (int i = 0; i < num_nodes; ++i) {
          vec3_t x;
          f >> x[0] >> x[1] >> x[2];
          rnodes.push_back(x);
        }
        for (int i = 0; i < num_faces; ++i) {
          int N;
          f >> N;
          QVector<int> face(N+1);
          face[0] = i_part;
          for (int j = 0; j < N; ++j) {
            f >> face[j+1];
          }
          rfaces.push_back(face);
        }
      }
      QVector<vec3_t> nodes(rnodes.size());
      qCopy(rnodes.begin(), rnodes.end(), nodes.begin());
      QVector<QVector<int> > faces(rfaces.size());
      qCopy(rfaces.begin(), rfaces.end(), faces.begin());

      // find smallest edge length
      double L = 1e99;
      foreach (QVector<int> face, faces) {
        for (int i = 1; i < face.size(); ++i) {
          int n1 = face[i];
          int n2 = face[1];
          if (i < face.size() - 1) {
            n2 = face[i+1];
          }
          double l = (nodes[n1] - nodes[n2]).abs();
          L = min(l, L);
        }
      }

      cout << "smallest edge length is " << L << endl;

      // delete duplicate nodes
      PointFinder finder;
      finder.setPoints(nodes);
      QList<vec3_t> non_dup;
      QVector<int> o2n(nodes.size());
      int num_non_dup = 0;
      for (int i = 0; i < nodes.size(); ++i) {
        o2n[i] = num_non_dup;
        bool dup = false;
        QVector<int> close_points;
        finder.getClosePoints(nodes[i], close_points);
        foreach (int j, close_points) {
          if (i > j) {
            double l = (nodes[i] - nodes[j]).abs();
            if (l < m_RelativeTolerance*L || l == 0) {
              o2n[i] = o2n[j];
              dup = true;
              break;
            }
          }
        }
        if (!dup) {
          non_dup.push_back(nodes[i]);
          ++num_non_dup;
        }
      }

      EG_VTKSP(vtkUnstructuredGrid, new_grid);

      allocateGrid(new_grid, faces.size(), non_dup.size());
      EG_VTKDCC(vtkIntArray, cell_code, new_grid, "cell_code");
      EG_VTKDCC(vtkIntArray, orgdir, new_grid, "cell_orgdir");
      EG_VTKDCC(vtkIntArray, voldir, new_grid, "cell_voldir");
      EG_VTKDCC(vtkIntArray, curdir, new_grid, "cell_curdir");
      vtkIdType id_node = 0;
      foreach (vec3_t x, non_dup) {
        new_grid->GetPoints()->SetPoint(id_node, x.data());
        ++id_node;
      }

      foreach (QVector<int> face, faces) {
        if (face.size() == 4) {
          vtkIdType pts[3];
          pts[0] = o2n[face[1]];
          pts[1] = o2n[face[2]];
          pts[2] = o2n[face[3]];
          vtkIdType id_cell = new_grid->InsertNextCell(VTK_TRIANGLE, 3, pts);
          cell_code->SetValue(id_cell, part_bc[face[0]]);
          orgdir->SetValue(id_cell, 0);
          voldir->SetValue(id_cell, 0);
          curdir->SetValue(id_cell, 0);
        }
        if (face.size() == 5) {
          vtkIdType pts[4];
          pts[0] = o2n[face[1]];
          pts[1] = o2n[face[2]];
          pts[2] = o2n[face[3]];
          pts[3] = o2n[face[4]];
          vtkIdType id_cell = new_grid->InsertNextCell(VTK_QUAD, 4, pts);
          cell_code->SetValue(id_cell, part_bc[face[0]]);
          orgdir->SetValue(id_cell, 0);
          voldir->SetValue(id_cell, 0);
          curdir->SetValue(id_cell, 0);
        }
      }

      if (m_Append) {
        EG_BUG;
        MeshPartition new_part(new_grid);
        new_part.setAllCells();
        m_Part.addPartition(new_part);
      } else {
        makeCopy(new_grid, m_Grid);
      }

      UpdateNodeIndex(m_Grid);
      UpdateCellIndex(m_Grid);

      // check and set the boundary names if required
      int update_required = true;
      QSet<int> old_bcs = GuiMainWindow::pointer()->getAllBoundaryCodes();
      if (old_bcs.size() == part_name.size()) {
        QSet<QString> old_names;
        foreach (int bc, old_bcs) {
          old_names.insert(GuiMainWindow::pointer()->getBC(bc).getName());
        }
        QSet<QString> new_names;
        foreach (QString name, part_name) {
          new_names.insert(name);
        }
        if (old_names == new_names) {
          update_required = false;
          cout << "no update required" << endl;
        }
      }
      if (update_required) {
        GuiMainWindow::pointer()->resetXmlDoc();
        GuiMainWindow::pointer()->clearBCs();
        for (int i_part = 0; i_part < part_name.size(); ++i_part) {
          GuiMainWindow::pointer()->setBC(part_bc[i_part], BoundaryCondition(part_name[i_part], "patch", part_bc[i_part]));
        }
      }

    }
  } catch (Error err) {
    err.display();
  }
}
