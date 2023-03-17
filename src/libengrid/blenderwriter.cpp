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
#include "blenderwriter.h"
#include "guimainwindow.h"

BlenderWriter::BlenderWriter()
{
  setFormat("Blender/Engrid files(*.begc)");
}

void BlenderWriter::operate()
{
  try {
    QFileInfo file_info(GuiMainWindow::pointer()->getFilename());
    readOutputFileName(file_info.completeBaseName() + ".begc");
    if (isValid()) {
      QFile file(getFileName());
      file.open(QIODevice::WriteOnly | QIODevice::Text);
      QTextStream f(&file);
      
      // write number of objects
      QSet<int> bcs = GuiMainWindow::pointer()->getAllBoundaryCodes();
      int N_BoundaryCodes = bcs.size();
      qWarning()<<"N_BoundaryCodes="<<N_BoundaryCodes;
      f << N_BoundaryCodes << "\n";
      
      // write names of objects
      foreach(int bc, bcs) {
        QString BCname = getBC(bc).getName();
        qWarning()<<"BCname="<<BCname;
        f << BCname << "\n";
      }
      
      // write objects
      int offset = 0;
      foreach(int bc, bcs) {
        
        QVector <vtkIdType> object_cells;
        QSet <int> object_bc_set;
        object_bc_set.insert(bc);
        getSurfaceCells(object_bc_set, object_cells, m_Grid);
        
        EG_VTKSP(vtkUnstructuredGrid,SubGrid);
        getSubGrid(m_Grid,object_cells,SubGrid);
        
        int N_verts = SubGrid->GetNumberOfPoints();
        int N_faces = SubGrid->GetNumberOfCells();
        qWarning()<<"N_verts="<<N_verts << " N_faces=" << N_faces;
        f << N_verts << " " << N_faces << "\n";
        
        // prepare to write vertices and faces
        QVector <vec3_t> subvertices_data(N_verts);
        QMap <vtkIdType, vtkIdType> subvertices_idx;
        int idx = 0;
        QVector <bool> vertex_written(m_Grid->GetNumberOfPoints(), false);
        foreach (vtkIdType cellId, object_cells) {
          EG_GET_CELL(cellId, m_Grid);
          for (int i = 0; i < num_pts; ++i) {
            if(!vertex_written[pts[i]]) {
              vec3_t x;
              m_Grid->GetPoints()->GetPoint(pts[i], x.data());
              if(idx<0 || idx>=N_verts) EG_BUG;
              subvertices_data[idx] = x;
              subvertices_idx[pts[i]]=idx;
              idx++;
              vertex_written[pts[i]] = true;
            }
          }
        }
        
        // write vertices
        foreach(vec3_t x, subvertices_data) {
              f << x[0] << ' ' << x[1] << ' ' << x[2] << '\n';
        }
        
        // write faces
        foreach (vtkIdType cellId, object_cells) {
          EG_GET_CELL(cellId, m_Grid);
          f << num_pts;
          for (int i = 0; i < num_pts; ++i) {
            f << ' ' << offset + subvertices_idx[pts[i]];
          }
          f << '\n';
        }
        
        offset+=N_verts;
        
      }// end of loop through objects
      
      file.close();
      
    }// end of if(isValid)
  } catch (Error err) {
    err.display();
  }
}
