//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008 Oliver Gloth                                          +
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
#include "vtkEgEliminateShortEdges.h"
#include "geometrytools.h"

vtkStandardNewMacro(vtkEgEliminateShortEdges);

vtkEgEliminateShortEdges::vtkEgEliminateShortEdges()
{
  max_ratio = 1000.0;
  max_length = 1e99;
};

void vtkEgEliminateShortEdges::CheckEdges()
{
  new_node_idx.fill(-1, m_Input->GetNumberOfPoints());
  delete_cell.fill(false, m_Input->GetNumberOfCells());
  QVector<bool> marked(m_Input->GetNumberOfPoints(), false);
  for (vtkIdType cellId = 0; cellId < m_Input->GetNumberOfCells(); ++cellId) {
    int cellType = m_Input->GetCellType(cellId);
    if ((cellType == VTK_TRIANGLE) || (cellType == VTK_QUAD)) {
      vtkIdType *pts, Npts;
      m_Input->GetCellPoints(cellId, Npts, pts);
      double L_av = 0;
      vector<vec3_t> x(Npts);
      int N_del = 0;
      double Lmin = 1e99;
      int imin = 0;
      for (int i = 0; i < Npts; ++i) {
        vtkIdType p1 = pts[i];
        vtkIdType p2;
        if (i < Npts-1) p2 = pts[i+1];
        else            p2 = pts[0];
        if (p1 == p2) {
          EG_ERR_RETURN("bug encountered");
        };
        vec3_t x1, x2;
        m_Input->GetPoints()->GetPoint(p1, x1.data());
        m_Input->GetPoints()->GetPoint(p2, x2.data());
        double L = (x2-x1).abs();
        L_av += L;
        if (L < Lmin) {
          Lmin = L;
          imin = i;
        };
        x[i] = x1;
        if (marked[pts[i]]) {
          ++N_del;
        };
      };
      if (N_del == 0) {
        double A = GeometryTools::triArea(x[0],x[1],x[2]);
        if (Npts == 4) {
          A = GeometryTools::quadArea(x[0],x[1],x[2],x[3]);
        };
        L_av /= Npts;
        for (int i = 0; i < Npts; ++i) {
          vtkIdType p1 = pts[i];
          vtkIdType p2;
          if (i < Npts-1) p2 = pts[i+1];
          else            p2 = pts[0];
          vec3_t x1, x2;
          m_Input->GetPoints()->GetPoint(p1, x1.data());
          m_Input->GetPoints()->GetPoint(p2, x2.data());
          double L = (x2-x1).abs();
          bool delete_edge = false;
          if (L == 0) {
            delete_edge = true;
          } else if ((L_av*L_av/A > max_ratio) && (i == imin)) {
            if (L < max_length) {
              delete_edge = true;
            };
          };
          if (delete_edge) {
            new_node_idx[p1] = p2;
            marked[p1] = true;
            marked[p2] = true;
            if (cellType != VTK_TRIANGLE) {
              EG_ERR_RETURN("The present configuration cannot be handled yet.");
            };
            delete_cell[cellId] = true;
            break;
          };
        };
      };
    };
  };
};

void vtkEgEliminateShortEdges::CheckCells()
{
  for (vtkIdType cellId = 0; cellId < m_Input->GetNumberOfCells(); ++cellId) {
    int cellType = m_Input->GetCellType(cellId);
    if (cellType == VTK_TRIANGLE) {
      vtkIdType *pts, Npts;
      m_Input->GetCellPoints(cellId, Npts, pts);
      for (int i = 0; i < Npts; ++i) {
        vtkIdType p1 = pts[i];
        vtkIdType p2;
        if (i < Npts-1) p2 = pts[i+1];
        else            p2 = pts[0];
        if (p1 == p2) {
          EG_ERR_RETURN("bug encountered");
        };
        if (new_node_idx[p1] == p2) {
          delete_cell[cellId] = true;
        };
        if (new_node_idx[p2] == p1) {
          delete_cell[cellId] = true;
        };
      };
    };
  };
};

void vtkEgEliminateShortEdges::CopyPoints()
{
  node_mapping.fill(-1, m_Input->GetNumberOfPoints());
  vtkIdType newPointId = 0;
  for (vtkIdType pointId = 0; pointId < m_Input->GetNumberOfPoints(); ++pointId) {
    if(new_node_idx[pointId] < 0) {
      node_mapping[pointId] = newPointId;
      vec3_t x;
      m_Input->GetPoints()->GetPoint(pointId,x.data());
      m_Output->GetPoints()->SetPoint(newPointId,x.data());
      ++newPointId;
    };
  };
  for (vtkIdType pointId = 0; pointId < m_Input->GetNumberOfPoints(); ++pointId) {
    if(new_node_idx[pointId] >= 0) {
      if (node_mapping[new_node_idx[pointId]] < 0) {
        EG_ERR_RETURN("bug encountered");
      };
      node_mapping[pointId] = node_mapping[new_node_idx[pointId]];
    };
  };
};

void vtkEgEliminateShortEdges::CopyCells()
{
  for (vtkIdType cellId = 0; cellId < m_Input->GetNumberOfCells(); ++cellId) {
    if(!delete_cell[cellId]) {
      vtkIdType *old_pts;
      vtkIdType  Npts;
      m_Input->GetCellPoints(cellId, Npts, old_pts);
      vtkIdType *new_pts = new vtkIdType[Npts];
      for (int i = 0; i < Npts; ++i) {
        new_pts[i] = node_mapping[old_pts[i]];
      };
      vtkIdType newCellId = m_Output->InsertNextCell(m_Input->GetCellType(cellId), Npts, new_pts);
      copyCellData(m_Input, cellId, m_Output, newCellId);
      delete [] new_pts;
    };
  };
};

void vtkEgEliminateShortEdges::ExecuteEg()
{
  N_eliminated = 0;
  N_new_points = 0;
  N_new_cells  = 0;
  CheckEdges();
  CheckCells();
  int N = 0;
  for (vtkIdType pointId = 0; pointId < m_Input->GetNumberOfPoints(); ++pointId) {
    if(new_node_idx[pointId] < 0) {
      ++N_new_points;
    };
  };
  for (vtkIdType cellId = 0; cellId < m_Input->GetNumberOfCells(); ++cellId) {
    if(!delete_cell[cellId]) {
      ++N_new_cells;
    } else {
      ++N;
    };
  };
  allocateGrid(m_Output, N_new_cells, N_new_points);
  CopyPoints();
  CopyCells();
  UpdateCellIndex(m_Output);
  N_eliminated = N;
};
