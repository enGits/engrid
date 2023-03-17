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
#include "vtkEgExtractVolumeCells.h"

#include <vtkIdList.h>
#include <vtkSmartPointer.h>
#include <vtkType.h>

vtkStandardNewMacro(vtkEgExtractVolumeCells)

vtkEgExtractVolumeCells::vtkEgExtractVolumeCells()
{
  //SetClippingOff uses an IF to check the value first, which may have unexpected behavior in some OSes.
  m_Clip = false;   SetClippingOff();
  SetX(vec3_t(0,0,0));
  SetN(vec3_t(1,0,0));
  m_ExtrTetras   = true;
  m_ExtrPyramids = true;
  m_ExtrWedges   = true;
  m_ExtrHexes    = true;
  m_ExtrPolys    = true;
}

void vtkEgExtractVolumeCells::SetX(vec3_t x)
{
  m_X = x;
  Modified();
}

void vtkEgExtractVolumeCells::SetN(vec3_t n)
{
  m_N = n;
  Modified();
}

void vtkEgExtractVolumeCells::SetClippingOn()
{
  if (!m_Clip) {
    m_Clip = true;
    Modified();
  }
}

void vtkEgExtractVolumeCells::SetClippingOff()
{
  if (m_Clip) {
    m_Clip = false;
    Modified();
  }
}


void vtkEgExtractVolumeCells::SetAllOn()
{
  SetTetrasOn();
  SetPyramidsOn();
  SetWedgesOn();
  SetHexesOn();
  SetPolysOn();
}

void vtkEgExtractVolumeCells::SetAllOff()
{
  SetTetrasOff();
  SetPyramidsOff();
  SetWedgesOff();
  SetHexesOff();
  SetPolysOff();
}

void vtkEgExtractVolumeCells::SetTetrasOn()
{
  if (!m_ExtrTetras) {
    m_ExtrTetras = true;
    Modified();
  }
}

void vtkEgExtractVolumeCells::SetTetrasOff()
{
  if (m_ExtrTetras) {
    m_ExtrTetras = false;
    Modified();
  }
}

void vtkEgExtractVolumeCells::SetPyramidsOn()
{
  if (!m_ExtrPyramids) {
    m_ExtrPyramids = true;
    Modified();
  }
}

void vtkEgExtractVolumeCells::SetPyramidsOff()
{
  if (m_ExtrPyramids) {
    m_ExtrPyramids = false;
    Modified();
  }
}

void vtkEgExtractVolumeCells::SetWedgesOn()
{
  if (!m_ExtrWedges) {
    m_ExtrWedges = true;
    Modified();
  }
}

void vtkEgExtractVolumeCells::SetWedgesOff()
{
  if (m_ExtrWedges) {
    m_ExtrWedges = false;
    Modified();
  }
}

void vtkEgExtractVolumeCells::SetHexesOn()
{
  if (!m_ExtrHexes) {
    m_ExtrHexes = true;
    Modified();
  }
}

void vtkEgExtractVolumeCells::SetHexesOff()
{
  if (m_ExtrHexes) {
    m_ExtrHexes = false;
    Modified();
  }
}

void vtkEgExtractVolumeCells::SetPolysOn()
{
  if (!m_ExtrPolys) {
    m_ExtrPolys = true;
    Modified();
  }
}

void vtkEgExtractVolumeCells::SetPolysOff()
{
  if (m_ExtrPolys) {
    m_ExtrPolys = false;
    Modified();
  }
}

void vtkEgExtractVolumeCells::Setx(double x)
{
  m_X[0] = x;
  Modified();
}

void vtkEgExtractVolumeCells::Sety(double y)
{
  m_X[1] = y;
  Modified();
}

void vtkEgExtractVolumeCells::Setz(double z)
{
  m_X[2] = z;
  Modified();
}

void vtkEgExtractVolumeCells::Setnx(double nx)
{
  m_N[0] = nx;
  Modified();
}

void vtkEgExtractVolumeCells::Setny(double ny)
{
  m_N[1] = ny;
  Modified();
}

void vtkEgExtractVolumeCells::Setnz(double nz)
{
  m_N[2] = nz;
  Modified();
}


void vtkEgExtractVolumeCells::ExecuteEg()
{
  QSet<vtkIdType> ex_cells;
  for (vtkIdType id_cell = 0; id_cell < m_Input->GetNumberOfCells(); ++id_cell) {
    if (isVolume(id_cell, m_Input)) {
      bool select = true;
      vtkIdType type_cell = m_Input->GetCellType(id_cell);
      if (!m_ExtrTetras && type_cell == VTK_TETRA) {
        select = false;
      }
      if (!m_ExtrPyramids && type_cell == VTK_PYRAMID) {
        select = false;
      }
      if (!m_ExtrWedges && type_cell == VTK_WEDGE) {
        select = false;
      }
      if (!m_ExtrHexes && type_cell == VTK_HEXAHEDRON) {
        select = false;
      }
      if (!m_ExtrPolys && type_cell == VTK_POLYHEDRON) {
        select = false;
      }
      if (m_Clip && select) {
        QList<vtkIdType> pts;
        getPointsOfCell(m_Input, id_cell, pts);
        foreach (vtkIdType id_node, pts) {
          vec3_t x;
          m_Input->GetPoints()->GetPoint(id_node, x.data());
          if ((x - m_X)*m_N < 0) {
            select = false;
            break;
          }
        }
      }
      if (select) {
        ex_cells.insert(id_cell);
      }
    }
  }
  QVector<vtkIdType> cells(ex_cells.size());
  qCopy(ex_cells.begin(), ex_cells.end(), cells.begin());
  QVector<vtkIdType> nodes;
  QVector<int>       _nodes;
  getNodesFromCells(cells, nodes, m_Input);
  createNodeMapping(nodes, _nodes, m_Input);
  allocateGrid(m_Output, cells.size(), nodes.size());
  foreach(vtkIdType id_node, nodes) {
    vec3_t x;
    m_Input->GetPoints()->GetPoint(id_node, x.data());
    m_Output->GetPoints()->SetPoint(_nodes[id_node], x.data());
  }
  foreach(vtkIdType id_cell, cells) {
    vtkSmartPointer<vtkIdList> stream     = vtkSmartPointer<vtkIdList>::New();
    vtkSmartPointer<vtkIdList> new_stream = vtkSmartPointer<vtkIdList>::New();
    vtkIdType type_cell = m_Input->GetCellType(id_cell);
    m_Input->GetFaceStream(id_cell, stream);
    new_stream->SetNumberOfIds(stream->GetNumberOfIds());
    if (type_cell == VTK_POLYHEDRON) {
      vtkIdType num = stream->GetId(0);
      {
        new_stream->SetId(0, stream->GetId(0));
        vtkIdType id = 1;
        for (int i = 0; i < num; ++i) {
          vtkIdType num_pts = stream->GetId(id);
          new_stream->SetId(id, stream->GetId(id));
          ++id;
          for (int j = 0; j < num_pts; ++j) {
            new_stream->SetId(id, _nodes[stream->GetId(id)]);
            ++id;
          }
        }
      }
    } else {
      for (vtkIdType i = 0; i < stream->GetNumberOfIds(); ++i) {
        new_stream->SetId(i, _nodes[stream->GetId(i)]);
      }
    }

    vtkIdType id_new_cell = m_Output->InsertNextCell(type_cell, new_stream);
    copyCellData(m_Input, id_cell, m_Output, id_new_cell);
  }
}

