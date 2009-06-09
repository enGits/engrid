//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008,2009 Oliver Gloth                                     +
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
#include "gmshreader.h"
#include "correctsurfaceorientation.h"

void GmshReader::readAscii1(vtkUnstructuredGrid *grid)
{
  vtkIdType Nnodes, Ncells;
  QFile file(getFileName());
  file.open(QIODevice::ReadOnly | QIODevice::Text);
  QTextStream f(&file);
  QString word;
  f >> word;
  if (word != "$NOD") EG_ERR_RETURN("$NOD expected");
  f >> Nnodes;
  EG_VTKSP(vtkUnstructuredGrid, ug);
  QVector<vtkIdType> idxmap(Nnodes + 1);
  QVector<vec3_t> x(Nnodes);
  for (vtkIdType i = 0; i < Nnodes; ++i) {
    int ir;
    f >> ir >> x[i][0] >> x[i][1] >> x[i][2];
    idxmap[ir] = i;
  }
  f >> word;
  if (word != "$ENDNOD") EG_ERR_RETURN("$ENDNOD expected");
  f >> word;
  if (word != "$ELM") EG_ERR_RETURN("$ELM expected");
  f >> Ncells;
  allocateGrid(ug, Ncells, Nnodes, false);
  for (vtkIdType i = 0; i < Nnodes; ++i) {
    ug->GetPoints()->SetPoint(i, x[i].data());
  }
  EG_VTKSP(vtkIntArray, cell_code);
  cell_code->SetName("cell_code");
  cell_code->SetNumberOfValues(Ncells);
  ug->GetCellData()->AddArray(cell_code);
  for (vtkIdType i = 0; i < Ncells; ++i) {
    f >> word;
    int elm_type, reg_phys;
    f >> elm_type;
    f >> reg_phys;
    f >> word;
    f >> word;
    cell_code->SetValue(i, reg_phys);
    vtkIdType pts[8];
    size_t node;
    if        (elm_type == 2) { // triangle
      for (vtkIdType j = 0; j < 3; ++j) { 
        f >> node; 
        pts[j] = idxmap[node]; 
      }
      ug->InsertNextCell(VTK_TRIANGLE,3,pts);
    } else if (elm_type == 3) { // quad
      for (vtkIdType j = 0; j < 4; ++j) { 
        f >> node; 
        pts[j] = idxmap[node]; 
      }
      ug->InsertNextCell(VTK_QUAD,4,pts);
    } else if (elm_type == 4) { // tetrahedron
      for (vtkIdType j = 0; j < 4; ++j) { 
        f >> node; 
        pts[j] = idxmap[node]; 
      }
      vtkIdType h = pts[0];
      pts[0] = pts[1];
      pts[1] = h;
      ug->InsertNextCell(VTK_TETRA,4,pts);
    } else if (elm_type == 5) { // hexhedron
      for (vtkIdType j = 0; j < 8; ++j) { 
        f >> node; 
        pts[j] = idxmap[node]; 
      }
      ug->InsertNextCell(VTK_HEXAHEDRON,8,pts);
    } else if (elm_type == 6) { // prism/wedge
      for (vtkIdType j = 0; j < 3; ++j) { 
        f >> node; 
        pts[j+3] = idxmap[node]; 
      }
      for (vtkIdType j = 0; j < 3; ++j) { 
        f >> node; 
        pts[j] = idxmap[node]; 
      }
      ug->InsertNextCell(VTK_WEDGE,6,pts);
    } else if (elm_type == 7) { // pyramid
      for (vtkIdType j = 0; j < 5; ++j) { 
        f >> node; 
        pts[j] = idxmap[node]; 
      }
      ug->InsertNextCell(VTK_PYRAMID,5,pts);
    }
  }
  ug->GetCellData()->AddArray(cell_code);
  grid->DeepCopy(ug);
}

void GmshReader::readAscii2(vtkUnstructuredGrid *grid)
{
  vtkIdType Nnodes, Ncells;
  QFile file(getFileName());
  file.open(QIODevice::ReadOnly | QIODevice::Text);
  QTextStream f(&file);
  QString word;
  f >> word;
  if (word != "$MeshFormat") EG_ERR_RETURN("$MeshFormat expected");
  f >> word;
  f >> word;
  f >> word;
  f >> word;
  if (word != "$EndMeshFormat") EG_ERR_RETURN("$EndMeshFormat expected");
  f >> word;
  if (word != "$Nodes") EG_ERR_RETURN("$Nodes expected");
  f >> Nnodes;
  EG_VTKSP(vtkUnstructuredGrid, ug);
  QVector<vtkIdType> idxmap(Nnodes + 1);
  QVector<vec3_t> x(Nnodes);
  for (vtkIdType i = 0; i < Nnodes; ++i) {
    int ir;
    f >> ir >> x[i][0] >> x[i][1] >> x[i][2];
    idxmap[ir] = i;
  }
  f >> word;
  if (word != "$EndNodes") EG_ERR_RETURN("$EndNotes expected");
  f >> word;
  if (word != "$Elements") EG_ERR_RETURN("$Elements expected");
  f >> Ncells;
  allocateGrid(ug, Ncells, Nnodes, false);
  for (vtkIdType i = 0; i < Nnodes; ++i) {
    ug->GetPoints()->SetPoint(i, x[i].data());
  }
  EG_VTKSP(vtkIntArray, cell_code);
  cell_code->SetName("cell_code");
  cell_code->SetNumberOfValues(Ncells);
  ug->GetCellData()->AddArray(cell_code);

  for (vtkIdType i = 0; i < Ncells; ++i) {
    f >> word;
    int elm_type, Ntags, bc;
    f >> elm_type;
    f >> Ntags;
    if (Ntags == 0) {
      bc = 1;
    } else {
      f >> bc;
      if (bc <= 0) {
        bc = 99;
      }
      for (int j = 1; j < Ntags; ++j) f >> word;
    }
    cell_code->SetValue(i, bc);
    vtkIdType pts[8];
    size_t node;
    if        (elm_type == 2) { // triangle
      for (vtkIdType j = 0; j < 3; ++j) { 
        f >> node; 
        pts[j] = idxmap[node]; 
      }
      ug->InsertNextCell(VTK_TRIANGLE,3,pts);
    } else if (elm_type == 3) { // quad
      for (vtkIdType j = 0; j < 4; ++j) { 
        f >> node; 
        pts[j] = idxmap[node]; 
      }
      ug->InsertNextCell(VTK_QUAD,4,pts);
    } else if (elm_type == 4) { // tetrahedron
      for (vtkIdType j = 0; j < 4; ++j) { 
        f >> node; 
        pts[j] = idxmap[node]; 
      }
      ug->InsertNextCell(VTK_TETRA,4,pts);
    } else if (elm_type == 5) { // hexhedron
      for (vtkIdType j = 0; j < 8; ++j) { 
        f >> node; 
        pts[j] = idxmap[node]; 
      }
      ug->InsertNextCell(VTK_HEXAHEDRON,8,pts);
    } else if (elm_type == 6) { // prism/wedge
      for (vtkIdType j = 0; j < 3; ++j) { 
        f >> node; 
        pts[j+3] = idxmap[node]; 
      }
      for (vtkIdType j = 0; j < 3; ++j) { 
        f >> node; 
        pts[j] = idxmap[node]; 
      }
      ug->InsertNextCell(VTK_WEDGE,6,pts);
    } else if (elm_type == 7) { // pyramid
      for (vtkIdType j = 0; j < 5; ++j) { 
        f >> node; 
        pts[j] = idxmap[node]; 
      }
      ug->InsertNextCell(VTK_PYRAMID,5,pts);
    }
  }
  ug->GetCellData()->AddArray(cell_code);
  grid->DeepCopy(ug);
}

void GmshReader::operate()
{
  try {
    readInputFileName();
    if (isValid()) {
      if (format == ascii1) {
        readAscii1(grid);
      }
      if (format == ascii2) {
        readAscii2(grid);
      }
      createBasicFields(grid, grid->GetNumberOfCells(), grid->GetNumberOfPoints());
      UpdateCellIndex(grid);
      CorrectSurfaceOrientation corr_surf;
      corr_surf.setGrid(grid);
      corr_surf();
    }
  } catch (Error err) {
    err.display();
  }
}
