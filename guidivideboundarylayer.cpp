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
#include "guidivideboundarylayer.h"
#include "math/linsolve.h"

void GuiDivideBoundaryLayer::findBoundaryLayer()
{
  pairs.clear();
  N_prisms = 0;
  N_quads = 0;
  for (int i_cells = 0; i_cells < cells.size(); ++i_cells) {
    if (grid->GetCellType(cells[i_cells]) == VTK_WEDGE) {
      ++N_prisms;
      vtkIdType N_pts, *pts;
      grid->GetCellPoints(cells[i_cells],N_pts,pts);
      for (int j = 0; j < 3; ++j) {
        pairs.insert(QPair<vtkIdType,vtkIdType>(pts[j],pts[j+3]));
      };
      for (int j = 2; j < 5; ++j) {
        if (c2c[i_cells][j] != -1) {
          vtkIdType type_ncell = grid->GetCellType(cells[c2c[i_cells][j]]);
          if ((type_ncell != VTK_WEDGE) && (type_ncell != VTK_QUAD)) {
            EG_ERR_RETURN("unable to identify boundary layer");
          };
        } else {
          EG_BUG;
        };
      };
      for (int j = 0; j < 2; ++j) {
        if (c2c[i_cells][j] != -1) {
          vtkIdType type_ncell = grid->GetCellType(cells[c2c[i_cells][j]]);
          if (type_ncell == VTK_WEDGE) {
            EG_ERR_RETURN("the boundary layer seems to have been split already");
          };
        } else {
          EG_BUG;
        };
      };
    };
    if (grid->GetCellType(cells[i_cells]) == VTK_QUAD) {
      ++N_quads;
    };
  };
  if (N_prisms == 0) {
    EG_ERR_RETURN("unable to identify boundary layer");
  };
  is_blayer_node.clear();
  is_blayer_node.fill(false, grid->GetNumberOfPoints());
  
  QPair<vtkIdType,vtkIdType> P;
  foreach (P, pairs) {
    is_blayer_node[P.second] = true;
  };
};

void GuiDivideBoundaryLayer::bisectF(double &f1, double &f2)
{
  f = 0.5*(f1+f2);
  computeY();
  if (y[y.size()-1] > 1) f2 = f;
  else                   f1 = f;
  f = 0.5*(f1+f2);
};


//begin GROSSER MURKS 

void GuiDivideBoundaryLayer::computeF()
{
  double f1 = 0;
  double f2 = 100;//1.0/h;
  double err = 0;
  y[1] = min(0.33,y[1]);
  do {
    //cout << f1 << ',' << f2 << endl;
    bisectF(f1,f2);
    //while (fabs(1-y[y.size()-1]) > 1e-6);
    err = fabs((y[y.size()-1]-1)/(y[y.size()-1]-y[y.size()-2]));
  } while (err > 0.01);
  f = 0.5*(f1+f2);
};

void GuiDivideBoundaryLayer::computeY()
{
  double C = F;
  for (int i = 2; i < y.size(); ++i) {
    y[i] = y[i-1] + C*(y[i-1]-y[i-2]);
    C *= f;
  };
};

//end GROSSER MURKS

void GuiDivideBoundaryLayer::createEdges(vtkUnstructuredGrid *new_grid)
{
  edges.fill(QVector<vtkIdType>(N_layers+1), pairs.size());
  old2edge.fill(-1, grid->GetNumberOfPoints());
  int N = 0;
  vtkIdType id_new_node = grid->GetNumberOfPoints();
  QPair<vtkIdType,vtkIdType> P;
  double max_step = 0;
  double ymax     = 0;
  double ymin     = 1e99;
  foreach (P, pairs) {
    edges[N][0] = P.first;
    edges[N][N_layers] = P.second;
    old2edge[P.first] = N;
    old2edge[P.second] = N;
    
    vec3_t x1,x2;
    grid->GetPoint(P.first, x1.data());
    grid->GetPoint(P.second, x2.data());
    vec3_t n = x2-x1;
    if (!ui.checkBoxH->isChecked() || !y_computed) {
      y.resize(N_layers + 1);
      x.resize(y.size());
      for (int i = 0; i < x.size(); ++i) {
        x[i] = i*1.0/(x.size() - 1);
      };
      double Dy1;
      if (ui.checkBoxH->isChecked()) {
        Dy1 = h;
      } else {
        Dy1 = h/n.abs();
      };
      y[0]   = 0;
      y[1]   = Dy1;
      computeF();
      computeY();
      y_computed = true;
    };
    ymin = min(ymin,y[1]*n.abs());
    ymax = max(ymax,y[1]*n.abs());
    for (int i = 1; i < N_layers; ++i) {
      vec3_t x = x1 + y[i]*n;
      max_step = max(max_step,(y[i+1]-y[i])/(y[i]-y[i-1]));
      new_grid->GetPoints()->SetPoint(id_new_node, x.data());
      edges[N][i] = id_new_node;
      ++id_new_node;
    };
    ++N;
  };
  cout << LINE;
  cout << "maximal increment : " << max_step << endl;
  cout << "min(y) : " << ymin << endl;
  cout << "max(y) : " << ymax << endl;
  cout << LINE;
};

void GuiDivideBoundaryLayer::operate()
{
  y_computed = false;
  N_layers = ui.spinBoxLayers->value();
  h = ui.lineEditH->text().toDouble();
  F = ui.doubleSpinBoxF->value();
  cout << "dividing boundary layer into " << N_layers << " layers:" << endl;
  findBoundaryLayer();
  
  EG_VTKSP(vtkUnstructuredGrid,new_grid);
  allocateGrid(new_grid, grid->GetNumberOfCells() + (N_prisms + N_quads)*(N_layers-1), grid->GetNumberOfPoints() + pairs.size()*(N_layers-1));
  
  // copy existing mesh without prisms and adjacent cells
  vtkIdType id_new_node = 0;
  for (vtkIdType id_node = 0; id_node < grid->GetNumberOfPoints(); ++id_node) {
    vec3_t x;
    grid->GetPoint(id_node, x.data());
    new_grid->GetPoints()->SetPoint(id_new_node, x.data());
    copyNodeData(grid, id_node, new_grid, id_new_node);
    ++id_new_node;
  };
  vtkIdType id_new_cell;
  for (vtkIdType id_cell = 0; id_cell < grid->GetNumberOfCells(); ++id_cell) {
    vtkIdType N_pts, *pts;
    grid->GetCellPoints(id_cell, N_pts, pts);
    bool insert_cell = true;
    if (grid->GetCellType(id_cell) == VTK_WEDGE) insert_cell = false;
    if (grid->GetCellType(id_cell) == VTK_QUAD) insert_cell = false;
    if (insert_cell) {
      id_new_cell = new_grid->InsertNextCell(grid->GetCellType(id_cell), N_pts, pts);
      copyCellData(grid, id_cell, new_grid, id_new_cell);
    };
  };
  
  // create divided boundary layer
  createEdges(new_grid);
  
  for (vtkIdType id_cell = 0; id_cell < grid->GetNumberOfCells(); ++id_cell) {
    if (grid->GetCellType(id_cell) == VTK_WEDGE) {
      vtkIdType N_pts, *pts;
      grid->GetCellPoints(id_cell, N_pts, pts);
      for (int i = 0; i < N_layers; ++i) {
        vtkIdType p[6];
        p[0] = edges[old2edge[pts[0]]][i];
        p[1] = edges[old2edge[pts[1]]][i];
        p[2] = edges[old2edge[pts[2]]][i];
        p[3] = edges[old2edge[pts[0]]][i+1];
        p[4] = edges[old2edge[pts[1]]][i+1];
        p[5] = edges[old2edge[pts[2]]][i+1];
        id_new_cell = new_grid->InsertNextCell(VTK_WEDGE, 6, p);
        copyCellData(grid, id_cell, new_grid, id_new_cell);
      };
    };
    if (grid->GetCellType(id_cell) == VTK_QUAD) {
      vtkIdType N_pts, *pts;
      grid->GetCellPoints(id_cell, N_pts, pts);
      if ((old2edge[pts[0]] != -1) && (old2edge[pts[1]] != -1) && (old2edge[pts[2]] != -1) && (old2edge[pts[3]] != -1)) {
        for (int i = 0; i < N_layers; ++i) {
          vtkIdType p[4];
          p[0] = edges[old2edge[pts[0]]][i];
          p[1] = edges[old2edge[pts[1]]][i];
          p[2] = edges[old2edge[pts[1]]][i+1];
          p[3] = edges[old2edge[pts[0]]][i+1];
          id_new_cell = new_grid->InsertNextCell(VTK_QUAD, 4, p);
          copyCellData(grid, id_cell, new_grid, id_new_cell);
        };
      };
    };
  };
  
  makeCopy(new_grid, grid);
};

