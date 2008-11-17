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
#include "gridsmoother.h"
#include "computewalldistance.h"
#include "computegradg.h"
#include "smoothwalldistance.h"

// TESTING
#include <vtkXMLUnstructuredGridWriter.h>

GridSmoother::GridSmoother()
{
  N_iterations           = 20;
  N_relaxations          = 10;
  N_boundary_corrections = 20;
  layer_g                = 0.0;
  layer_dg               = 1.0;
  bg_grid                = vtkUnstructuredGrid::New();
  relax                  = 1.0;
  N_smooth_layers        = 0;
  ortho_factor           = 0.5;
};

bool GridSmoother::setNewPosition(vtkIdType id_node, vec3_t x_new)
{
  vec3_t x_old;
  grid->GetPoint(id_node, x_old.data());
  grid->GetPoints()->SetPoint(id_node, x_new.data());
  bool positive = true;
  foreach (int i_cells, n2c[id_node]) {
    vtkIdType id_cell = cells[i_cells];
    vtkIdType type_cell = grid->GetCellType(id_cell);
    if (isVolume(id_cell, grid)) {
      if (GeometryTools::cellVA(grid, id_cell) <= 0) {
        /*
        cout << "rejected move " << id_node << ',' << x_old << " -> " << x_new << endl;
        cout << id_cell << ',' << grid->GetCellType(id_cell) << endl << endl;
        */
        positive = false;
        grid->GetPoints()->SetPoint(id_node, x_old.data());
        break;
      };
    };
  };
  return positive;
};

void GridSmoother::setBackgroundGrid(vtkUnstructuredGrid *bg) 
{
  makeCopy(bg, bg_grid);
  gf.setNumSmoothLayers(N_smooth_layers);
  gf.setBoundaryCodes(boundary_codes);
  gf.setGrid(bg);
  //gf.setGrid(bg_grid);
};

/*
void GridSmoother::operate()
{
  EG_VTKDCC(vtkIntArray, bc, grid, "cell_code");
  for (int i_iterations = 0; i_iterations < N_iterations; ++i_iterations) {
    for (int i_nodes = 0; i_nodes < nodes.size(); ++i_nodes) {
      vec3_t x_old;
      vtkIdType id_node = nodes[i_nodes];
      grid->GetPoint(id_node, x_old.data());
      bool pure_tetra = true;
      bool pure_prism = true;
      vec3_t x_tetra(0,0,0);
      vec3_t x_prism(0,0,0);
      vec3_t x_tri  (0,0,0);
      int N_tetra = 0;
      int N_prism = 0;
      int N_tri   = 0;
      int N_quad  = 0;
      
      bool debug = false;
      
      foreach (vtkIdType i_cells, n2c[i_nodes]) {
        vtkIdType id_cell = cells[i_cells];
        vtkIdType type_cell = grid->GetCellType(id_cell);
        if ((type_cell == VTK_WEDGE) || (type_cell == VTK_HEXAHEDRON)) {
          pure_tetra = false;
          vtkIdType *pts, N_pts;
          grid->GetCellPoints(id_cell, N_pts, pts);
          for (int i_pts = 0; i_pts < N_pts; ++i_pts) {
            vec3_t x;
            grid->GetPoint(pts[i_pts], x.data());
          };
          N_pts = 3;
          if (type_cell == VTK_HEXAHEDRON) {
            N_pts = 4;
          };
          vector<vec3_t> x(N_pts);
          for (int i_pts = 0; i_pts < N_pts; ++i_pts) {
            grid->GetPoint(pts[i_pts], x[i_pts].data());
          };
          for (int i_pts = 0; i_pts < N_pts; ++i_pts) {
            if (pts[i_pts + N_pts] == id_node) {
              vec3_t n(0,0,0);
              double G;
              gf(x_old, G, n);
              n.normalise();
              double step = layer_g - G;
              x_prism += x_old + (relax*step)*n;
              ++N_prism;
            };
          };
        } else if (type_cell == VTK_TETRA) {
          pure_prism = false;
          x_tetra += cellCentre(grid, id_cell);
          ++N_tetra;
        } else if (type_cell == VTK_TRIANGLE) {
          x_tri += cellCentre(grid, id_cell);
          ++N_tri;
        } else if (type_cell == VTK_QUAD) {
          ++N_quad;
        };
      };
      vec3_t x_new;
      bool move = true;
      if ((N_tri > 0) && (N_quad == 0) && (N_prism == 0)) {
        x_new = 1.0/N_tri*x_tri;
      } else if (pure_tetra) {
        if (N_tetra == 0) EG_BUG;
        x_new = 1.0/N_tetra*x_tetra;
      } else {
        if (pure_prism) {
          x_new = x_old;
          move = false;
        } else {
          x_new = 1.0/N_prism*x_prism;
        };
      };
      if (move) {
        vec3_t dx = x_new - x_old;
        for (int i_boundary_correction = 0; i_boundary_correction < N_boundary_corrections; ++i_boundary_correction) {
          foreach (vtkIdType id_cell, n2c[i_nodes]) {
            vtkIdType type_cell = grid->GetCellType(id_cell);
            if (isSurface(id_cell, grid)) {
              if (type_cell == VTK_TRIANGLE) {
                vec3_t n = GeometryTools::cellNormal(grid, id_cell);
                n.normalise();
                dx -= (n*dx)*n;
              };
            };
          };
        };
        for (int i_relaxation = 0; i_relaxation < N_relaxations; ++i_relaxation) {
          if (debug) cout << i_relaxation << "  " << dx << endl;
          if (setNewPosition(id_node, x_old + dx)) {
            break;
          };
          dx *= 0.5;
        };
      };
    };
  };
};
*/

void GridSmoother::operate()
{
  EG_VTKDCC(vtkIntArray,    bc,          grid, "cell_code");
  EG_VTKDCN(vtkDoubleArray, node_g,      grid, "node_g");
  EG_VTKDCN(vtkDoubleArray, node_G,      grid, "node_G");
  EG_VTKDCN(vtkDoubleArray, node_vx,     grid, "node_vx");
  EG_VTKDCN(vtkDoubleArray, node_vy,     grid, "node_vy");
  EG_VTKDCN(vtkDoubleArray, node_vz,     grid, "node_vz");
  EG_VTKDCN(vtkIntArray,    node_status, grid, "node_status");
  EG_VTKDCN(vtkIntArray,    node_layer,  grid, "node_layer");
  
  for (int i_iterations = 0; i_iterations < N_iterations; ++i_iterations) {
    for (int i_nodes = 0; i_nodes < nodes.size(); ++i_nodes) {
      vec3_t x_old;
      vtkIdType id_node = nodes[i_nodes];
      grid->GetPoint(id_node, x_old.data());
      vec3_t x_tet(0,0,0);
      vec3_t x_pri(0,0,0);
      vec3_t x_tri(0,0,0);
      vec3_t x_qua(0,0,0);
      bool move  = true;
      int  N_tet = 0;
      int  N_pri = 0;
      int  N_tri = 0;
      int  N_qua = 0;
      bool f1 = false;
      bool f2 = false;
      bool f3 = false;
      vtkIdType id_lower_node = -1;
      foreach (int i_neigh_cells, n2c[i_nodes]) {
        vtkIdType id_neigh_cell = cells[i_neigh_cells];
        if (isSurface(grid, id_neigh_cell)) {
          if (boundary_codes.contains(bc->GetValue(id_neigh_cell))) {
            move = false;
            //break;
          };
          if (bc->GetValue(id_neigh_cell) == 1) f1 = true;
          if (bc->GetValue(id_neigh_cell) == 2) f2 = true;
          if (bc->GetValue(id_neigh_cell) == 3) f3 = true;
        };
        vtkIdType type_neigh_cell = grid->GetCellType(id_neigh_cell);
        vec3_t xc = cellCentre(grid, id_neigh_cell);
        if (type_neigh_cell == VTK_TRIANGLE) {
          x_tri += xc;
          ++N_tri;
        } else if (type_neigh_cell == VTK_QUAD) {
          x_qua += xc;
          ++N_qua;
        } else if (type_neigh_cell == VTK_TETRA) {
          x_tet += xc;
          ++N_tet;
        } else if (type_neigh_cell == VTK_WEDGE) {
          x_pri += xc;
          vtkIdType N_pts, *pts;
          grid->GetCellPoints(id_neigh_cell, N_pts, pts);
          for (int i_pts = 0; i_pts < 3; ++i_pts) {
            if (pts[i_pts + 3] == id_node) {
              if ((id_lower_node != -1) && (id_lower_node != pts[i_pts])) {
                EG_BUG;
              };
              id_lower_node = pts[i_pts];
            };
          };
          ++N_pri;
        };
      };
      vec3_t x_lower(0,0,0);
      if (id_lower_node != -1) {
        grid->GetPoint(id_lower_node, x_lower.data());
      };
      if (f1 && f2 && f3) {
        //cout << "found: " << id_node << ',' << x_old << endl;
      };
      if (node_status->GetValue(id_node) & 4) {
        move = false;
      };
      if (move) {
        vec3_t x1(0,0,0);
        
        if      (N_qua > 0) x1 = (1.0/N_qua)*x_qua;
        else if (N_tri > 0) x1 = (1.0/N_tri)*x_tri;
        else if (N_pri > 0) x1 = (1.0/N_pri)*x_pri;
        else if (N_tet > 0) x1 = (1.0/N_tet)*x_tet;
        else EG_BUG;
        
        if (node_layer->GetValue(id_node) > 0) {
          x1 = vec3_t(0,0,0);
          int N = 0;
          /*
          foreach (int i_neigh_nodes, n2n[i_nodes]) {
            vtkIdType id_neigh_node = nodes[i_neigh_nodes];
            if (node_layer->GetValue(id_node) == node_layer->GetValue(id_neigh_node)) {
              vec3_t xn;
              grid->GetPoint(id_neigh_node, xn.data());
              x1 += xn;
              ++N;
            };
          };
          */
          foreach (int i_neigh_cells, n2c[i_nodes]) {
            vtkIdType id_neigh_cell   = cells[i_neigh_cells];
            vtkIdType type_neigh_cell = grid->GetCellType(id_neigh_cell);
            if (type_neigh_cell == VTK_WEDGE) {
              vtkIdType N_pts, *pts;
              grid->GetCellPoints(id_neigh_cell, N_pts, pts);
              for (int i_pts = 0; i_pts < N_pts; ++i_pts) {
                vtkIdType id_neigh_node = pts[i_pts];
                if (node_layer->GetValue(id_node) == node_layer->GetValue(id_neigh_node)) {
                  vec3_t xn;
                  grid->GetPoint(id_neigh_node, xn.data());
                  x1 += xn;
                  ++N;
                };
              };
            };
          };
          if (N < 2) {
            cout << id_node << ',' << node_layer->GetValue(id_node) << ',' << N << ',' << x_old << ',' << n2n[i_nodes].size() << endl;
            //EG_BUG;
            x1 = x_old;
          } else {
            x1 *= 1.0/N;
          };
        };
        
        vec3_t Dx = x1 - x_old;
        double g = node_g->GetValue(id_node);
        if ((g > 0) && (id_lower_node != -1)) {
          /*
          x1[0] += node_vx->GetValue(id_node);
          x1[1] += node_vy->GetValue(id_node);
          x1[2] += node_vz->GetValue(id_node);
          */
          
          vec3_t n1(0,0,0);
          double G1;
          gf(x1, G1, n1);
          n1.normalise();
          vec3_t Dx1 = Dx;
          Dx1 -= (Dx1*n1)*n1;
          double step1 = g - G1;
          Dx1 += (relax*step1)*n1;
          
          /*
          vec3_t Dx1(0,0,0);
          {
            vec3_t n1(0,0,0);
            double G1;
            vec3_t n_lower;
            double G_lower;
            gf(x_lower, G_lower, n_lower);
            n_lower.normalise();
            double scale = max(0.0,(x1 - x_lower)*n_lower);
            vec3_t x_virt = x_lower + scale*n_lower;
            gf(x_virt, G1, n1);
            n1.normalise();
            double step1 = g - G1;
            Dx1 = (relax*step1)*n_lower;
          };
          */
          vec3_t n2(0,0,0);
          double G2;
          vec3_t n_lower;
          double G_lower;
          gf(x_lower, G_lower, n_lower);
          n_lower.normalise();
          double scale = max(0.0,(x_old - x_lower)*n_lower);
          vec3_t x_virt = x_lower + scale*n_lower;
          gf(x_virt, G2, n2);
          n2.normalise();
          double step2 = g - G2;
          vec3_t Dx2 = (relax*step2)*n_lower;
          double f1 = ortho_factor;
          if (node_status->GetValue(id_node) & 1) {
            f1 = 0.0;
          };
          x1 = f1*(x_old + Dx1) + (1-f1)*(x_virt + Dx2);
          
          Dx = x1 - x_old;
        };
        for (int i_boundary_correction = 0; i_boundary_correction < N_boundary_corrections; ++i_boundary_correction) {
          foreach (vtkIdType id_cell, n2c[i_nodes]) {
            if (isSurface(id_cell, grid)) {
              double A = GeometryTools::cellVA(grid, id_cell);
              if (A > 1e-20) {
                vec3_t n = GeometryTools::cellNormal(grid, id_cell);
                n.normalise();
                Dx -= (n*Dx)*n;
              };
            };
          };
        };
        for (int i_relaxation = 0; i_relaxation < N_relaxations; ++i_relaxation) {
          if (setNewPosition(id_node, x_old + Dx)) {
            break;
          };
          Dx *= 0.5;
        };
        {
          vec3_t x_new = x_old + Dx;
          vec3_t n;
          double G;
          gf(x_new, G, n);
          node_G->SetValue(id_node, G);
        };
        
      };
    };
  };
};

