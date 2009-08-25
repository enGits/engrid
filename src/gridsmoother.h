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

#ifndef gridsmoother_H
#define gridsmoother_H

class GridSmoother;

#include "operation.h"
#include "optimisation.h"

#include <vtkCellLocator.h>
#include <QSettings>

class GridSmoother : public Operation, public Optimisation
{
  
private: // attributes
  
  bool smooth_prisms;
  QVector<bool> node_marked;
  int N_marked_nodes;
  bool dbg;
  
protected: // attributes
  
  int    N_iterations;
  int    N_relaxations;
  int    N_boundary_corrections;
  int    N_search;
  
  double L_search;
  double F_old;
  double F_new;
  double F_max_old;
  double F_max_new;
  
  double w_tet;
  double w_tet_save;
  double w_h;
  double e_h;
  double w_par;
  double w_n;
  double w_A;
  double w_skew;
  double w_orth;
  double w_sharp1;
  double e_sharp1;
  double w_sharp2;
  double e_sharp2;
  double H;

  double under_relaxation;
  
  struct stencil_node_t {
    vec3_t x;
    double C;
  };
  double V0;
  double L0;
  double sum_C;
  int i_nodes_opt;
  QList<stencil_node_t> stencil;

  QVector<vtkIdType> m_IdFoot;
  QVector<double> m_L;
  
protected: // methods
  
  virtual void operate();
  virtual double func(vec3_t x);
  
  double errThickness(double x);
  
  bool setNewPosition(vtkIdType id_node, vec3_t x_new);
  void resetStencil();
  void addToStencil(double C, vec3_t x);
  void correctDx(int i_nodes, vec3_t &Dx);
  bool moveNode(int i_nodes, vec3_t &Dx);
  void markNodes();
  void setPrismWeighting() { w_tet_save = w_tet; w_tet = 0; };
  void setAllWeighting() { w_tet = w_tet_save; };
    
public: // methods
  
  GridSmoother();
  void setNumIterations         (int N) { N_iterations  = N; };
  void setNumRelaxations        (int N) { N_relaxations = N; };
  void setNumBoundaryCorrections(int N) { N_boundary_corrections = N; };
  
  void prismsOn() { smooth_prisms = true; };
  void prismsOff() { smooth_prisms = false; };
  
  double improvement();
  double lastTotalError() { return F_new; }
  
};


typedef GridSmoother SmoothVolumeGrid;

#endif
