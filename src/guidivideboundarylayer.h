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
#ifndef guidivideboundarylayer_H
#define guidivideboundarylayer_H

class GuiDivideBoundaryLayer;

#include "dialogoperation.h"
#include "ui_guidivideboundarylayer.h"

class GuiDivideBoundaryLayer : public DialogOperation<Ui::GuiDivideBoundaryLayer, Operation>
{
  
private: // attributes
  
  int    N_layers;
  double h;
  double F;
  double f;
  double H;
  int    N_prisms;
  int    N_hexes;
  int    N_quads;
  bool   y_computed;
  
  QSet<QPair<vtkIdType,vtkIdType> > pairs;
  QVector<QVector<vtkIdType> > edges;
  QVector<bool> is_blayer_node;
  QVector<int> old2edge;
  QVector<double> x;
  QVector<double> y;
  QVector<bool> insert_cell;
  
  vtkUnstructuredGrid* m_rest_grid;///< used to store unselected volumes
  
private: // methods
  
  void findBoundaryLayer1();/// \todo Unused function. Delete if deprecated.
  bool findBoundaryLayer();
  void createEdges(vtkUnstructuredGrid *new_grid);
  inline double sech(double x) { return 1.0/cosh(x); };
  void bisectF(double &f1, double &f2);
  void computeF();
  void computeY();
  
protected: // methods
  
  virtual void before();
  virtual void operate1();/// \todo Unused function. Delete if deprecated.
  virtual void operate();
  void after();
};

#endif
