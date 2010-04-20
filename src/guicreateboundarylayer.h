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
#ifndef guicreateboundarylayer_H
#define guicreateboundarylayer_H

class GuiCreateBoundaryLayer;

#include "dialogoperation.h"
#include "ui_guicreateboundarylayer.h"

#include <QProgressDialog>

class GuiCreateBoundaryLayer : public DialogOperation<Ui::GuiCreateBoundaryLayer, Operation>
{
  
  Q_OBJECT;
  
private: // attributes
  
  QVector<vtkIdType> layer_cells;

  int    m_NumIterations;
  int    m_NumPreSteps;
  int    m_NumPostSteps;
  bool   m_WriteDebugFile;
  double m_PostStrength;
  QSet<int> m_LayerAdjacentBoundaryCodes; /// Boundary codes of the surface we want to remove points on. Normally the one next to the prismatic boundary layer.

private: // methods
  
  void deleteTouchingPrisms(int layer, double L);
  void dump(vtkUnstructuredGrid *grid, QString name);
  
protected: // methods
  
  virtual void before();
  virtual void operate();
  
public: // methods
  
  GuiCreateBoundaryLayer();

private slots:
  void SelectAll_BC();
  void ClearAll_BC();
};

#endif
