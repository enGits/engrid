// 
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008-2012 enGits GmbH                                     +
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
  
  int    m_NumLayers;
  int    m_NumPrisms;
  int    m_NumQuads;
  double m_RelativeHeight;
  double m_AbsoluteHeight;
  double m_Blending;
  double m_DesiredStretching;
  
  QSet<QPair<vtkIdType,vtkIdType> > m_Pairs;
  QVector<QVector<vtkIdType> >      m_Edges;
  QVector<bool>                     m_IsBlayerNode;
  QVector<int>                      m_Old2Edge;
  QVector<double>                   m_Y;
  QVector<bool>                     m_InsertCell;
  vtkUnstructuredGrid*              m_RestGrid;///< used to store unselected volumes
  
private: // methods
  
  bool findBoundaryLayer();
  void createEdges(vtkUnstructuredGrid *new_grid);
  void computeY();
  void finalise();

protected: // methods
  
  virtual void before();
  virtual void operate();

};

#endif
