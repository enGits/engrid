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

#ifndef GUICREATEHEXSHELL_H
#define GUICREATEHEXSHELL_H

#include "ui_guicreatehexshell.h"
#include "dialogoperation.h"

#include <QVector>

class GuiCreateHexShell : public DialogOperation<Ui::GuiCreateHexShell, Operation>
{

  Q_OBJECT

private: // attributes

  int m_NumICells;
  int m_NumJCells;
  int m_NumKCells;
  int m_NumINodes;
  int m_NumJNodes;
  int m_NumKNodes;
  int m_TotalNumberOfCells;
  int m_InnerBC;
  int m_OuterBC;

  double m_Dx, m_Dy, m_Dz;
  double m_X0, m_Y0, m_Z0;
  double m_H;
  QVector<vtkIdType> m_NodeIDs;
  QVector<vtkIdType> m_CellIDs;


protected: // methods

  virtual void before();
  virtual void operate();

  int indexCell(int i, int j, int k) { return i*m_NumJCells*m_NumKCells + j*m_NumKCells + k; }
  int indexNode(int i, int j, int k) { return i*m_NumJNodes*m_NumKNodes + j*m_NumKNodes + k; }
  void ijkCell(int idx, int& i, int& j, int& k);
  void ijkNode(int idx, int& i, int& j, int& k);
  void defineBoundaryCodes();
  void createGridWithNodes(vtkUnstructuredGrid *grid);
  void createHexCells(vtkUnstructuredGrid *grid);
  void createOuterBoundary(vtkUnstructuredGrid *grid);
  void createInnerBoundary(vtkUnstructuredGrid *grid);
  void createSourceBox();
  double x(int i) { return m_X0 + i*m_Dx; }
  double y(int j) { return m_Y0 + j*m_Dy; }
  double z(int k) { return m_Z0 + k*m_Dz; }



public: // methods

  GuiCreateHexShell();


public slots:

  void updateNumberOfCells(int = 0);

};

#endif // GUICREATEHEXSHELL_H
