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
#ifndef guieditboundaryconditions_H
#define guieditboundaryconditions_H

class GuiEditBoundaryConditions;

#include "boundarycondition.h"
#include "physicalboundaryconditions.h"
#include "dialogoperation.h"
#include "guivolumedelegate.h"
#include "filetemplate.h"
#include "multipagewidgetpage.h"
#include "volumedefinition.h"
#include "multipagewidget.h"

#include "ui_guieditboundaryconditions.h"

class GuiEditBoundaryConditions : public DialogOperation<Ui::GuiEditBoundaryConditions, Operation>
{
  
  Q_OBJECT;
  
private: // attributes
  
  // variables to store settings locally while changing them. They will be copied over to their GuiMainWindow counterparts once OK is clicked.
  QMap<int,BoundaryCondition> *m_BcMap;
  QMap<QString, VolumeDefinition> m_VolMap;
  QMap<QString,PhysicalBoundaryConditions> m_PhysicalBoundaryConditionsMap;
  
private: // utility attributes
  GuiVolumeDelegate *delegate;
  QVector <MultiPageWidgetPage*> m_page_vector;
  
  int m_PreviousSelected;
  QString m_PreviousSelectedName;
  int m_PreviousSelectedIndex;
  MultiPageWidget* m_multipagewidget_Solver;
  
  /// vector to hold the binaries
  QVector <QString> m_SolverBinary;
  
protected: // methods
  
  virtual void operate();
  void updateVol();
  void updatePhysicalBoundaryConditions();
  void setupSolvers();
  void saveSolverParameters();
  void saveMpiParameters();
  void loadMpiParameters();
  
protected slots:

  void addVol();
  void delVol();
  void addBoundaryType();
  void deleteBoundaryType();
  void changePhysicalValues();
  void loadPhysicalValues(QString name);
  void savePhysicalValues(QString name, int index);
  
public: // methods
  
  GuiEditBoundaryConditions();
  virtual ~GuiEditBoundaryConditions();

  virtual void before();
  void setMap(QMap<int,BoundaryCondition> *a_bcmap) { m_BcMap = a_bcmap; }
  
};

#endif
