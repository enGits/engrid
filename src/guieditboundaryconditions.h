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
#include "physicalboundarycondition.h"
#include "dialogoperation.h"
#include "guivolumedelegate.h"
#include "filetemplate.h"
#include "multipagewidgetpage.h"
#include "volumedefinition.h"
#include "multipagewidget.h"

#include "ui_guieditboundaryconditions.h"

// tabs:
// boundary conditions -> add/del/update
// boundary types -> add/del/update + change/load/save
// solver ->setup/save
// MPI -> add/del/update

class GuiEditBoundaryConditions : public DialogOperation<Ui::GuiEditBoundaryConditions, Operation>
{

    Q_OBJECT;

  private: // attributes

    // variables to store settings locally while changing them. They will be copied over to their GuiMainWindow counterparts once OK is clicked.
    QMap<int, BoundaryCondition> *m_BcMap;
    QMap<QString, VolumeDefinition> m_VolMap;
    QMap<QString, PhysicalBoundaryCondition> m_PhysicalBoundaryConditionsMap;

  private: // utility attributes

    GuiVolumeDelegate *delegate;
    QVector <MultiPageWidgetPage*> m_page_vector;

    int m_PreviousSelected;
    QString m_PreviousSelectedName;
    int m_PreviousSelectedIndex;
    MultiPageWidget* m_multipagewidget_Solver;

    /// vector to hold the binaries
    QVector <QString> m_SolverBinary;

  public:
    GuiEditBoundaryConditions();
    virtual ~GuiEditBoundaryConditions();
    void setMap(QMap<int, BoundaryCondition> *a_bcmap) { m_BcMap = a_bcmap; }

  private:
    virtual void before();

  protected: // methods
    virtual void operate();

    // Boundary conditions tab
  protected:
    void updateVol();
  protected slots:
    void addVol();
    void delVol();

    // Boundary types tab
  protected:
    void updatePhysicalBoundaryConditions();
    void loadPhysicalValues(QString name);
    void savePhysicalValues(QString name, int index);
  protected slots:
    void addBoundaryType();
    void deleteBoundaryType();
    void changePhysicalValues();

    // Solver tab
  protected:
    void setupSolvers();
    void saveSolverParameters();
  protected slots:

    // MPI configuration tab
  protected:
    void loadMpiParameters();
    void saveMpiParameters();
    QString tableToString();
    void stringToTable(QString hostfile_txt);
  protected slots:
    void addProcess();
    void deleteProcess();
    void importHostFile();
    void exportHostFile();
};

#endif
