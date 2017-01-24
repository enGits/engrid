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

#ifndef DRNUMWRITER_H
#define DRNUMWRITER_H

#include "iooperation.h"
#include "solverobject.h"

#include <QFile>
#include <QTextStream>
#include <QMap>

/**
 * @brief A very experimental export function for DrNUM grids.
 */
class DrNumWriter : public IOOperation, public SolverObject
{

protected: // data types


protected: // attributes

  double m_MaximalEdgeLength;
  double m_MinimalEdgeLength;
  double m_GrowthFactor;
  vec3_t m_X1;
  vec3_t m_X2;

  vtkSmartPointer<vtkUnstructuredGrid> m_BackupGrid;

protected: // methods

  QList<BoundaryCondition> getBcsOfType(QString type);
  void                     prepareLevelSets(QList<BoundaryCondition> bc, double distance);
  void                     prepareWallLevelSets(QList<BoundaryCondition> bc);
  void                     readSettings();
  double                   edgeLength(QString bc_name);
  void                     computeBBox();
  void                     writeGlobals();
  void                     backup();
  void                     restore();

  virtual void operate();


public: // methods

  DrNumWriter();

};

#endif // DRNUMWRITER_H
