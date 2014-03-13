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
// 
#ifndef TETGENOPERATION_H
#define TETGENOPERATION_H

#include "operation.h"
#include "tetgen.h"
#include "edgelengthsourcemanager.h"

class TetGenOperation;

class TetGenOperation : public Operation
{

protected: // attributes

  double m_MinimalEdgeLength;
  double m_MaximalEdgeLength;
  double m_GrowthFactor;
  double m_NodesPerQuarterCircle;
  double m_2dFeatureResolution;
  double m_3dFeatureResolution;
  int    m_OrgDir;
  int    m_CurDir;
  int    m_VolDir;

  EdgeLengthSourceManager m_ELSManager;

  QString m_TetGenPath;


protected: // methods

  void copyToTetGen(tetgenio &tgio, bool set_edge_lengths = false);
  void copyFromTetGen(tetgenio &tgio);

  //void copyToTetGen();
  //void copyFromTetGen();

  void tetgen(QString flags);
  void readSettings();


public:

  TetGenOperation();

};

#endif // TETGENOPERATION_H
