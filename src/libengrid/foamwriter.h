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
#ifndef foamwriter_H
#define foamwriter_H

class FoamWriter;

#include "iooperation.h"
#include "polymesh.h"

/**
 * Writer for OpenFOAM poly-cell grids
 */
class FoamWriter : public IOOperation
{

protected: // attributes
  
  QString m_Path;
  QMap<int, QList<QString> > m_Bc2Vol;
  QString m_CurrentVolume;

protected: // methods
  
  void writePoints(const PolyMesh &poly);
  void writeFaces(const PolyMesh &poly);
  void writeOwner(const PolyMesh &poly);
  void writeNeighbour(const PolyMesh &poly);
  void writeBoundary(const PolyMesh &poly);

  bool    hasNeighbour(int bc);
  QString getNeighbourName(int bc);

  void writeSingleVolume();
  void writeMultipleVolumes();

  virtual void operate();
  
public: // methods
  
  FoamWriter();
  
};

#endif
