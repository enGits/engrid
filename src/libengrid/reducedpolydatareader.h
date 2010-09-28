//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008-2010 enGits GmbH                                     +
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

#ifndef REDUCEDPOLYDATAREADER_H
#define REDUCEDPOLYDATAREADER_H

class ReducedPolyDataReader;

#include "iooperation.h"

/**
 * Reader for VTK legacy files containing vtkPolyData.
 * This reader will perform a reduction of the number of nodes and faces.
 */
class ReducedPolyDataReader : public IOOperation
{

private: // data-types

  struct Triangle
  {
    vtkIdType id_a, id_b, id_c;
    vec3_t a, b, c;
    vec3_t g1, g2, g3;
    mat3_t G, GI;
    double A;
    double smallest_length;
  };


private: // attributes

  double m_MaxEdgeLength;
  double m_MinEdgeLength;


protected: // methods

  void computeLevelSet(vtkUnstructuredGrid *grid, vtkPolyData* poly);
  virtual void operate();


public: // methods

  ReducedPolyDataReader();

};

#endif // REDUCEDPOLYDATAREADER_H
