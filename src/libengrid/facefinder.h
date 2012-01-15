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
#ifndef FACEFINDER_H
#define FACEFINDER_H

#include "octree.h"
#include "triangle.h"
#include "timer.h"

#include <QVector>
#include <QList>

class FaceFinder : public EgVtkObject
{

  Octree m_Octree;
  vtkUnstructuredGrid *m_Grid;
  QVector<QList<vtkIdType> > m_Faces;
  double m_MinSize;
  int    m_MaxFaces;
  QVector<Triangle> m_Triangles;
  QVector<vec3_t>   m_Centres;
  QVector<double>   m_CritLength;
  Timer m_Timer;


private: // methods

  double calcCritLength(vtkIdType id_cell);
  int refine();


public: // methods

  FaceFinder();

  void setGrid(vtkUnstructuredGrid *grid);
  void setMaxNumFaces(int N) { m_MaxFaces = N; }
  void getCloseFaces(vec3_t x, QVector<vtkIdType> &faces);
  vtkIdType getClosestFace(vec3_t x, double &L);

};

#endif // FACEFINDER_H
