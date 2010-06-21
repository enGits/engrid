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
#ifndef LAPLACESMOOTHER_H
#define LAPLACESMOOTHER_H

#include "surfaceoperation.h"
#include "surfaceprojection.h"

class LaplaceSmoother : public SurfaceOperation
{

private:

  QSet<int> m_BCs;
  int       m_NumberOfIterations;
  bool      m_UseProjection;
  bool      m_UseNormalCorrection;
  double    m_UnderRelaxation;
  bool      m_Success;
  int       m_ProjectionIterations;
  bool      m_FreeProjectionForEdges;

  QVector<QVector<int> > m_NodeToBc;

  bool      m_correctCurvature;
  bool      m_NoCheck;
  
private: // methods

  bool setNewPosition(vtkIdType id_node, vec3_t x_new);
  bool moveNode(vtkIdType id_node, vec3_t &Dx);


public:

  LaplaceSmoother(); ///< default constructor
  virtual void operate(); ///< Run operation
  void setNumberOfIterations(int N) { m_NumberOfIterations = N;} ///< Set number of iterations
  void setProjectionOn() { m_UseProjection = true; }
  void setProjectionOff() { m_UseProjection = false; }
  void setNormalCorrectionOn() { m_UseNormalCorrection = true; }
  void setNormalCorrectionOff() { m_UseNormalCorrection = false; }
  bool succeeded() { return m_Success; }

public:

  void setCorrectCurvature(bool b) { m_correctCurvature = b; }
  bool getCorrectCurvature() { return m_correctCurvature; }
  void setNoCheck(bool b) { m_NoCheck = b; }
  bool getNoCheck() { return m_NoCheck; }
  void setProjectionIterations(int n) { m_ProjectionIterations = n; }
  void setFreeProjectionForEdgesOn() { m_FreeProjectionForEdges = true; }
  void setFreeProjectionForEdgesOff() { m_FreeProjectionForEdges = false; }

};

#endif
