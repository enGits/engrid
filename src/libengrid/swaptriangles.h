// 
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008-2013 enGits GmbH                                      +
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
#ifndef swaptriangles_H
#define swaptriangles_H

class SwapTriangles;

#include "surfaceoperation.h"
#include "cadinterface.h"

/**
  * \todo This class desperately needs a clean-up and optimisation!
  */
class SwapTriangles : public SurfaceOperation
{
  
private: // attributes
  
  QVector<bool> m_Swapped;
  bool          m_RespectBC;
  bool          m_FeatureSwap;
  bool          m_SmallAreaSwap;
  bool          m_Verbose;
  int           m_MaxNumLoops;
  double        m_SmallAreaRatio;
  double        m_SurfErrorThreshold;
  double        m_SurfErrorRatio;
  double        m_AverageSurfaceError;
  double        m_SurfaceErrorDeviation;
  double        m_DelaunayThreshold;
  int           m_NumSwaps;

private: // methods
  
  ///returns true if performing a swap on the stencil does not change the orientation of the cells (tetra volume test)
  bool testOrientation(stencil_t S);
  
  ///returns true if id_node1 is linked to id_node2
  bool isEdge(vtkIdType id_node1, vtkIdType id_node2);

protected: // methods
  
  int swap();
  double computeSurfaceDistance(vec3_t x1, vec3_t x2, vec3_t x3, CadInterface* cad_interface);
  double edgeAngle(vec3_t x1, vec3_t x2, vec3_t x3, vec3_t x4);
  vtkIdType neighbourNode(vtkIdType id_node0, vtkIdType id_node1, vtkIdType id_node2);
  bool swapDueToSurfaceNoise(stencil_t S);
  void computeSurfaceErrors(const QVector<vec3_t> &x, int bc, double &err1, double &err2);
  void computeAverageSurfaceError();
  virtual void operate();

public:

  SwapTriangles();

  void setRespectBC(bool b)      { m_RespectBC   = b; }
  void setFeatureSwap(bool b)    { m_FeatureSwap = b; }
  void setMaxNumLoops(int n)     { m_MaxNumLoops = n; }
  void setSmallAreaSwap(bool b)  { m_SmallAreaSwap = b; }
  void setVerboseOn() { m_Verbose = true; }
  void setVerboseOff() { m_Verbose = false; }
  void setDelaunayThreshold(double t) { m_DelaunayThreshold = t; }
  int  getNumSwaps() { return m_NumSwaps; }

};

#endif
