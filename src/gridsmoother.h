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

#ifndef gridsmoother_H
#define gridsmoother_H

class GridSmoother;

#include "surfaceoperation.h"

#include <vtkCellLocator.h>
#include <QSettings>

class GridSmoother : public SurfaceOperation
{
  
private: // attributes
  
  //bool          m_SmoothPrisms;
  QVector<bool> m_NodeMarked;
  QVector<bool> m_SurfNode;
  int           m_NumMarkedNodes;
  
protected: // attributes
  
  int m_NumIterations;
  int m_NumRelaxations;
  int m_NumBoundaryCorrections;
  int m_NumSearch;
  int m_NumNormalRelaxations;
  int m_NumHeightRelaxations;

  //double m_LSearch;
  //double m_FOld;
  //double m_FNew;
  //double m_FMaxOld;
  //double m_FMaxNew;
  //double m_ReductionFactor;
  double m_Blending;
  double m_AbsoluteHeight;
  double m_RelativeHeight;
  double m_CritAngle;

  bool m_StrictPrismChecking;

  QVector<vtkIdType> m_FootToField;
  //QVector<bool>      m_IsSharpNode;
  //MeshPartition      m_BPart;

  //double m_V0;
  //double m_L0;
  //double m_SumC;
  //int    m_INodesOpt;

  QVector<vtkIdType> m_IdFoot;
  //QVector<double>    m_L;
  QVector<double>    m_Height;
  QVector<vec3_t>    m_NodeNormal;

protected: // methods
  
  virtual void operate();
  
  bool setNewPosition(vtkIdType id_node, vec3_t x_new);
  void correctDx(int i_nodes, vec3_t &Dx);
  bool moveNode(int i_nodes, vec3_t &Dx);
  void markNodes();
  void computeNormals();
  void relaxNormalVectors();
  void correctNormalVectors();
  void computeHeights();
  void computeFeet();
  void simpleNodeMovement(int i_nodes);
  bool noCollision(vtkIdType id_node);

  void writeDebugFile(QString file_name);

public: // methods
  
  GridSmoother();
  void setNumIterations         (int N)    { m_NumIterations  = N; };
  void setNumRelaxations        (int N)    { m_NumRelaxations = N; };
  void setNumBoundaryCorrections(int N)    { m_NumBoundaryCorrections = N; };
  void setRelativeHeight        (double h) { m_RelativeHeight = h; }
  void setAbsoluteHeight        (double h) { m_AbsoluteHeight = h; }
  void setBlending              (double b) { m_Blending = b; }

  //void prismsOn()  { m_SmoothPrisms = true; };
  //void prismsOff() { m_SmoothPrisms = false; };
  
};


typedef GridSmoother SmoothVolumeGrid;

#endif
