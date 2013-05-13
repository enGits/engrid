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

#ifndef gridsmoother_H
#define gridsmoother_H

class GridSmoother;

#include "surfaceoperation.h"

#include <vtkCellLocator.h>
#include <QSettings>

class GridSmoother : public SurfaceOperation
{

private: // types

  struct rule_t {
    QSet<int> bcs;
    double h;
  };
  
private: // attributes
  
  QVector<bool>   m_NodeMarked;
  QVector<bool>   m_SurfNode;
  int             m_NumMarkedNodes;
  QList<rule_t>   m_Rules;
  bool            m_FirstCall;
  QVector<vec3_t> m_GeoNormal;
  
protected: // attributes
  
  int m_NumIterations;
  int m_NumRelaxations;
  int m_NumBoundaryCorrections;
  int m_NumSearch;
  int m_NumNormalRelaxations;
  int m_NumHeightRelaxations;

  double m_Blending;
  double m_AbsoluteHeight;
  double m_RelativeHeight;
  double m_CritAngle;
  double m_RadarAngle;
  double m_MaxHeightInGaps;
  double m_DesiredStretching;
  double m_FarRatio;
  int    m_NumLayers;

  bool m_StrictPrismChecking;
  bool m_CollisionDetected;

  QVector<vtkIdType> m_FootToField;

  QVector<vtkIdType> m_IdFoot;
  QVector<double>    m_Height;
  QSet<int>          m_LayerAdjacentBoundaryCodes;

protected: // methods
  
  virtual void operate();
  
  bool setNewPosition(vtkIdType id_node, vec3_t x_new);
  void correctDx(int i_nodes, vec3_t &Dx);
  bool moveNode(int i_nodes, vec3_t &Dx);
  void markNodes();
  void computeNormals();
  void relaxNormalVectors();
  void correctNormalVectors();
  void computeDesiredHeights();
  void computeHeights();
  void computeFeet();
  void simpleNodeMovement(int i_nodes);
  void getRules();

  void writeDebugFile(QString file_name);

public: // methods
  
  GridSmoother();
  void setNumIterations         (int N)    { m_NumIterations  = N; }
  void setNumRelaxations        (int N)    { m_NumRelaxations = N; }
  void setNumHeightRelaxations  (int N)    { m_NumHeightRelaxations = N; }
  void setNumNormalRelaxations  (int N)    { m_NumNormalRelaxations = N; }
  void setNumBoundaryCorrections(int N)    { m_NumBoundaryCorrections = N; }
  void setRelativeHeight        (double h) { m_RelativeHeight = h; }
  void setAbsoluteHeight        (double h) { m_AbsoluteHeight = h; }
  void setBlending              (double b) { m_Blending = b; }
  void setDesiredStretching     (double s) { m_DesiredStretching = s; }
  void setFarRatio              (double r) { m_FarRatio = r; }
  void forceNormalCalculation   ()         { m_FirstCall = true; }

  void setLayerAdjacentBoundaryCodes(const QSet<int> &abcs) { m_LayerAdjacentBoundaryCodes = abcs; }

};


typedef GridSmoother SmoothVolumeGrid;

#endif
