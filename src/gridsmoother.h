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

#include "operation.h"
#include "optimisation.h"

#include <vtkCellLocator.h>
#include <QSettings>

class GridSmoother : public Operation, public Optimisation
{
  
private: // attributes
  
  bool          m_SmoothPrisms;
  QVector<bool> m_NodeMarked;
  int           m_NumMarkedNodes;
  
protected: // attributes
  
  int    m_NumIterations;
  int    m_NumRelaxations;
  int    m_NumBoundaryCorrections;
  int    m_NumSearch;
  
  double m_LSearch;
  double m_FOld;
  double m_FNew;
  double m_FMaxOld;
  double m_FMaxNew;
  double m_ReductionFactor;
  
  double m_WTet;
  double m_WTetSave;
  double m_WH;
  double m_WPar;
  double m_WN;
  double m_WA;
  double m_WSkew;
  double m_WOrth;
  double m_WSharp1;
  double m_ESharp1;
  double m_WSharp2;
  double m_ESharp2;
  double m_H;

  bool m_StrictPrismChecking;

  QVector<vtkIdType>  m_FootToField;
  QVector<bool>       m_IsSharpNode;
  QVector<bool>       m_IsTripleNode;
  MeshPartition       m_BPart;

  double m_RelativeHeight;
  double m_CritAngle;
  
  bool m_SimpleOperation;

  struct stencil_node_t {
    vec3_t x;
    double C;
  };
  double m_V0;
  double m_L0;
  double m_SumC;
  int    m_INodesOpt;

  QList<stencil_node_t> m_Stencil;
  QVector<vtkIdType>    m_IdFoot;
  QVector<double>       m_L;
  QVector<vec3_t>       m_NodeNormal;

protected: // methods
  
  virtual void operate();
  virtual double func(vec3_t x);
  
  double errThickness(double x);
  
  bool setNewPosition(vtkIdType id_node, vec3_t x_new);
  void resetStencil();
  void addToStencil(double C, vec3_t x);
  void correctDx(int i_nodes, vec3_t &Dx);
  bool moveNode(int i_nodes, vec3_t &Dx);
  void markNodes();
  void setPrismWeighting() { m_WTetSave = m_WTet; m_WTet = 0; };
  void setAllWeighting() { m_WTet = m_WTetSave; };
  void computeNormals();
  void computeFeet();
  void simpleNodeMovement(int i_nodes);

  void operateOptimisation();
  void operateSimple();

public: // methods
  
  GridSmoother();
  void setNumIterations         (int N)    { m_NumIterations  = N; };
  void setNumRelaxations        (int N)    { m_NumRelaxations = N; };
  void setNumBoundaryCorrections(int N)    { m_NumBoundaryCorrections = N; };
  void setRelativeHeight        (double h) { m_RelativeHeight = h; }

  void prismsOn()  { m_SmoothPrisms = true; };
  void prismsOff() { m_SmoothPrisms = false; };
  void simpleOn()  { m_SimpleOperation = true; }
  void simpleOff() { m_SimpleOperation = false; }

  double improvement();
  
};


typedef GridSmoother SmoothVolumeGrid;

#endif
