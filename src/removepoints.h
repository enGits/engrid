//
// C++ Interface: removepoints
//
// Description: 
//
//
// Author: Mike Taverne <mtaverne@engits.com>, (C) 2009
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef REMOVEPOINTS_H
#define REMOVEPOINTS_H

#include "operation.h"

#include <vtkUnstructuredGrid.h>
#include <vtkPolyData.h>
#include <vtkCharArray.h>

#include <QSet>
#include <QVector>
#include <QString>
#include <QTextStream>
#include <QTime>

#include "egvtkobject.h"
#include "operation.h"
#include "vertexmeshdensity.h"
#include "smoothingutilities.h"
#include "swaptriangles.h"
#include "laplacesmoother.h"
#include "guimainwindow.h"

#include "geometrytools.h"
using namespace GeometryTools;

#include <cmath>
using namespace std;

#include <iostream>
using namespace std;

/**
	@author Mike Taverne <mtaverne@engits.com>
*/
class RemovePoints : public Operation
{
private:
  QMap <vtkIdType,bool> m_marked_cells;
  QMap <vtkIdType,bool> m_marked_nodes;
  
  QVector<vtkIdType> m_SelectedCells;
  QVector<vtkIdType> m_AllCells;
  QVector <vtkIdType> m_SelectedNodes;
  QVector <vtkIdType> m_AllNodes;
  QVector <int> m_hitlist;//Elements to be terminated (0=keep alive, 1=field agent to eliminate, 2=border agent to eliminate)
  QVector <int> m_offset;//offset caused by terminated elements
  
  int N_points;
  int N_cells;
  int N_newpoints;
  int N_newcells;
  int m_total_N_newpoints;
  int m_total_N_newcells;
  vtkIdType m_newNodeId;
  
  //attributes with setter functions
public:
  QSet<int> m_bcs;
  void SetBCS(QSet<int> a_bcs) {m_bcs=a_bcs;};
  bool remove_FP;
  void Set_remove_FP(bool B){remove_FP=B;};
  bool remove_EP;
  void Set_remove_EP(bool B){remove_EP=B;};
  
public:
  RemovePoints();
  
  void operate();
  
  
  int remove_FP_counter();
  int remove_EP_counter();
  int remove_EP_all();
  int remove_FP_all();
  
  ///Check if a field point needs to be removed
  bool remove_fieldpoint(vtkIdType P);
  ///Check if an edge point needs to be removed
  bool remove_edgepoint(vtkIdType P);
  
};

#endif
