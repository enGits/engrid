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
#ifndef REMOVEPOINTS_H
#define REMOVEPOINTS_H

#include "surfaceoperation.h"

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

class RemovePoints : public SurfaceOperation
{
  private:
    int m_NumRemoved;

    //attributes with setter functions
  public:
    QSet<int> m_bcs;
    void setBCS( QSet<int> a_bcs ) { m_bcs = a_bcs;}
    bool remove_FP;
    void set_remove_FP( bool B ) { remove_FP = B;}
    bool remove_EP;
    void set_remove_EP( bool B ) { remove_EP = B;}

  public:

    RemovePoints();

    virtual void operate();

    bool removePointCriteria( vtkIdType id_node ); ///< Check if a point needs to be removed

    int getNumRemoved() { return m_NumRemoved; }

    /// deletes set of points DeadNodes
  bool DeleteSetOfPoints( QVector <vtkIdType> deadnode_vector, QVector <vtkIdType> snappoint_vector, QSet <vtkIdType> & all_deadcells, QSet <vtkIdType> & all_mutatedcells, QSet <vtkIdType> & all_mutilatedcells, int& num_newpoints, int& num_newcells);
  
    /// returns a valid potential snappoint (checks for flipped cells, etc). If none is found, returns -1.
    vtkIdType FindSnapPoint( vtkIdType DeadNode, QSet <vtkIdType> & DeadCells, QSet <vtkIdType> & MutatedCells, QSet <vtkIdType> & MutilatedCells, int& N_newpoints, int& N_newcells );
  
    ///returns true if moving id_node to position P leads to flipped cells
    bool FlippedCells( vtkIdType id_node, vec3_t P );
  
    /// returns number of common neighbour nodes of id_node1 and id_node2. IsTetra becomes true if id_node1 and id_node2 belong to the edge of a tetrahedron.
    int NumberOfCommonPoints( vtkIdType id_node1, vtkIdType id_node2, bool& IsTetra );
  
};

#endif
