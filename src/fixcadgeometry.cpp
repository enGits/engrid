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
#include "fixcadgeometry.h"
#include "surfacemesher.h"
#include "vertexmeshdensity.h"
#include "guimainwindow.h"

fixCadGeometry::fixCadGeometry()
{
  EG_TYPENAME;
}

void fixCadGeometry::operate()
{
  qDebug()<<"==>fixing CAD geometry...";
  
  //update node info
  updateNodeInfo(true);
  
  //prepare BCmap
  QSet <int> bcs;
  GuiMainWindow::pointer()->getAllBoundaryCodes(bcs);
  QMap <int,int> BCmap;
  foreach(int bc, bcs) BCmap[bc]=1;
  
  //set density infinite
  QVector <VertexMeshDensity> VMDvector;
  VertexMeshDensity VMD;
  VMD.density = 9000;
  VMD.BCmap = BCmap;
  qDebug()<<"VMD="<<VMD;
  VMDvector.push_back(VMD);
  
  //call surface mesher
  SurfaceMesher mesher;
  mesher.setGrid(grid);
  mesher.setBoundaryCodes(bcs);
  mesher.setVertexMeshDensityVector(VMDvector);
  mesher();

  // finalize
  createIndices(grid);
  updateNodeInfo(false);
  computeMeshDensity();
}
