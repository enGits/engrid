// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008-2014 enGits GmbH                                      +
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
#ifndef createvolumemesh_H
#define createvolumemesh_H

class CreateVolumeMesh;

#include "tetgenoperation.h"
#include "edgelengthsourcemanager.h"


class CreateVolumeMesh : public TetGenOperation
{
  
private: // attributes

  bool m_CreateBoundaryLayer;
  bool m_CreateVolumeMesh;
  bool m_FirstCall;

  vtkSmartPointer<vtkUnstructuredGrid> m_BackgroundGrid;
  

private: // methods
  
  int  numVolumeCells();
  

protected: // methods
  
  void createTetMesh(int max_num_passes, bool preserve_surface);
  void reduceSurface(QSet<int> boundary_codes);

  virtual void operate();
  

public: // methods
  
  CreateVolumeMesh();
  void setBoundaryLayerOn();
  void setVolumeMeshOn();

};

#endif
