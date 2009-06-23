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
#include "laplacesmoother.h"
#include <vtkCellLocator.h>
#include <vtkCharArray.h>
#include <vtkGenericCell.h>
#include "guimainwindow.h"

using namespace GeometryTools;

LaplaceSmoother::LaplaceSmoother() : SurfaceOperation()
{
  DebugLevel = 0;
  setQuickSave(true);
}

void LaplaceSmoother::operate()
{
  m_CartMesh.setBounds(vec3_t(-3,-3,-3), vec3_t(3,3,3));

  m_CartMesh.markToRefine(0);
  m_CartMesh.refineAll();

  for (int level = 0; level < 5; ++level) {
    for (int i = 0; i < m_CartMesh.getNumCells(); ++i) {
      double r_min =  1e99;
      double r_max = -1e99;
      for (int j = 0; j < 8; ++j) {
        r_min = min(r_min, m_CartMesh.getNodePoition(i,j).abs());
        r_max = max(r_max, m_CartMesh.getNodePoition(i,j).abs());
      }
      r_min -= 1;
      r_max -= 1;
      if (r_min*r_max <= 0) {
        m_CartMesh.markToRefine(i);
      }
    }
    m_CartMesh.refineAll();
    cout << m_CartMesh.getNumCells() << " cells" << endl;
  }
  EG_VTKSP(vtkUnstructuredGrid, grid);
  m_CartMesh.toVtkGrid(grid);
  writeGrid(grid, "octree");

  //UpdateNodeType();
}
