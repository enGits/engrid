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
#include "surfacealgorithm.h"

#include "insertpoints.h"
#include "removepoints.h"
#include "updatedesiredmeshdensity.h"
#include "smoothingutilities.h"
#include "swaptriangles.h"
#include "laplacesmoother.h"
#include "guimainwindow.h"


SurfaceAlgorithm::SurfaceAlgorithm()
{
  EG_TYPENAME;
  getSet("surface meshing", "maximal number of iterations",  5, m_NumMaxIter);
  getSet("surface meshing", "number of smoothing steps"   ,  2, m_NumSmoothSteps);
  getSet("surface meshing", "number of Delaunay sweeps"   ,  1, m_NumDelaunaySweeps);
  m_NodesPerQuarterCircle = 0;
  m_RespectFeatureEdgesForDeleteNodes = false;
  m_FeatureAngleForDeleteNodes = deg2rad(45);
  m_PerformGeometricTests = true;
  m_UseProjectionForSmoothing = true;
  m_UseNormalCorrectionForSmoothing = false;
  m_AllowFeatureEdgeSwapping = true;
  m_AllowSmallAreaSwapping = false;
  m_GrowthFactor = 1.5;
  m_FeatureResolution2D = 0;
  m_FeatureResolution3D = 0;
  setDeleteNodesOn();
  setInsertNodesOn();
}

void SurfaceAlgorithm::readSettings()
{
  QString buffer = GuiMainWindow::pointer()->getXmlSection("engrid/surface/settings").replace("\n", " ");
  QTextStream in(&buffer, QIODevice::ReadOnly);
  in >> m_MaxEdgeLength;
  in >> m_MinEdgeLength;
  in >> m_GrowthFactor;
  in >> m_NodesPerQuarterCircle;
  int num_bcs;
  in >> num_bcs;
  QVector<int> tmp_bcs;
  GuiMainWindow::pointer()->getAllBoundaryCodes(tmp_bcs);
  m_BoundaryCodes.clear();
  if (num_bcs == tmp_bcs.size()) {
    foreach (int bc, tmp_bcs) {
      int check_state;
      in >> check_state;
      if (check_state == 1) {
        m_BoundaryCodes.insert(bc);
      }
    }
  }
  if (!in.atEnd()) {
    in >> m_FeatureResolution2D;
    in >> m_FeatureResolution3D;
  }
}

void SurfaceAlgorithm::prepare()
{
  setAllCells();
  readSettings();
  
  EG_VTKDCN(vtkCharArray, node_type, m_Grid, "node_type");//node type
  
  updateNodeInfo();

}

void SurfaceAlgorithm::computeMeshDensity()
{
  ///\todo Optimize by using only one loop through nodes!
  UpdateDesiredMeshDensity update_desired_mesh_density;
  update_desired_mesh_density.setGrid(m_Grid);
  update_desired_mesh_density.setVertexMeshDensityVector(m_VMDvector);
  update_desired_mesh_density.setMaxEdgeLength(m_MaxEdgeLength);
  update_desired_mesh_density.setMinEdgeLength(m_MinEdgeLength);
  update_desired_mesh_density.setNodesPerQuarterCircle(m_NodesPerQuarterCircle);
  update_desired_mesh_density.setCellGrowthFactor(m_GrowthFactor);
  update_desired_mesh_density.setBoundaryCodes(m_BoundaryCodes);
  update_desired_mesh_density.setFeatureResolution2D(m_FeatureResolution2D);
  update_desired_mesh_density.setFeatureResolution3D(m_FeatureResolution3D);
  update_desired_mesh_density();
}

void SurfaceAlgorithm::swap(double delaunay_threshold, bool verbose)
{
  SwapTriangles swap;
  swap.setGrid(m_Grid);
  swap.setRespectBC(true);
  swap.setFeatureSwap(m_AllowFeatureEdgeSwapping);
  swap.setFeatureAngle(m_FeatureAngle);
  swap.setMaxNumLoops(m_NumDelaunaySweeps);
  swap.setSmallAreaSwap(m_AllowSmallAreaSwapping);
  swap.setBCodesFeatureDefinition(m_BCodeFeatureDefinition);
  swap.setDelaunayThreshold(delaunay_threshold);
  if (verbose) {
    swap.setVerboseOn();
  }
  QSet<int> rest_bcs = GuiMainWindow::pointer()->getAllBoundaryCodes();
  rest_bcs -= m_BoundaryCodes;
  swap.setBoundaryCodes(rest_bcs);
  swap();
}

void SurfaceAlgorithm::smooth(int N_iter, bool correct_curveture)
{
  LaplaceSmoother lap;
  lap.setGrid(m_Grid);
  QVector<vtkIdType> cls;
  getSurfaceCells(m_BoundaryCodes, cls, m_Grid);
  lap.setCells(cls);
  lap.setNumberOfIterations(N_iter);
  lap.setBoundaryCodes(m_BoundaryCodes);//IMPORTANT: so that unselected nodes become fixed when node types are updated!
  lap.setCorrectCurvature(correct_curveture);
  lap.setBCodesFeatureDefinition(m_BCodeFeatureDefinition);
  if (m_UseProjectionForSmoothing) {
    lap.setProjectionOn();
  } else {
    lap.setProjectionOff();
  }
  if (m_UseNormalCorrectionForSmoothing) {
    lap.setNormalCorrectionOn();
  } else {
    lap.setNormalCorrectionOff();
  }
  lap();
  m_SmoothSuccess = lap.succeeded();
}

int SurfaceAlgorithm::insertNodes()
{
  if (m_InsertNodes) {
    InsertPoints insert_points;
    insert_points.setGrid(m_Grid);
    insert_points.setBoundaryCodes(m_BoundaryCodes);
    insert_points.setBCodesFeatureDefinition(m_BCodeFeatureDefinition);
    insert_points();
    return insert_points.getNumInserted();
  }
  return 0;
}

int SurfaceAlgorithm::deleteNodes()
{
  if (m_DeleteNodes) {
    RemovePoints remove_points;
    remove_points.setGrid(m_Grid);
    remove_points.setBoundaryCodes(m_BoundaryCodes);
    remove_points.setStretchingFactor(m_StretchingFactor);
    remove_points.setFeatureAngle(m_FeatureAngle);
    remove_points.setBCodesFeatureDefinition(m_BCodeFeatureDefinition);
    if (m_RespectFeatureEdgesForDeleteNodes) {
      remove_points.setProtectFeatureEdgesOn();
    } else {
      remove_points.setProtectFeatureEdgesOff();
    }
    //remove_points.setFeatureAngle(m_FeatureAngleForDeleteNodes);
    if (m_PerformGeometricTests) {
      remove_points.setPerformGeometricChecksOn();
    } else {
      remove_points.setPerformGeometricChecksOff();
    }
    remove_points();
    return remove_points.getNumRemoved();
  }
  return 0;
}
