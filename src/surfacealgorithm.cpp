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
  getSet("surface meshing", "maximal number of iterations", 20, m_NumMaxIter);
  getSet("surface meshing", "number of smoothing steps"   ,  1, m_NumSmoothSteps);
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
}

void SurfaceAlgorithm::readVMD()
{
  QString buffer = GuiMainWindow::pointer()->getXmlSection("engrid/surface/table").replace("\n", " ");
  QTextStream in(&buffer, QIODevice::ReadOnly);
  int row_count = 0;
  int column_count = 0;
  in >> row_count >> column_count;
  QSet<int> tmp_bcs;
  GuiMainWindow::pointer()->getAllBoundaryCodes(tmp_bcs);
  m_VMDvector.clear();
  if (column_count == tmp_bcs.size() + 3) {
    m_VMDvector.fill(VertexMeshDensity(), row_count);
    for (int i = 0; i < row_count; ++i) {
      int row, column;
      QString formula;
      foreach (int bc, tmp_bcs) {
        in >> row >> column >> formula;
        m_VMDvector[row].BCmap[bc] = formula.toInt();
      }
      in >> row >> column >> formula;
      m_VMDvector[row].type = Str2VertexType(formula);
      in >> row >> column >> formula;
      if (formula == "{{{empty}}}") {
        formula = "";
      }
      m_VMDvector[i].setNodes(formula);
      in >> row >> column >> formula;
      m_VMDvector[i].density = formula.toDouble();
      cout << m_VMDvector[i] << endl;
    }
  } else {
    EG_ERR_RETURN("The number of boundary conditions don't match between the mesh and the settings table");
  }
}

void SurfaceAlgorithm::readSettings()
{
  QString buffer = GuiMainWindow::pointer()->getXmlSection("engrid/surface/settings").replace("\n", " ");
  QTextStream in(&buffer, QIODevice::ReadOnly);
  in >> m_MaxEdgeLength;
  in >> m_GrowthFactor;
  in >> m_NodesPerQuarterCircle;
  int num_bcs;
  in >> num_bcs;
  QSet<int> tmp_bcs;
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
}

void SurfaceAlgorithm::prepare()
{
  setAllCells();
  readSettings();
  readVMD();
  updateNodeInfo(true);
}

void SurfaceAlgorithm::computeMeshDensity()
{
  ///\todo Optimize by using only one loop through nodes!
  UpdateDesiredMeshDensity update_desired_mesh_density;
  update_desired_mesh_density.setGrid(grid);
  update_desired_mesh_density.setVertexMeshDensityVector(m_VMDvector);
  update_desired_mesh_density.setMaxEdgeLength(m_MaxEdgeLength);
  update_desired_mesh_density.setNodesPerQuarterCircle(m_NodesPerQuarterCircle);
  update_desired_mesh_density.setCellGrowthFactor(m_GrowthFactor);
  update_desired_mesh_density.setBoundaryCodes(m_BoundaryCodes);
  update_desired_mesh_density();
}

void SurfaceAlgorithm::updateNodeInfo(bool update_type)
{
  setAllCells();
  l2g_t nodes = getPartNodes();
  foreach (vtkIdType id_node, nodes) {
    if(update_type) {
      EG_VTKDCN(vtkCharArray, node_type, grid, "node_type");//node type
      node_type->SetValue(id_node, getNodeType(id_node, true));
    }
    EG_VTKDCN(vtkDoubleArray, node_meshdensity_current, grid, "node_meshdensity_current");//what we have
    node_meshdensity_current->SetValue(id_node, CurrentVertexAvgDist(id_node));

    EG_VTKDCN(vtkIntArray, node_specified_density, grid, "node_specified_density");//density index from table
    VertexMeshDensity nodeVMD = getVMD(id_node);
    int idx = m_VMDvector.indexOf(nodeVMD);
//     int idx = nodeVMD.findSmallestVMD(m_VMDvector);
//     qWarning()<<"idx="<<idx;
    node_specified_density->SetValue(id_node, idx);
  }
  writeGrid(grid, "info");
}

void SurfaceAlgorithm::swap()
{
  SwapTriangles swap;
  swap.setGrid(grid);
  swap.setRespectBC(true);
  swap.setFeatureSwap(m_AllowFeatureEdgeSwapping);
  swap.setFeatureAngle(m_FeatureAngle);
  swap.setMaxNumLoops(m_NumDelaunaySweeps);
  swap.setSmallAreaSwap(m_AllowSmallAreaSwapping);
  QSet<int> rest_bcs;
  GuiMainWindow::pointer()->getAllBoundaryCodes(rest_bcs);
  rest_bcs -= m_BoundaryCodes;
  swap.setBoundaryCodes(rest_bcs);
  swap();
}

void SurfaceAlgorithm::smooth(int N_iter)
{
  LaplaceSmoother lap;
  lap.setGrid(grid);
  QVector<vtkIdType> cls;
  getSurfaceCells(m_BoundaryCodes, cls, grid);
  lap.setCells(cls);
  lap.setNumberOfIterations(N_iter);
  lap.setBoundaryCodes(m_BoundaryCodes);//IMPORTANT: so that unselected nodes become fixed when node types are updated!
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
  InsertPoints insert_points;
  insert_points.setGrid(grid);
  insert_points.setBoundaryCodes(m_BoundaryCodes);
  insert_points();
  return insert_points.getNumInserted();
}

int SurfaceAlgorithm::deleteNodes()
{
  RemovePoints remove_points;
  remove_points.setGrid(grid);
  remove_points.setBoundaryCodes(m_BoundaryCodes);
  if (m_RespectFeatureEdgesForDeleteNodes) {
    remove_points.setProtectFeatureEdgesOn();
  } else {
    remove_points.setProtectFeatureEdgesOff();
  }
  remove_points.setFeatureAngle(m_FeatureAngleForDeleteNodes);
  if (m_PerformGeometricTests) {
    remove_points.setPerformGeometricChecksOn();
  } else {
    remove_points.setPerformGeometricChecksOff();
  }
  remove_points();
  return remove_points.getNumRemoved();
}

