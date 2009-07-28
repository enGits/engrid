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
  m_NodesPerQuarterCircle = 0;
  m_RespectFeatureEdgesForDeleteNodes = false;
  m_FeatureAngleForDeleteNodes = deg2rad(45);
  m_PerformGeometricTests = true;
  m_UseProjectionForSmoothing = true;
}

void SurfaceAlgorithm::readVMD()
{
  QString buffer = GuiMainWindow::pointer()->getXmlSection("engrid/surface/table");
  QTextStream in(&buffer, QIODevice::ReadOnly);
  int row_count = 0;
  int column_count = 0;
  in >> row_count >> column_count;
  VMDvector.clear();
  if (column_count == m_BoundaryCodes.size() + 3) {
    VMDvector.fill(VertexMeshDensity(), row_count);
    for (int i = 0; i < row_count; ++i) {
      int row, column;
      QString formula;
      foreach (int bc, m_BoundaryCodes) {
        in >> row >> column >> formula;
        VMDvector[row].BCmap[bc] = formula.toInt();
      }
      in >> row >> column >> formula;
      VMDvector[row].type = Str2VertexType(formula);
      in >> row >> column >> formula;
      if (formula == "{{{empty}}}") {
        formula = "";
      }
      VMDvector[i].setNodes(formula);
      in >> row >> column >> formula;
      VMDvector[i].density = formula.toDouble();
    }
  }
}

void SurfaceAlgorithm::readSettings()
{
  QString buffer = GuiMainWindow::pointer()->getXmlSection("engrid/surface/settings");
  QTextStream in(&buffer, QIODevice::ReadOnly);
  QString str;
  in >> str;
  m_MaxEdgeLength = str.toDouble();
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
  readSettings();
  readVMD();
}

void SurfaceAlgorithm::computeMeshDensity()
{
  ///@@@  TODO: Optimize by using only one loop through nodes!
  UpdateDesiredMeshDensity update_desired_mesh_density;
  update_desired_mesh_density.setGrid(grid);
  update_desired_mesh_density.setVertexMeshDensityVector(VMDvector);
  update_desired_mesh_density.setMaxEdgeLength(m_MaxEdgeLength);
  update_desired_mesh_density.setNodesPerQuarterCircle(m_NodesPerQuarterCircle);
  update_desired_mesh_density();
}

void SurfaceAlgorithm::updateNodeInfo(bool update_type)
{
  setAllCells();
  l2g_t nodes = getPartNodes();
  foreach (vtkIdType id_node, nodes) {
    if(update_type) {
      EG_VTKDCN(vtkCharArray, node_type, grid, "node_type");//node type
      node_type->SetValue(id_node, getNodeType(id_node));
    }
    EG_VTKDCN(vtkDoubleArray, node_meshdensity_current, grid, "node_meshdensity_current");//what we have
    node_meshdensity_current->SetValue(id_node, CurrentVertexAvgDist(id_node));

    EG_VTKDCN(vtkIntArray, node_specified_density, grid, "node_specified_density");//density index from table
    VertexMeshDensity nodeVMD = getVMD(id_node);
    int idx=VMDvector.indexOf(nodeVMD);
    node_specified_density->SetValue(id_node, idx);
  }
}

void SurfaceAlgorithm::swap()
{
  SwapTriangles swap;
  swap.setGrid(grid);
  swap.setRespectBC(true);
  swap.setFeatureSwap(true);
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
  if (m_UseProjectionForSmoothing) {
    lap.setUseProjectionForSmoothingOn();
  } else {
    lap.setUseProjectionForSmoothingOff();
  }
  lap();
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

