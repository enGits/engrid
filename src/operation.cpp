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
#include "operation.h"

#include "guimainwindow.h"
#include "egvtkobject.h"

#include <vtkTriangleFilter.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolygon.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkCellLocator.h>
#include <vtkFloatArray.h>
#include <vtkMath.h>
#include <vtkCharArray.h>

#include <QApplication>

#include "geometrytools.h"
using namespace GeometryTools;

QSet<Operation*> Operation::garbage_operations;

void Operation::collectGarbage()
{
  QSet<Operation*> delete_operations;

  foreach (Operation *op, garbage_operations)
  {
    if (!op->getThread().isRunning()) {
      delete_operations.insert(op);
      cout << "deleting Operation " << op << endl;
      delete op;
    }
  }

  foreach (Operation *op, delete_operations)
  {
    garbage_operations.remove(op);
  }
}

Operation::Operation()
{
  grid = NULL;
  gui = false;
  m_quicksave = false;
  m_resetoperationcounter = false;
  err = NULL;
  autoset = true;
  m_CellLocator = NULL;
  m_ProjectionSurface = NULL;
  
  //default values for determining node types and for smoothing operations
  ///@@@ TODO: Remove useless attributes
  Convergence=0;
  NumberOfIterations=20;
  RelaxationFactor=0.01;
  FeatureEdgeSmoothing=1;//0 by default in VTK, but we need 1 to avoid the "potatoe effect" ^^
  FeatureAngle=45;
  EdgeAngle=15;
  BoundarySmoothing=1;
  GenerateErrorScalars=0;
  GenerateErrorVectors=0;
}

Operation::~Operation()
{
  if (err) {
    err->display();
    delete err;
  }
}

void Operation::del() 
{ 
  garbage_operations.insert(this); 
}

void OperationThread::run()
{
  try {
    GuiMainWindow::lock();
    GuiMainWindow::pointer()->setBusy();
    op->operate();
  } catch (Error err) {
    op->err = new Error();
    *(op->err) = err;
  }
  GuiMainWindow::unlock();
  GuiMainWindow::pointer()->setIdle();
}

void Operation::operator()()
{
  if (gui) {
    if (GuiMainWindow::tryLock()) {
      checkGrid();
      thread.setOperation(this);
      GuiMainWindow::unlock();
      thread.start(QThread::LowPriority);
    } else {
      QMessageBox::warning(NULL, "not permitted", "Operation is not permitted while background process is running!");
    }
  } else {
    checkGrid();
    try {
      operate();
    } catch (Error err) {
      err.display();
    }
    if(m_resetoperationcounter) GuiMainWindow::pointer()->ResetOperationCounter();
    if(m_quicksave) GuiMainWindow::pointer()->QuickSave();
  }
}

void Operation::setAllCells()
{
  QVector<vtkIdType> all_cells;
  getAllCells(all_cells, grid);
  setCells(all_cells);
}

void Operation::setAllVolumeCells()
{
  QVector<vtkIdType> cells;
  getAllVolumeCells(cells, grid);
  setCells(cells);
}

void Operation::setAllSurfaceCells()
{
  QVector<vtkIdType> cells;
  getAllSurfaceCells(cells, grid);
  setCells(cells);
}

void Operation::initMapping()
{
  nodes_map.resize(nodes.size());
  for (int i_nodes = 0; i_nodes < nodes.size(); ++i_nodes) {
    nodes_map[i_nodes] = nodes[i_nodes];
  }
  cells_map.resize(cells.size());
  for (int i_cells = 0; i_cells < cells.size(); ++i_cells) {
    cells_map[i_cells] = cells[i_cells];
  }
}

void Operation::checkGrid()
{
  if (grid == NULL) {
    grid = GuiMainWindow::pointer()->getGrid();
  }
  if ((cells.size() == 0) && autoset) {
    setAllCells();
  }
}

void Operation::updateActors()
{
  mainWindow()->updateActors();
}

GuiMainWindow* Operation::mainWindow()
{
  return GuiMainWindow::pointer();
}

void Operation::populateBoundaryCodes(QListWidget *lw)
{
  QSet<int> bcs;
  mainWindow()->getAllBoundaryCodes(bcs);
  foreach(int bc, bcs) {
    QListWidgetItem *lwi = new QListWidgetItem(lw);
    lwi->setCheckState(Qt::Unchecked);
    QString text = "";
    QTextStream ts(&text);
    ts << bc;
    lwi->setText(text);
    lwi->setFlags(Qt::ItemIsUserCheckable | Qt::ItemIsEnabled);
  }
}

void Operation::populateVolumes(QListWidget *lw)
{
  QList<VolumeDefinition> vols = mainWindow()->getAllVols();
  foreach (VolumeDefinition V, vols) {
    QListWidgetItem *lwi = new QListWidgetItem(lw);
    lwi->setText(V.getName());
  }
}

stencil_t Operation::getStencil(vtkIdType id_cell1, int j1)
{
  if(grid->GetCellType(id_cell1)!=VTK_TRIANGLE) {
    cout<<"CELL IS NOT A TRIANGLE"<<endl;
    EG_BUG;
  }
  
  //return variable
  stencil_t S;
  
  //default values:
  S.sameBC = false;
  S.twocells = false;
  S.neighbour_type = -1;
  
  //initialize first cell
  S.id_cell1 = id_cell1;
  vtkIdType N1, *pts1;
  grid->GetCellPoints(S.id_cell1, N1, pts1);
  //place points 0,1,3
  if      (j1 == 0) { S.p[0] = pts1[2]; S.p[1] = pts1[0]; S.p[3] = pts1[1]; }
  else if (j1 == 1) { S.p[0] = pts1[0]; S.p[1] = pts1[1]; S.p[3] = pts1[2]; }
  else if (j1 == 2) { S.p[0] = pts1[1]; S.p[1] = pts1[2]; S.p[3] = pts1[0]; };
  
  //initialize second cell
  S.id_cell2 = -1;
  S.p[2]=-1;
  
  //twocells
  if (c2c[_cells[id_cell1]][j1] != -1) {//if neighbour cell
    
    //twocells
    S.twocells = true;
    S.id_cell2 = cells[c2c[_cells[id_cell1]][j1]];
    
    //sameBC
    EG_VTKDCC(vtkIntArray, cell_code, grid, "cell_code");
    if(cell_code->GetValue(S.id_cell1)==cell_code->GetValue(S.id_cell2)) S.sameBC = true;
    
    //neighbour_type
    S.neighbour_type = grid->GetCellType(S.id_cell2);
    if ( S.neighbour_type == VTK_TRIANGLE) {//if neighbour cell is a triangle
      vtkIdType N2, *pts2;
      grid->GetCellPoints(S.id_cell2, N2, pts2);
      
      //place point 2
      bool p2 = false;
      if (c2c[_cells[S.id_cell2]][0] != -1) {
        if (cells[c2c[_cells[S.id_cell2]][0]] == S.id_cell1) {
          S.p[2] = pts2[2];
          p2 = true;
        }
      }
      if (c2c[_cells[S.id_cell2]][1] != -1) {
        if (cells[c2c[_cells[S.id_cell2]][1]] == S.id_cell1) {
          S.p[2] = pts2[0];
          p2 = true;
        }
      }
      if (c2c[_cells[S.id_cell2]][2] != -1) {
        if (cells[c2c[_cells[S.id_cell2]][2]] == S.id_cell1) {
          S.p[2] = pts2[1];
          p2 = true;
        }
      }
      
      if (!p2) {//failed to place point 2, appears when cell1 is linked to cell2, but cell2 not to cell1
        cout<<"S.id_cell1="<<S.id_cell1<<endl;
        cout<<"S.id_cell2="<<S.id_cell2<<endl;
        createNodeToCell(cells, nodes, _nodes, n2c, grid);
        EG_BUG;
      }
    }
  }//end of if neighbour cell
  return S;
}

ostream& operator<<(ostream &out, stencil_t S)
{
  out<<"S.id_cell1="<<S.id_cell1<<" ";
  out<<"S.id_cell2="<<S.id_cell2<<" ";
  out<<"S.sameBC="<<S.sameBC<<" ";
  out<<"S.twocells="<<S.twocells<<" ";
  out<<"S.neighbour_type="<<S.neighbour_type<<" ";
  out<<"[";
  for(int i=0;i<4;i++){
    out<<S.p[i];
    if(i!=3) out<<",";
  }
  out<<"]";
  return(out);
}

vtkIdType Operation::getClosestNode(vtkIdType a_id_node,vtkUnstructuredGrid* a_grid)
{
  vec3_t C;
  a_grid->GetPoint(a_id_node,C.data());
  vtkIdType id_minlen=-1;
  double minlen=-1;
  foreach(vtkIdType neighbour,n2n[a_id_node])
  {
    vec3_t M;
    a_grid->GetPoint(neighbour,M.data());
    double len=(M-C).abs();
    if(minlen<0 or len<minlen)
    {
      minlen=len;
      id_minlen=neighbour;
    }
  }
  return(id_minlen);
}

vtkIdType Operation::getFarthestNode(vtkIdType a_id_node,vtkUnstructuredGrid* a_grid)
{
  vec3_t C;
  a_grid->GetPoint(a_id_node,C.data());
  vtkIdType id_maxlen=-1;
  double maxlen=-1;
  foreach(vtkIdType neighbour,n2n[a_id_node])
  {
    vec3_t M;
    a_grid->GetPoint(neighbour,M.data());
    double len=(M-C).abs();
    if(maxlen<0 or len>maxlen)
    {
      maxlen=len;
      id_maxlen=neighbour;
    }
  }
  return(id_maxlen);
}

void Operation::quad2triangle(vtkUnstructuredGrid* src,vtkIdType quadcell,vtkIdType MovingPoint)
{
  vtkIdType type_cell = grid->GetCellType(quadcell);
  cout<<"It's a "<<type_cell<<endl;
  if(type_cell==VTK_QUAD)
  {
    cout_grid(cout,src,true,true,true,true);
    EG_VTKSP(vtkUnstructuredGrid, dst);
    //src grid info
    int N_points=src->GetNumberOfPoints();
    int N_cells=src->GetNumberOfCells();
    allocateGrid(dst,N_cells+1,N_points);
    
    for (vtkIdType id_node = 0; id_node < src->GetNumberOfPoints(); ++id_node) {
      vec3_t x;
      src->GetPoints()->GetPoint(id_node, x.data());
      dst->GetPoints()->SetPoint(id_node, x.data());
      copyNodeData(src, id_node, dst, id_node);
    }
    for (vtkIdType id_cell = 0; id_cell < src->GetNumberOfCells(); ++id_cell) {
      vtkIdType N_pts, *pts;
      src->GetCellPoints(id_cell, N_pts, pts);
      vtkIdType type_cell = src->GetCellType(id_cell);
      vtkIdType id_new_cell;
      vtkIdType id_new_cell1;
      vtkIdType id_new_cell2;
      if(id_cell!=quadcell)
      {
        id_new_cell = dst->InsertNextCell(type_cell, N_pts, pts);
        copyCellData(src, id_cell, dst, id_new_cell);
      }
      else
      {
        vtkIdType triangle1[3];
        vtkIdType triangle2[3];
        if(MovingPoint==pts[0] || MovingPoint==pts[2])
        {
          triangle1[0]=pts[1];
          triangle1[1]=pts[3];
          triangle1[2]=pts[0];
          triangle2[0]=pts[3];
          triangle2[1]=pts[1];
          triangle2[2]=pts[2];
        }
        else
        {
          triangle1[0]=pts[2];
          triangle1[1]=pts[0];
          triangle1[2]=pts[1];
          triangle2[0]=pts[0];
          triangle2[1]=pts[2];
          triangle2[2]=pts[3];
        }
        id_new_cell1 = dst->InsertNextCell(VTK_TRIANGLE, 3, triangle1);
        copyCellData(src, id_cell, dst, id_new_cell1);
        id_new_cell2 = dst->InsertNextCell(VTK_TRIANGLE, 3, triangle2);
        copyCellData(src, id_cell, dst, id_new_cell2);
      }
    }
    cout_grid(cout,dst,true,true,true,true);
    makeCopy(dst, src);
  }//end of if quad
}
