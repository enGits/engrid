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
  foreach (Operation *op, garbage_operations) {
    if (!op->getThread().isRunning()) {
      delete_operations.insert(op);
      cout << "deleting Operation " << op << endl;
      delete op;
    };
  };
  foreach (Operation *op, delete_operations) {
    garbage_operations.remove(op);
  };
};

Operation::Operation()
{
  grid = NULL;
  gui = false;
  m_quicksave = false;
  m_resetoperationcounter = false;
  err = NULL;
  autoset = true;

  //default values for determining node types and for smoothing operations
  Convergence=0;
  NumberOfIterations=20;
  RelaxationFactor=0.01;
  FeatureEdgeSmoothing=1;//0 by default in VTK, but we need 1 to avoid the "potatoe effect" ^^
  FeatureAngle=45;
  EdgeAngle=15;
  BoundarySmoothing=1;
  GenerateErrorScalars=0;
  GenerateErrorVectors=0;
};

Operation::~Operation()
{
  if (err) {
    err->display();
    delete err;
  };
};

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
  };
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
    };
  } else {
    checkGrid();
    try {
      operate();
    } catch (Error err) {
      err.display();
    }
    if(m_resetoperationcounter) GuiMainWindow::pointer()->ResetOperationCounter();
    if(m_quicksave) GuiMainWindow::pointer()->QuickSave();
  };
}

void Operation::setAllCells()
{
  QVector<vtkIdType> all_cells;
  getAllCells(all_cells, grid);
  setCells(all_cells);
};

void Operation::setAllVolumeCells()
{
  QVector<vtkIdType> cells;
  getAllVolumeCells(cells, grid);
  setCells(cells);
};

void Operation::setAllSurfaceCells()
{
  QVector<vtkIdType> cells;
  getAllSurfaceCells(cells, grid);
  setCells(cells);
};

void Operation::initMapping()
{
  nodes_map.resize(nodes.size());
  for (int i_nodes = 0; i_nodes < nodes.size(); ++i_nodes) {
    nodes_map[i_nodes] = nodes[i_nodes];
  };
  cells_map.resize(cells.size());
  for (int i_cells = 0; i_cells < cells.size(); ++i_cells) {
    cells_map[i_cells] = cells[i_cells];
  };
};

void Operation::checkGrid()
{
  if (grid == NULL) {
    grid = GuiMainWindow::pointer()->getGrid();
  };
  if ((cells.size() == 0) && autoset) {
    setAllCells();
  };
};

void Operation::updateActors()
{
  mainWindow()->updateActors();
};

GuiMainWindow* Operation::mainWindow()
{
  return GuiMainWindow::pointer();
};

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
  };
};

//TODO: Get the EG_BUG error again and figure out where it came from
stencil_t Operation::getStencil(vtkIdType id_cell1, int j1)
{
  stencil_t S;
  S.valid = true;
  S.id_cell1 = id_cell1;
  if (c2c[_cells[id_cell1]][j1] != -1) {
    S.id_cell2 = cells[c2c[_cells[id_cell1]][j1]];
/*    if(S.id_cell1==138)
    {
      cout<<"id_cell1="<<id_cell1<<endl;
      cout<<"j1="<<j1<<endl;
      cout<<"cells[c2c[_cells[id_cell1]][j1]]="<<cells[c2c[_cells[id_cell1]][j1]]<<endl;
      cout<<"c2c[id_cell1][j1]="<<c2c[id_cell1][j1]<<endl;
      cout<<"c2c[_cells[id_cell1]]="<<c2c[_cells[id_cell1]]<<endl;
      cout<<"c2c[id_cell1]="<<c2c[id_cell1]<<endl;
      cout<<"S.id_cell1="<<S.id_cell1<<endl;
      cout<<"S.id_cell2="<<S.id_cell2<<endl;
    }*/
    if (grid->GetCellType(S.id_cell2) != VTK_TRIANGLE) {
      EG_BUG;
    };
    vtkIdType N1, N2, *pts1, *pts2;
    grid->GetCellPoints(S.id_cell1, N1, pts1);
    grid->GetCellPoints(S.id_cell2, N2, pts2);
    if      (j1 == 0) { S.p[0] = pts1[2]; S.p[1] = pts1[0]; S.p[3] = pts1[1]; }
    else if (j1 == 1) { S.p[0] = pts1[0]; S.p[1] = pts1[1]; S.p[3] = pts1[2]; }
    else if (j1 == 2) { S.p[0] = pts1[1]; S.p[1] = pts1[2]; S.p[3] = pts1[0]; };
    bool p2 = false;
    if (c2c[_cells[S.id_cell2]][0] != -1) {
      if (cells[c2c[_cells[S.id_cell2]][0]] == S.id_cell1) {
        S.p[2] = pts2[2];
        p2 = true;
      };
    };
    if (c2c[_cells[S.id_cell2]][1] != -1) {
      if (cells[c2c[_cells[S.id_cell2]][1]] == S.id_cell1) {
        S.p[2] = pts2[0];
        p2 = true;
      };
    };
    if (c2c[_cells[S.id_cell2]][2] != -1) {
      if (cells[c2c[_cells[S.id_cell2]][2]] == S.id_cell1) {
        S.p[2] = pts2[1];
        p2 = true;
      };
    };
    if (!p2) {
      DualSave("/data1/home/mtaverne/Geometries/simulations/SurfaceTests/abort");
      cout<<"S.id_cell1="<<S.id_cell1<<endl;
      cout<<"S.id_cell2="<<S.id_cell2<<endl;
      EG_BUG;
    };
  } else {
    S.valid = false;
    S.id_cell2 = -1;
    vtkIdType N1, *pts1;
    grid->GetCellPoints(S.id_cell1, N1, pts1);
    if      (j1 == 0) { S.p[0] = pts1[2]; S.p[1] = pts1[0]; S.p[3] = pts1[1]; }
    else if (j1 == 1) { S.p[0] = pts1[0]; S.p[1] = pts1[1]; S.p[3] = pts1[2]; }
    else if (j1 == 2) { S.p[0] = pts1[1]; S.p[1] = pts1[2]; S.p[3] = pts1[0]; };
  };
  return S;
};

ostream& operator<<(ostream &out, stencil_t S)
{
  out<<"S.id_cell1="<<S.id_cell1<<" ";
  out<<"S.id_cell2="<<S.id_cell2<<" ";
  out<<"S.valid="<<S.valid<<" ";
  out<<"[";
  for(int i=0;i<4;i++){
    out<<S.p[i];
    if(i!=3) out<<",";
  }
  out<<"]";
  return(out);
}

//////////////////////////////////////////////
double Operation::CurrentVertexAvgDist(vtkIdType a_vertex)
{
  double total_dist=0;
  double avg_dist=0;
  int N=n2n_func(a_vertex).size();
  vec3_t C;
  grid->GetPoint(a_vertex, C.data());
  foreach(int i,n2n_func(a_vertex))
  {
    vec3_t M;
    grid->GetPoint(i, M.data());
    total_dist+=(M-C).abs();
  }
  avg_dist=total_dist/(double)N;
  return(avg_dist);
}

double Operation::CurrentMeshDensity(vtkIdType a_vertex)
{
  double total_dist=0;
  double avg_dist=0;
  int N=n2n_func(a_vertex).size();
  vec3_t C;
  grid->GetPoint(a_vertex, C.data());
  foreach(int i,n2n_func(a_vertex))
  {
    vec3_t M;
    grid->GetPoint(i, M.data());
    total_dist+=(M-C).abs();
  }
  avg_dist=total_dist/(double)N;
  double avg_density=1./avg_dist;
  return(avg_density);
}

double Operation::DesiredVertexAvgDist(vtkIdType a_vertex)
{
  double total_dist=0;
  double avg_dist=0;
  EG_VTKDCN(vtkDoubleArray, node_meshdensity, grid, "node_meshdensity");
  int N=n2n_func(a_vertex).size();
  foreach(int i,n2n_func(a_vertex))
  {
    total_dist+=1./node_meshdensity->GetValue(i);
  }
  avg_dist=total_dist/(double)N;
  return(avg_dist);
}

double Operation::DesiredMeshDensity(vtkIdType a_vertex)
{
  double total_density=0;
  double avg_density=0;
  EG_VTKDCN(vtkDoubleArray, node_meshdensity, grid, "node_meshdensity");
  int N=n2n_func(a_vertex).size();
  foreach(int i,n2n_func(a_vertex))
  {
    total_density+=node_meshdensity->GetValue(i);
  }
  avg_density=total_density/(double)N;
  return(avg_density);
}

///////////////////////////////////////////
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

bool Operation::SwapCells(vtkUnstructuredGrid* a_grid, stencil_t S)
{
  bool swap = false;
  if (S.valid) {
    vec3_t x3[4], x3_0(0,0,0);
    vec2_t x[4];
    for (int k = 0; k < 4; ++k) {
      a_grid->GetPoints()->GetPoint(S.p[k], x3[k].data());
      x3_0 += x3[k];
    };
    vec3_t n1 = triNormal(x3[0], x3[1], x3[3]);
    vec3_t n2 = triNormal(x3[1], x3[2], x3[3]);
    n1.normalise();
    n2.normalise();
    if ( (n1*n2) > 0.8) {
      vec3_t n = n1 + n2;
      n.normalise();
      vec3_t ex = orthogonalVector(n);
      vec3_t ey = ex.cross(n);
      for (int k = 0; k < 4; ++k) {
        x[k] = vec2_t(x3[k]*ex, x3[k]*ey);
      };
      vec2_t r1, r2, r3, u1, u2, u3;
      r1 = 0.5*(x[0] + x[1]); u1 = turnLeft(x[1] - x[0]);
      r2 = 0.5*(x[1] + x[2]); u2 = turnLeft(x[2] - x[1]);
      r3 = 0.5*(x[1] + x[3]); u3 = turnLeft(x[3] - x[1]);
      double k, l;
      vec2_t xm1, xm2;
      bool ok = true;
      if (intersection(k, l, r1, u1, r3, u3)) {
        xm1 = r1 + k*u1;
        if (intersection(k, l, r2, u2, r3, u3)) {
          xm2 = r2 + k*u2;
        } else {
          ok = false;
        };
      } else {
        ok = false;
        swap = true;
      };
      if (ok) {
        if ((xm1 - x[2]).abs() < (xm1 - x[0]).abs()) {
          swap = true;
        };
        if ((xm2 - x[0]).abs() < (xm2 - x[2]).abs()) {
          swap = true;
        };
      };
    };
  };
  if (swap) {
    vtkIdType new_pts1[3], new_pts2[3];
    new_pts1[0] = S.p[1];
    new_pts1[1] = S.p[2];
    new_pts1[2] = S.p[0];
    new_pts2[0] = S.p[2];
    new_pts2[1] = S.p[3];
    new_pts2[2] = S.p[0];
    a_grid->ReplaceCell(S.id_cell1, 3, new_pts1);
    a_grid->ReplaceCell(S.id_cell2, 3, new_pts2);
  };
  return(swap);
}

void Operation::quad2triangle(vtkUnstructuredGrid* src,vtkIdType quadcell)
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
    };
    for (vtkIdType id_cell = 0; id_cell < src->GetNumberOfCells(); ++id_cell) {
      vtkIdType N_pts, *pts;
      vtkIdType type_cell = src->GetCellType(id_cell);
      src->GetCellPoints(id_cell, N_pts, pts);
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
        triangle1[0]=pts[1];
        triangle1[1]=pts[3];
        triangle1[2]=pts[0];
        triangle2[0]=pts[3];
        triangle2[1]=pts[1];
        triangle2[2]=pts[2];
        id_new_cell1 = dst->InsertNextCell(VTK_TRIANGLE, 3, triangle1);
        copyCellData(src, id_cell, dst, id_new_cell1);
        id_new_cell2 = dst->InsertNextCell(VTK_TRIANGLE, 3, triangle2);
        copyCellData(src, id_cell, dst, id_new_cell2);
        stencil_t S;
        S.id_cell1=id_new_cell1;
        S.id_cell2=id_new_cell2;
        S.p[0]=pts[0];
        S.p[1]=pts[1];
        S.p[2]=pts[2];
        S.p[3]=pts[3];
        S.valid=true;
        SwapCells(dst,S);
      }
    };
    cout_grid(cout,dst,true,true,true,true);
    makeCopy(dst, src);
  }//end of if quad
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
    };
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
    };
    cout_grid(cout,dst,true,true,true,true);
    makeCopy(dst, src);
  }//end of if quad
}

int Operation::NumberOfCommonPoints(vtkIdType node1, vtkIdType node2, bool& IsTetra)
{
//   QVector< QSet< int > > 	n2n
  QSet <int> node1_neighbours=n2n[node1];
  QSet <int> node2_neighbours=n2n[node2];
  QSet <int> intersection=node1_neighbours.intersect(node2_neighbours);
  int N=intersection.size();
  IsTetra=false;
  if(N==2)
  {
    QSet<int>::const_iterator p1=intersection.begin();
    QSet<int>::const_iterator p2=p1+1;
    vtkIdType intersection1=_nodes[*p1];
    vtkIdType intersection2=_nodes[*p2];
    if(n2n[intersection1].contains(intersection2))//if there's an edge between intersection1 and intersection2
    {
      //check if (node1,intersection1,intersection2) and (node2,intersection1,intersection2) are defined as cells!
  //     QVector< QSet< int > > 	n2c
      QSet< int > S1=n2c[intersection1];
      QSet< int > S2=n2c[intersection2];
      QSet< int > Si=S1.intersect(S2);
      int counter=0;
      foreach(vtkIdType C,Si){
        vtkIdType N_pts, *pts;
        grid->GetCellPoints(C, N_pts, pts);
        for(int i=0;i<N_pts;i++)
        {
          if(pts[i]==node1 || pts[i]==node2) counter++;
        }
      }
      if(counter>=2) IsTetra=true;
    }
  }
  return(N);
}

//Remove or finish???
//Function to check if empty volumes appear when moving DeadNode tp PSP
bool Operation::EmptyVolume(vtkIdType DeadNode, vtkIdType PSP)
{
  c2c[DeadNode];
  c2c[PSP];
  return(true);
}

vec3_t Operation::GetCenter(vtkIdType cellId, double& R)
{
  vtkIdType *pts, Npts;
  grid->GetCellPoints(cellId, Npts, pts);
  
  vec3_t x(0,0,0);
  for (vtkIdType i = 0; i < Npts; ++i) {
    vec3_t xp;
    grid->GetPoints()->GetPoint(pts[i], xp.data());
    x += double(1)/Npts * xp;
  };
  
  R = 1e99;
  for (vtkIdType i = 0; i < Npts; ++i) {
    vec3_t xp;
    grid->GetPoints()->GetPoint(pts[i], xp.data());
    R = min(R, 0.25*(xp-x).abs());
  };
  
  return(x);
}

bool Operation::getNeighbours(vtkIdType Boss, QVector <vtkIdType>& Peons, int BC)
{
//   QVector <vtkIdType> Peons;
  
  QSet <int> S1=n2c[Boss];
//   cout<<"S1="<<S1<<endl;
  foreach(vtkIdType PN,n2n[Boss])
  {
//     cout<<"PN="<<PN<<endl;
    QSet <int> S2=n2c[PN];
//     cout<<"S2="<<S2<<endl;
    QSet <int> Si=S2.intersect(S1);
//     cout<<"PN="<<PN<<" Si="<<Si<<endl;
    if(Si.size()<2)//only one common cell
    {
      Peons.push_back(PN);
    }
    else
    {
      QSet <int> bc_set;
      foreach(vtkIdType C,Si)
      {
        EG_VTKDCC(vtkIntArray, cell_code, grid, "cell_code");
        int bc=cell_code->GetValue(C);
//         cout<<"C="<<C<<" bc="<<bc<<endl;
        bc_set.insert(bc);
      }
      if(bc_set.size()>1)//2 different boundary codes
      {
        Peons.push_back(PN);
      }
    }
  }
  if(Peons.size()==2)
  {
/*    Peon1=Peons[0];
    Peon2=Peons[1];*/
    return(true);
  }
  else
  {
    int N=n2n[Boss].size();
    QVector <vtkIdType> neighbours(N);
    qCopy(n2n[Boss].begin(), n2n[Boss].end(), neighbours.begin());
    
    double alphamin_value;
    vtkIdType alphamin_i;
    vtkIdType alphamin_j;
    bool first=true;
    
    for(int i=0;i<N;i++)
    {
      for(int j=i+1;j<N;j++)
      {
        double alpha=deviation(grid,neighbours[i],Boss,neighbours[j]);
//         cout<<"alpha("<<neighbours[i]<<","<<Boss<<","<<neighbours[j]<<")="<<alpha<<endl;
        if(first) {
          alphamin_value=alpha;
          alphamin_i=i;
          alphamin_j=j;
          first=false;
        }
        else
        {
          if(alpha<alphamin_value)
          {
            alphamin_value=alpha;
            alphamin_i=i;
            alphamin_j=j;
          }
        }
      }
    }
//     cout<<"alphamin_value="<<alphamin_value<<endl;
    
    Peons.resize(2);
    Peons[0]=neighbours[alphamin_i];
    Peons[1]=neighbours[alphamin_j];
    return(true);
/*    cout<<"FATAL ERROR: number of neighbours != 2"<<endl;
    EG_BUG;*/
  }
  return(false);//should never happen
}

bool Operation::getNeighbours_BC(vtkIdType Boss, QVector <vtkIdType>& Peons)
{
//   QVector <vtkIdType> Peons;
  
  QSet <int> S1=n2c[Boss];
//   cout<<"S1="<<S1<<endl;
  foreach(vtkIdType PN,n2n[Boss])
  {
//     cout<<"PN="<<PN<<endl;
    QSet <int> S2=n2c[PN];
//     cout<<"S2="<<S2<<endl;
    QSet <int> Si=S2.intersect(S1);
//     cout<<"PN="<<PN<<" Si="<<Si<<endl;
    if(Si.size()<2)//only one common cell
    {
      Peons.push_back(PN);
    }
    else
    {
      QSet <int> bc_set;
      foreach(vtkIdType C,Si)
      {
        EG_VTKDCC(vtkIntArray, cell_code, grid, "cell_code");
        int bc=cell_code->GetValue(C);
//         cout<<"C="<<C<<" bc="<<bc<<endl;
        bc_set.insert(bc);
      }
      if(bc_set.size()>1)//2 different boundary codes
      {
        Peons.push_back(PN);
      }
    }
  }
  if(Peons.size()==2)
  {
/*    Peon1=Peons[0];
    Peon2=Peons[1];*/
    return(true);
  }
  else
  {
    int N=n2n[Boss].size();
    QVector <vtkIdType> neighbours(N);
    qCopy(n2n[Boss].begin(), n2n[Boss].end(), neighbours.begin());
    
    double alphamin_value;
    vtkIdType alphamin_i;
    vtkIdType alphamin_j;
    bool first=true;
    
    for(int i=0;i<N;i++)
    {
      for(int j=i+1;j<N;j++)
      {
        double alpha=deviation(grid,neighbours[i],Boss,neighbours[j]);
//         cout<<"alpha("<<neighbours[i]<<","<<Boss<<","<<neighbours[j]<<")="<<alpha<<endl;
        if(first) {
          alphamin_value=alpha;
          alphamin_i=i;
          alphamin_j=j;
          first=false;
        }
        else
        {
          if(alpha<alphamin_value)
          {
            alphamin_value=alpha;
            alphamin_i=i;
            alphamin_j=j;
          }
        }
      }
    }
//     cout<<"alphamin_value="<<alphamin_value<<endl;
    
    Peons.resize(2);
    Peons[0]=neighbours[alphamin_i];
    Peons[1]=neighbours[alphamin_j];
    return(true);
/*    cout<<"FATAL ERROR: number of neighbours != 2"<<endl;
    EG_BUG;*/
  }
  return(false);//should never happen
}

int Operation::UpdateMeshDensity()
{
  if(DebugLevel>0) cout<<"===UpdateMeshDensity START==="<<endl;
  
  getAllSurfaceCells(cells,grid);
  EG_VTKDCC(vtkIntArray, cell_code, grid, "cell_code");
  EG_VTKDCN(vtkDoubleArray, node_meshdensity, grid, "node_meshdensity");
  getNodesFromCells(cells, nodes, grid);
  setGrid(grid);
  setCells(cells);
  
  if(DebugLevel>5) cout<<"cells.size()="<<cells.size()<<endl;
  
  EG_VTKDCN(vtkDoubleArray, node_meshdensity_current, grid, "node_meshdensity_current");
  foreach(vtkIdType node,nodes)
  {
    double L=CurrentVertexAvgDist(node);
    double D=1./L;
    node_meshdensity_current->SetValue(node, D);
  }
  if(DebugLevel>0) cout<<"===UpdateMeshDensity END==="<<endl;
  return(0);
}

// Special structure for marking vertices
typedef struct _vtkMeshVertex 
{
  char      type;
  vtkIdList *edges; // connected edges (list of connected point ids)
} vtkMeshVertex, *vtkMeshVertexPtr;

int Operation::UpdateNodeType_all()
{
  if(DebugLevel>0) cout<<"===UpdateNodeType_all START==="<<endl;
  if(DebugLevel>47) cout<<"this->FeatureAngle="<<this->FeatureAngle<<endl;
  if(DebugLevel>47) cout<<"this->EdgeAngle="<<this->EdgeAngle<<endl;
  
  getAllSurfaceCells(cells,grid);
  if(DebugLevel>5) cout<<"cells.size()="<<cells.size()<<endl;
  
  EG_VTKSP(vtkPolyData, pdata);
  //   addToPolyData(m_SelectedCells, pdata, grid);
  addToPolyData(cells, pdata, grid);
  
  vtkPolyData* input=pdata;
  
  vtkPolyData *source = 0;
  
  vtkIdType numPts, numCells, i, numPolys;
  int j, k;
  vtkIdType npts = 0;
  vtkIdType *pts = 0;
  vtkIdType p1, p2;
  double x[3], y[3], deltaX[3], xNew[3], conv, maxDist, dist, factor;
  double x1[3], x2[3], x3[3], l1[3], l2[3];
  double CosFeatureAngle; //Cosine of angle between adjacent polys
  double CosEdgeAngle; // Cosine of angle between adjacent edges
  double closestPt[3], dist2, *w = NULL;
  int iterationNumber, abortExecute;
  vtkIdType numSimple=0, numBEdges=0, numFixed=0, numFEdges=0;
  vtkPolyData *inMesh, *Mesh;
  vtkPoints *inPts;
  vtkCellArray *inVerts, *inLines, *inPolys;
  vtkPoints *newPts;
  vtkMeshVertexPtr Verts;
  vtkCellLocator *cellLocator=NULL;
  
    // Check input
    //
  numPts=input->GetNumberOfPoints();
  numCells=input->GetNumberOfCells();
  if (numPts < 1 || numCells < 1)
  {
    cout<<"No data to smooth!"<<endl;
    return 1;
  }
  
  CosFeatureAngle = 
    cos((double) vtkMath::DegreesToRadians() * this->FeatureAngle);
  CosEdgeAngle = cos((double) vtkMath::DegreesToRadians() * this->EdgeAngle);
  
  if(DebugLevel>5) {
    cout<<"Smoothing " << numPts << " vertices, " << numCells 
      << " cells with:\n"
      << "\tConvergence= " << this->Convergence << "\n"
      << "\tIterations= " << this->NumberOfIterations << "\n"
      << "\tRelaxation Factor= " << this->RelaxationFactor << "\n"
      << "\tEdge Angle= " << this->EdgeAngle << "\n"
      << "\tBoundary Smoothing " << (this->BoundarySmoothing ? "On\n" : "Off\n")
      << "\tFeature Edge Smoothing " << (this->FeatureEdgeSmoothing ? "On\n" : "Off\n")
      << "\tError Scalars " << (this->GenerateErrorScalars ? "On\n" : "Off\n")
      << "\tError Vectors " << (this->GenerateErrorVectors ? "On\n" : "Off\n")<<endl;
  }
    // Peform topological analysis. What we're gonna do is build a connectivity
    // array of connected vertices. The outcome will be one of three
    // classifications for a vertex: VTK_SIMPLE_VERTEX, VTK_FIXED_VERTEX. or
    // VTK_EDGE_VERTEX. Simple vertices are smoothed using all connected 
    // vertices. FIXED vertices are never smoothed. Edge vertices are smoothed
    // using a subset of the attached vertices.
    //
  if(DebugLevel>5) cout<<"===>Analyze topology==="<<endl;
  if(DebugLevel>5) cout<<"Analyzing topology..."<<endl;
  if(DebugLevel>5) cout<<"0:numPts="<<numPts<<endl;
  Verts = new vtkMeshVertex[numPts];
  for (i=0; i<numPts; i++)
  {
    if(DebugLevel>5) cout<<"0:VTK_SIMPLE_VERTEX"<<endl;
    Verts[i].type = VTK_SIMPLE_VERTEX; //can smooth
    Verts[i].edges = NULL;
  }
  
  inPts = input->GetPoints();
  conv = this->Convergence * input->GetLength();
  
  if(DebugLevel>5) cout<<"==polygons and triangle strips=="<<endl;
    // now polygons and triangle strips-------------------------------
  inPolys=input->GetPolys();
  numPolys = inPolys->GetNumberOfCells();
  
  if(DebugLevel>5) cout<<"numPolys="<<numPolys<<endl;
  
  if ( numPolys > 0 )
  { //build cell structure
    vtkCellArray *polys;
    vtkIdType cellId;
    int numNei, nei, edge;
    vtkIdType numNeiPts;
    vtkIdType *neiPts;
    double normal[3], neiNormal[3];
    vtkIdList *neighbors;
    
    neighbors = vtkIdList::New();
    neighbors->Allocate(VTK_CELL_SIZE);
    
    inMesh = vtkPolyData::New();
    inMesh->SetPoints(inPts);
    inMesh->SetPolys(inPolys);
    Mesh = inMesh;
    
    Mesh->BuildLinks(); //to do neighborhood searching
    polys = Mesh->GetPolys();
    
    for (cellId=0, polys->InitTraversal(); polys->GetNextCell(npts,pts); 
         cellId++)
    {
      if(DebugLevel>5) cout<<"->cellId="<<cellId<<endl;
      for (i=0; i < npts; i++) 
      {
        if(DebugLevel>5) cout<<"-->i="<<i<<endl;
        p1 = pts[i];
        p2 = pts[(i+1)%npts];
        
        if ( Verts[p1].edges == NULL )
        {
          Verts[p1].edges = vtkIdList::New();
          Verts[p1].edges->Allocate(16,6);
        }
        if ( Verts[p2].edges == NULL )
        {
          Verts[p2].edges = vtkIdList::New();
          Verts[p2].edges->Allocate(16,6);
        }
        
        Mesh->GetCellEdgeNeighbors(cellId,p1,p2,neighbors);
        numNei = neighbors->GetNumberOfIds();
        if(DebugLevel>5) cout<<"-->numNei="<<numNei<<endl;
        
        edge = VTK_SIMPLE_VERTEX;
        if ( numNei == 0 )
        {
          edge = VTK_BOUNDARY_EDGE_VERTEX;
        }
        
        else if ( numNei >= 2 )
        {
            // check to make sure that this edge hasn't been marked already
          for (j=0; j < numNei; j++)
          {
            if ( neighbors->GetId(j) < cellId )
            {
              break;
            }
          }
          if ( j >= numNei )
          {
            edge = VTK_FEATURE_EDGE_VERTEX;
          }
        }
        
        else if ( numNei == 1 && (nei=neighbors->GetId(0)) > cellId ) 
        {
          vtkPolygon::ComputeNormal(inPts,npts,pts,normal);
          Mesh->GetCellPoints(nei,numNeiPts,neiPts);
          vtkPolygon::ComputeNormal(inPts,numNeiPts,neiPts,neiNormal);
          
          if ( this->FeatureEdgeSmoothing &&
               vtkMath::Dot(normal,neiNormal) <= CosFeatureAngle ) 
          {
            edge = VTK_FEATURE_EDGE_VERTEX;
          }
        }
        else // a visited edge; skip rest of analysis
        {
          continue;
        }
        
        if ( edge && Verts[p1].type == VTK_SIMPLE_VERTEX )
        {
          Verts[p1].edges->Reset();
          Verts[p1].edges->InsertNextId(p2);
          Verts[p1].type = edge;
        }
        else if ( (edge && Verts[p1].type == VTK_BOUNDARY_EDGE_VERTEX) ||
                  (edge && Verts[p1].type == VTK_FEATURE_EDGE_VERTEX) ||
                  (!edge && Verts[p1].type == VTK_SIMPLE_VERTEX ) )
        {
          Verts[p1].edges->InsertNextId(p2);
          if ( Verts[p1].type && edge == VTK_BOUNDARY_EDGE_VERTEX )
          {
            Verts[p1].type = VTK_BOUNDARY_EDGE_VERTEX;
          }
        }
        
        if ( edge && Verts[p2].type == VTK_SIMPLE_VERTEX )
        {
          Verts[p2].edges->Reset();
          Verts[p2].edges->InsertNextId(p1);
          Verts[p2].type = edge;
        }
        else if ( (edge && Verts[p2].type == VTK_BOUNDARY_EDGE_VERTEX ) ||
                  (edge && Verts[p2].type == VTK_FEATURE_EDGE_VERTEX) ||
                  (!edge && Verts[p2].type == VTK_SIMPLE_VERTEX ) )
        {
          Verts[p2].edges->InsertNextId(p1);
          if ( Verts[p2].type && edge == VTK_BOUNDARY_EDGE_VERTEX )
          {
            Verts[p2].type = VTK_BOUNDARY_EDGE_VERTEX;
          }
        }
      }
    }
    
    inMesh->Delete();
    
    neighbors->Delete();
  }//if strips or polys
  
    //post-process edge vertices to make sure we can smooth them
  for (i=0; i<numPts; i++)
  {
    if ( Verts[i].type == VTK_SIMPLE_VERTEX )
    {
      numSimple++;
    }
    
    else if ( Verts[i].type == VTK_FIXED_VERTEX )
    {
      numFixed++;
    }
    
    else if ( Verts[i].type == VTK_FEATURE_EDGE_VERTEX ||
              Verts[i].type == VTK_BOUNDARY_EDGE_VERTEX )
    { //see how many edges; if two, what the angle is
      
      if ( !this->BoundarySmoothing && 
           Verts[i].type == VTK_BOUNDARY_EDGE_VERTEX )
      {
        if(DebugLevel>5) cout<<"Verts[i].type = VTK_FIXED_VERTEX; 4"<<endl;
        Verts[i].type = VTK_FIXED_VERTEX;
        numBEdges++;
      }
      
      else if ( (npts = Verts[i].edges->GetNumberOfIds()) != 2 )
      {
        if(DebugLevel>5) cout<<"Verts["<<i<<"].type = VTK_FIXED_VERTEX; 5"<<endl;
        Verts[i].type = VTK_FIXED_VERTEX;
        numFixed++;
      }
      
      else //check angle between edges
      {
        inPts->GetPoint(Verts[i].edges->GetId(0),x1);
        inPts->GetPoint(i,x2);
        inPts->GetPoint(Verts[i].edges->GetId(1),x3);
        
        for (k=0; k<3; k++)
        {
          l1[k] = x2[k] - x1[k];
          l2[k] = x3[k] - x2[k];
        }
        if ( vtkMath::Normalize(l1) >= 0.0 &&
             vtkMath::Normalize(l2) >= 0.0 &&
             vtkMath::Dot(l1,l2) < CosEdgeAngle)
        {
          if(DebugLevel>5) cout<<"Verts["<<i<<"].type = VTK_FIXED_VERTEX; 6"<<endl;
          Verts[i].type = VTK_FIXED_VERTEX;
          numFixed++;
        }
        else
        {
          if ( Verts[i].type == VTK_FEATURE_EDGE_VERTEX )
          {
            numFEdges++;
          }
          else
          {
            numBEdges++;
          }
        }
      }//if along edge
    }//if edge vertex
  }//for all points
  
  if(DebugLevel>5) {
    cout<<"Found\n\t" << numSimple << " simple vertices\n\t"
      << numFEdges << " feature edge vertices\n\t"
      << numBEdges << " boundary edge vertices\n\t"
      << numFixed << " fixed vertices\n\t"<<endl;
    cout<<"1:numPts="<<numPts<<endl;
  }
  
  for (i=0; i<numPts; i++) 
  {
    if(DebugLevel>5) cout<<"Verts["<<i<<"].type="<<VertexType2Str(Verts[i].type)<<endl;
    if(Verts[i].edges != NULL && (npts = Verts[i].edges->GetNumberOfIds()) > 0)
    {
      for (j=0; j<npts; j++)
      {
        if(DebugLevel>5) cout<<"Verts["<<i<<"].edges->GetId("<<j<<")="<<Verts[i].edges->GetId(j)<<endl;
      }//for all connected points
    }
  }
  
  //Copy node type info from Verts
  EG_VTKDCN(vtkCharArray, node_type, grid, "node_type");
  if(DebugLevel>5) cout<<"nodes.size()="<<nodes.size()<<endl;
  foreach(vtkIdType node,nodes)
  {
    if(DebugLevel>5) cout<<"Verts["<<node<<"].type="<<VertexType2Str(Verts[node].type)<<endl;
    char T=Verts[node].type;
    int N=N_neighbour_BCs(node);
    //TODO: There could be more cases. Either define new node types or create a node field containing the number of BCs.
    if(N>2) node_type->SetValue(node,BC_FIXED_VERTEX);
    else if(N==2){
      if(T==VTK_FIXED_VERTEX) node_type->SetValue(node,BC_FIXED_VERTEX);
      else if(T==VTK_BOUNDARY_EDGE_VERTEX) node_type->SetValue(node,BC_BOUNDARY_EDGE_VERTEX);
      else node_type->SetValue(node,BC_FEATURE_EDGE_VERTEX);
    }
    else node_type->SetValue(node,T);
  }
  
  //free up connectivity storage
  for (int i=0; i<numPts; i++)
  {
    if ( Verts[i].edges != NULL )
    {
      Verts[i].edges->Delete();
      Verts[i].edges = NULL;
    }
  }
  delete [] Verts;
  
  if(DebugLevel>0) cout<<"===UpdateNodeType_all END==="<<endl;
  return(0);
}
//End of UpdateNodeType_all

int Operation::UpdateNodeType()
{
  if(DebugLevel>47) cout<<"this->FeatureAngle="<<this->FeatureAngle<<endl;
  if(DebugLevel>47) cout<<"this->EdgeAngle="<<this->EdgeAngle<<endl;
  cout<<"===UpdateNodeType START==="<<endl;
  
  getAllSurfaceCells(cells,grid);
  if(DebugLevel>5) cout<<"cells.size()="<<cells.size()<<endl;
  
  EG_VTKSP(vtkPolyData, pdata);
  //   addToPolyData(m_SelectedCells, pdata, grid);
  addToPolyData(cells, pdata, grid);
  
  vtkPolyData* input=pdata;
  
  vtkPolyData *source = 0;
  
  vtkIdType numPts, numCells, i, numPolys;
  int j, k;
  vtkIdType npts = 0;
  vtkIdType *pts = 0;
  vtkIdType p1, p2;
  double x[3], y[3], deltaX[3], xNew[3], conv, maxDist, dist, factor;
  double x1[3], x2[3], x3[3], l1[3], l2[3];
  double CosFeatureAngle; //Cosine of angle between adjacent polys
  double CosEdgeAngle; // Cosine of angle between adjacent edges
  double closestPt[3], dist2, *w = NULL;
  int iterationNumber, abortExecute;
  vtkIdType numSimple=0, numBEdges=0, numFixed=0, numFEdges=0;
  vtkPolyData *inMesh, *Mesh;
  vtkPoints *inPts;
  vtkCellArray *inVerts, *inLines, *inPolys;
  vtkPoints *newPts;
  vtkMeshVertexPtr Verts;
  vtkCellLocator *cellLocator=NULL;
  
    // Check input
    //
  numPts=input->GetNumberOfPoints();
  numCells=input->GetNumberOfCells();
  if (numPts < 1 || numCells < 1)
  {
    cout<<"No data to smooth!"<<endl;
    return 1;
  }
  
  CosFeatureAngle = 
    cos((double) vtkMath::DegreesToRadians() * this->FeatureAngle);
  CosEdgeAngle = cos((double) vtkMath::DegreesToRadians() * this->EdgeAngle);
  
  if(DebugLevel>5) {
    cout<<"Smoothing " << numPts << " vertices, " << numCells 
      << " cells with:\n"
      << "\tConvergence= " << this->Convergence << "\n"
      << "\tIterations= " << this->NumberOfIterations << "\n"
      << "\tRelaxation Factor= " << this->RelaxationFactor << "\n"
      << "\tEdge Angle= " << this->EdgeAngle << "\n"
      << "\tBoundary Smoothing " << (this->BoundarySmoothing ? "On\n" : "Off\n")
      << "\tFeature Edge Smoothing " << (this->FeatureEdgeSmoothing ? "On\n" : "Off\n")
      << "\tError Scalars " << (this->GenerateErrorScalars ? "On\n" : "Off\n")
      << "\tError Vectors " << (this->GenerateErrorVectors ? "On\n" : "Off\n")<<endl;
  }
    // Peform topological analysis. What we're gonna do is build a connectivity
    // array of connected vertices. The outcome will be one of three
    // classifications for a vertex: VTK_SIMPLE_VERTEX, VTK_FIXED_VERTEX. or
    // VTK_EDGE_VERTEX. Simple vertices are smoothed using all connected 
    // vertices. FIXED vertices are never smoothed. Edge vertices are smoothed
    // using a subset of the attached vertices.
    //
  if(DebugLevel>5) cout<<"===>Analyze topology==="<<endl;
  if(DebugLevel>5) cout<<"Analyzing topology..."<<endl;
  if(DebugLevel>5) cout<<"0:numPts="<<numPts<<endl;
  Verts = new vtkMeshVertex[numPts];
  for (i=0; i<numPts; i++)
  {
    if(DebugLevel>5) cout<<"0:VTK_SIMPLE_VERTEX"<<endl;
    Verts[i].type = VTK_SIMPLE_VERTEX; //can smooth
    Verts[i].edges = NULL;
  }
  
  inPts = input->GetPoints();
  conv = this->Convergence * input->GetLength();
  
  if(DebugLevel>5) cout<<"==polygons and triangle strips=="<<endl;
    // now polygons and triangle strips-------------------------------
  inPolys=input->GetPolys();
  numPolys = inPolys->GetNumberOfCells();
  
  if(DebugLevel>5) cout<<"numPolys="<<numPolys<<endl;
  
  if ( numPolys > 0 )
  { //build cell structure
    vtkCellArray *polys;
    vtkIdType cellId;
    int numNei, nei, edge;
    vtkIdType numNeiPts;
    vtkIdType *neiPts;
    double normal[3], neiNormal[3];
    vtkIdList *neighbors;
    
    neighbors = vtkIdList::New();
    neighbors->Allocate(VTK_CELL_SIZE);
    
    inMesh = vtkPolyData::New();
    inMesh->SetPoints(inPts);
    inMesh->SetPolys(inPolys);
    Mesh = inMesh;
    
    Mesh->BuildLinks(); //to do neighborhood searching
    polys = Mesh->GetPolys();
    
    for (cellId=0, polys->InitTraversal(); polys->GetNextCell(npts,pts); 
         cellId++)
    {
      if(DebugLevel>5) cout<<"->cellId="<<cellId<<endl;
      for (i=0; i < npts; i++) 
      {
        if(DebugLevel>5) cout<<"-->i="<<i<<endl;
        p1 = pts[i];
        p2 = pts[(i+1)%npts];
        
        if ( Verts[p1].edges == NULL )
        {
          Verts[p1].edges = vtkIdList::New();
          Verts[p1].edges->Allocate(16,6);
        }
        if ( Verts[p2].edges == NULL )
        {
          Verts[p2].edges = vtkIdList::New();
          Verts[p2].edges->Allocate(16,6);
        }
        
        Mesh->GetCellEdgeNeighbors(cellId,p1,p2,neighbors);
        numNei = neighbors->GetNumberOfIds();
        if(DebugLevel>5) cout<<"-->numNei="<<numNei<<endl;
        
        edge = VTK_SIMPLE_VERTEX;
        if ( numNei == 0 )
        {
          edge = VTK_BOUNDARY_EDGE_VERTEX;
        }
        
        else if ( numNei >= 2 )
        {
            // check to make sure that this edge hasn't been marked already
          for (j=0; j < numNei; j++)
          {
            if ( neighbors->GetId(j) < cellId )
            {
              break;
            }
          }
          if ( j >= numNei )
          {
            edge = VTK_FEATURE_EDGE_VERTEX;
          }
        }
        
        else if ( numNei == 1 && (nei=neighbors->GetId(0)) > cellId ) 
        {
          vtkPolygon::ComputeNormal(inPts,npts,pts,normal);
          Mesh->GetCellPoints(nei,numNeiPts,neiPts);
          vtkPolygon::ComputeNormal(inPts,numNeiPts,neiPts,neiNormal);
          
          if ( this->FeatureEdgeSmoothing &&
               vtkMath::Dot(normal,neiNormal) <= CosFeatureAngle ) 
          {
            edge = VTK_FEATURE_EDGE_VERTEX;
          }
        }
        else // a visited edge; skip rest of analysis
        {
          continue;
        }
        
        if ( edge && Verts[p1].type == VTK_SIMPLE_VERTEX )
        {
          Verts[p1].edges->Reset();
          Verts[p1].edges->InsertNextId(p2);
          Verts[p1].type = edge;
        }
        else if ( (edge && Verts[p1].type == VTK_BOUNDARY_EDGE_VERTEX) ||
                  (edge && Verts[p1].type == VTK_FEATURE_EDGE_VERTEX) ||
                  (!edge && Verts[p1].type == VTK_SIMPLE_VERTEX ) )
        {
          Verts[p1].edges->InsertNextId(p2);
          if ( Verts[p1].type && edge == VTK_BOUNDARY_EDGE_VERTEX )
          {
            Verts[p1].type = VTK_BOUNDARY_EDGE_VERTEX;
          }
        }
        
        if ( edge && Verts[p2].type == VTK_SIMPLE_VERTEX )
        {
          Verts[p2].edges->Reset();
          Verts[p2].edges->InsertNextId(p1);
          Verts[p2].type = edge;
        }
        else if ( (edge && Verts[p2].type == VTK_BOUNDARY_EDGE_VERTEX ) ||
                  (edge && Verts[p2].type == VTK_FEATURE_EDGE_VERTEX) ||
                  (!edge && Verts[p2].type == VTK_SIMPLE_VERTEX ) )
        {
          Verts[p2].edges->InsertNextId(p1);
          if ( Verts[p2].type && edge == VTK_BOUNDARY_EDGE_VERTEX )
          {
            Verts[p2].type = VTK_BOUNDARY_EDGE_VERTEX;
          }
        }
      }
    }
    
    inMesh->Delete();
    
    neighbors->Delete();
  }//if strips or polys
  
    //post-process edge vertices to make sure we can smooth them
  for (i=0; i<numPts; i++)
  {
    if ( Verts[i].type == VTK_SIMPLE_VERTEX )
    {
      numSimple++;
    }
    
    else if ( Verts[i].type == VTK_FIXED_VERTEX )
    {
      numFixed++;
    }
    
    else if ( Verts[i].type == VTK_FEATURE_EDGE_VERTEX ||
              Verts[i].type == VTK_BOUNDARY_EDGE_VERTEX )
    { //see how many edges; if two, what the angle is
      
      if ( !this->BoundarySmoothing && 
           Verts[i].type == VTK_BOUNDARY_EDGE_VERTEX )
      {
        if(DebugLevel>5) cout<<"Verts[i].type = VTK_FIXED_VERTEX; 4"<<endl;
        Verts[i].type = VTK_FIXED_VERTEX;
        numBEdges++;
      }
      
      else if ( (npts = Verts[i].edges->GetNumberOfIds()) != 2 )
      {
        if(DebugLevel>5) cout<<"Verts["<<i<<"].type = VTK_FIXED_VERTEX; 5"<<endl;
        Verts[i].type = VTK_FIXED_VERTEX;
        numFixed++;
      }
      
      else //check angle between edges
      {
        inPts->GetPoint(Verts[i].edges->GetId(0),x1);
        inPts->GetPoint(i,x2);
        inPts->GetPoint(Verts[i].edges->GetId(1),x3);
        
        for (k=0; k<3; k++)
        {
          l1[k] = x2[k] - x1[k];
          l2[k] = x3[k] - x2[k];
        }
        if ( vtkMath::Normalize(l1) >= 0.0 &&
             vtkMath::Normalize(l2) >= 0.0 &&
             vtkMath::Dot(l1,l2) < CosEdgeAngle)
        {
          if(DebugLevel>5) cout<<"Verts["<<i<<"].type = VTK_FIXED_VERTEX; 6"<<endl;
          Verts[i].type = VTK_FIXED_VERTEX;
          numFixed++;
        }
        else
        {
          if ( Verts[i].type == VTK_FEATURE_EDGE_VERTEX )
          {
            numFEdges++;
          }
          else
          {
            numBEdges++;
          }
        }
      }//if along edge
    }//if edge vertex
  }//for all points
  
  if(DebugLevel>5) {
    cout<<"Found\n\t" << numSimple << " simple vertices\n\t"
      << numFEdges << " feature edge vertices\n\t"
      << numBEdges << " boundary edge vertices\n\t"
      << numFixed << " fixed vertices\n\t"<<endl;
    cout<<"1:numPts="<<numPts<<endl;
  }
  
  for (i=0; i<numPts; i++) 
  {
    if(DebugLevel>5) cout<<"Verts["<<i<<"].type="<<VertexType2Str(Verts[i].type)<<endl;
    if(Verts[i].edges != NULL && (npts = Verts[i].edges->GetNumberOfIds()) > 0)
    {
      for (j=0; j<npts; j++)
      {
        if(DebugLevel>5) cout<<"Verts["<<i<<"].edges->GetId("<<j<<")="<<Verts[i].edges->GetId(j)<<endl;
      }//for all connected points
    }
  }
  
  //Copy node type info from Verts
  EG_VTKDCN(vtkCharArray, node_type, grid, "node_type");
  if(DebugLevel>5) cout<<"nodes.size()="<<nodes.size()<<endl;
  foreach(vtkIdType node,nodes)
  {
    if(DebugLevel>5) cout<<"Verts["<<node<<"].type="<<VertexType2Str(Verts[node].type)<<endl;
    node_type->SetValue(node,Verts[node].type);
  }
  
  //free up connectivity storage
  for (int i=0; i<numPts; i++)
  {
    if ( Verts[i].edges != NULL )
    {
      Verts[i].edges->Delete();
      Verts[i].edges = NULL;
    }
  }
  delete [] Verts;
  
  cout<<"===UpdateNodeType END==="<<endl;
  return(0);
}
//End of UpdateNodeType

// DEFINITIONS:
// Normal cell: nothing has changed
// Dead cell: the cell does not exist anymore
// Mutated cell: the cell's form has changed
// Mutilated cell: the cell has less points than before

vtkIdType Operation::FindSnapPoint(vtkUnstructuredGrid *src, vtkIdType DeadNode,QSet <vtkIdType> & DeadCells,QSet <vtkIdType> & MutatedCells,QSet <vtkIdType> & MutilatedCells, int& N_newpoints, int& N_newcells)
{
  //TODO: Organize cases and make sure all are considered if possible. It's the final countdown!!!
  getAllSurfaceCells(cells,src);
  getNodesFromCells(cells, nodes, src);
  setGrid(src);
  setCells(cells);
  
  UpdateNodeType_all();
  
  setDebugLevel(20);
  
  EG_VTKDCN(vtkCharArray, node_type, src, "node_type");
  if(node_type->GetValue(DeadNode)==VTK_FIXED_VERTEX)
  {
    cout<<"Sorry, unable to remove fixed vertex."<<endl;
    return(-1);
  }
  
    //src grid info
  int N_points=src->GetNumberOfPoints();
  int N_cells=src->GetNumberOfCells();
  N_newpoints=-1;
  N_newcells=0;
  
  vtkIdType SnapPoint=-1;
  
  foreach(vtkIdType PSP, n2n[DeadNode])
  {
    bool IsValidSnapPoint=true;
    
    if(DebugLevel>10) cout<<"====>PSP="<<PSP<<endl;
    bool IsTetra=true;
    if(NumberOfCommonPoints(DeadNode,PSP,IsTetra)>2)//common point check
    {
      if(DebugLevel>10) cout<<"Sorry, but you are not allowed to move point "<<DeadNode<<" to point "<<PSP<<"."<<endl;
      IsValidSnapPoint=false;
    }
    if(IsTetra)//tetra check
    {
      if(DebugLevel>10) cout<<"Sorry, but you are not allowed to move point "<<DeadNode<<" to point "<<PSP<<"."<<endl;
      IsValidSnapPoint=false;
    }
    
    //count number of points and cells to remove + analyse cell transformations
    N_newpoints=-1;
    N_newcells=0;
    DeadCells.clear();
    MutatedCells.clear();
    MutilatedCells.clear();
    foreach(vtkIdType C, n2c[DeadNode])//loop through potentially dead cells
    {
      //get points around cell
      vtkIdType N_pts, *pts;
      src->GetCellPoints(C, N_pts, pts);
      
      bool ContainsSnapPoint=false;
      bool invincible=false;
      for(int i=0;i<N_pts;i++)
      {
        if(DebugLevel>10) cout<<"pts["<<i<<"]="<<pts[i]<<" and PSP="<<PSP<<endl;
        if(pts[i]==PSP) {ContainsSnapPoint=true;}
        if(pts[i]!=DeadNode && pts[i]!=PSP &&  n2c[pts[i]].size()<=1) invincible=true;
      }
      if(ContainsSnapPoint)
      {
        if(N_pts==3)//dead cell
        {
          if(invincible)//Check that empty lines aren't left behind when a cell is killed
          {
            if(DebugLevel>10) cout<<"Sorry, but you are not allowed to move point "<<DeadNode<<" to point "<<PSP<<"."<<endl;
            IsValidSnapPoint=false;
          }
          else
          {
            DeadCells.insert(C);
            N_newcells-=1;
            if(DebugLevel>10) cout<<"cell "<<C<<" has been pwned!"<<endl;
          }
        }
        else
        {
          cout<<"RED ALERT: Xenomorph detected!"<<endl;
          EG_BUG;
        }
      }
      else
      {
        vtkIdType src_N_pts, *src_pts;
        src->GetCellPoints(C, src_N_pts, src_pts);
        
        if(src_N_pts!=3)
        {
          cout<<"RED ALERT: Xenomorph detected!"<<endl;
          EG_BUG;
        }
        
        vtkIdType OldTriangle[3];
        vtkIdType NewTriangle[3];
        
        for(int i=0;i<src_N_pts;i++)
        {
          OldTriangle[i]=src_pts[i];
          NewTriangle[i]=( (src_pts[i]==DeadNode) ? PSP : src_pts[i] );
        }
        vec3_t Old_N= triNormal(src, OldTriangle[0], OldTriangle[1], OldTriangle[2]);
        vec3_t New_N= triNormal(src, NewTriangle[0], NewTriangle[1], NewTriangle[2]);
        double OldArea=Old_N.abs();
        double NewArea=New_N.abs();
        double scal=Old_N*New_N;
        double cross=(Old_N.cross(New_N)).abs();//double-cross on Nar Shadaa B-)
        
        if(DebugLevel>10) {
          cout<<"OldArea="<<OldArea<<endl;
          cout<<"NewArea="<<NewArea<<endl;
          cout<<"scal="<<scal<<endl;
          cout<<"cross="<<cross<<endl;
        }
        
        if(Old_N*New_N<Old_N*Old_N*1./100.)//area + inversion check
        {
          if(DebugLevel>10) cout<<"Sorry, but you are not allowed to move point "<<DeadNode<<" to point "<<PSP<<"."<<endl;
          IsValidSnapPoint=false;
        }
        
        //mutated cell
        MutatedCells.insert(C);
        if(DebugLevel>10) cout<<"cell "<<C<<" has been infected!"<<endl;
      }
    }
    
    if(N_cells+N_newcells<=0)//survivor check
    {
      if(DebugLevel>10) cout<<"Sorry, but you are not allowed to move point "<<DeadNode<<" to point "<<PSP<<"."<<endl;
      IsValidSnapPoint=false;
    }
    
    if(node_type->GetValue(DeadNode)==VTK_BOUNDARY_EDGE_VERTEX && node_type->GetValue(PSP)==VTK_SIMPLE_VERTEX)
    {
      if(DebugLevel>10) cout<<"Sorry, but you are not allowed to move point "<<DeadNode<<" to point "<<PSP<<"."<<endl;
      IsValidSnapPoint=false;
    }
    
    if(node_type->GetValue(DeadNode)==VTK_BOUNDARY_EDGE_VERTEX)
    {
      int BC=0;
      QVector <vtkIdType> Peons;
      getNeighbours(DeadNode, Peons, BC);
      if(!Peons.contains(PSP))
      {
        if(DebugLevel>0) cout<<"Sorry, but you are not allowed to move point "<<DeadNode<<" to point "<<PSP<<"."<<endl;
        IsValidSnapPoint=false;
      }
    }
    
    if(node_type->GetValue(DeadNode)==VTK_FEATURE_EDGE_VERTEX && node_type->GetValue(PSP)==VTK_SIMPLE_VERTEX)
    {
      if(DebugLevel>10) cout<<"Sorry, but you are not allowed to move point "<<DeadNode<<" to point "<<PSP<<"."<<endl;
      IsValidSnapPoint=false;
    }
    
    if(node_type->GetValue(DeadNode)==BC_FEATURE_EDGE_VERTEX && node_type->GetValue(PSP)==VTK_SIMPLE_VERTEX)
    {
      if(DebugLevel>10) cout<<"Sorry, but you are not allowed to move point "<<DeadNode<<" to point "<<PSP<<"."<<endl;
      IsValidSnapPoint=false;
    }
    
    if(node_type->GetValue(DeadNode)==VTK_FEATURE_EDGE_VERTEX)
    {
      int BC=0;
      QVector <vtkIdType> Peons;
      getNeighbours(DeadNode, Peons, BC);
      if(!Peons.contains(PSP))
      {
        if(DebugLevel>0) cout<<"Sorry, but you are not allowed to move point "<<DeadNode<<" to point "<<PSP<<"."<<endl;
        IsValidSnapPoint=false;
      }
    }
    
    //TODO: merge with previous case if possible
    if(node_type->GetValue(DeadNode)==BC_FEATURE_EDGE_VERTEX)
    {
      int BC=0;
      QVector <vtkIdType> Peons;
      getNeighbours_BC(DeadNode, Peons);
      if(!Peons.contains(PSP))
      {
        if(DebugLevel>0) cout<<"Sorry, but you are not allowed to move point "<<DeadNode<<" to point "<<PSP<<"."<<endl;
        IsValidSnapPoint=false;
      }
    }
    
    if(IsValidSnapPoint) {SnapPoint=PSP; break;}
  }//end of loop through potential SnapPoints
  
  if(DebugLevel>10)
  {
    cout<<"AT FINDSNAPPOINT EXIT"<<endl;
    cout<<"N_points="<<N_points<<endl;
    cout<<"N_cells="<<N_cells<<endl;
    cout<<"N_newpoints="<<N_newpoints<<endl;
    cout<<"N_newcells="<<N_newcells<<endl;
  }
  cout<<"MutilatedCells.size()="<<MutilatedCells.size()<<endl;
  cout<<"MutatedCells.size()="<<MutatedCells.size()<<endl;
  cout<<"DeadCells.size()="<<DeadCells.size()<<endl;
  return(SnapPoint);
  
  setDebugLevel(0);
}
//End of FindSnapPoint

bool Operation::DeletePoint(vtkUnstructuredGrid *src, vtkIdType DeadNode, int& N_newpoints, int& N_newcells)
{
  QSet <vtkIdType> DeadNodes;
  DeadNodes.insert(DeadNode);
  bool ret = DeleteSetOfPoints(src,DeadNodes, N_newpoints, N_newcells);
  return(ret);
}
//End of DeletePoint

bool Operation::DeleteSetOfPoints(vtkUnstructuredGrid *src, QSet <vtkIdType> DeadNodes, int& N_newpoints, int& N_newcells)
{
  QVector <vtkIdType> DeadNode_vector=Set2Vector(DeadNodes,false);
  
  getAllSurfaceCells(cells,src);
//   getNodesFromCells(cells, nodes, src);
  setGrid(src);
  setCells(cells);
  UpdateNodeType_all();
  
  //src grid info
  int N_points=src->GetNumberOfPoints();
  int N_cells=src->GetNumberOfCells();
  
  QSet <vtkIdType> DeadCells;
  QSet <vtkIdType> MutatedCells;
  QSet <vtkIdType> MutilatedCells;
  QVector <vtkIdType> SnapPoint(DeadNode_vector.size());
  
    //counter init
  N_newpoints=0;
  N_newcells=0;
  
  for(int i=0;i<DeadNode_vector.size();i++)
  {
    if(DeadNode_vector[i]<0 || DeadNode_vector[i]>=N_points)
    {
      cout<<"Warning: Point out of range: DeadNode_vector[i]="<<DeadNode_vector[i]<<" N_points="<<N_points<<endl;
      return(false);
    }
    
    if(DebugLevel>10) {
      cout<<"BEFORE FINDSNAPPOINT"<<endl;
      cout<<"N_points="<<N_points<<endl;
      cout<<"N_cells="<<N_cells<<endl;
      cout<<"N_newpoints="<<N_newpoints<<endl;
      cout<<"N_newcells="<<N_newcells<<endl;
    }
    
    //local values
    int l_N_newpoints;
    int l_N_newcells;
    QSet <vtkIdType> l_DeadCells;
    QSet <vtkIdType> l_MutatedCells;
    QSet <vtkIdType> l_MutilatedCells;
    
    SnapPoint[i]=FindSnapPoint(src,DeadNode_vector[i], l_DeadCells, l_MutatedCells, l_MutilatedCells, l_N_newpoints, l_N_newcells);
    //global values
    N_newpoints+=l_N_newpoints;
    N_newcells+=l_N_newcells;
    DeadCells.unite(l_DeadCells);//DeadCells unite! Kill the living! :D
    MutatedCells.unite(l_MutatedCells);
    MutilatedCells.unite(l_MutilatedCells);
    
    if(DebugLevel>0) cout<<"===>DeadNode_vector[i]="<<DeadNode_vector[i]<<" moving to SNAPPOINT="<<SnapPoint[i]<<" DebugLevel="<<DebugLevel<<endl;
    if(SnapPoint[i]<0) {cout<<"Sorry no possible SnapPoint found."<<endl; return(false);}
    
  }
  //allocate
  if(DebugLevel>10) {
    cout<<"BEFORE ALLOCATION"<<endl;
    cout<<"N_points="<<N_points<<endl;
    cout<<"N_cells="<<N_cells<<endl;
    cout<<"N_newpoints="<<N_newpoints<<endl;
    cout<<"N_newcells="<<N_newcells<<endl;
  }
//   N_points=src->GetNumberOfPoints();
//   N_cells=src->GetNumberOfCells();
  
  if(DebugLevel>47) {
    cout<<"N_points="<<N_points<<endl;
    cout<<"N_cells="<<N_cells<<endl;
    cout<<"N_newpoints="<<N_newpoints<<endl;
    cout<<"N_newcells="<<N_newcells<<endl;
  }
  EG_VTKSP(vtkUnstructuredGrid,dst);
  allocateGrid(dst,N_cells+N_newcells,N_points+N_newpoints);
  
  //vector used to redefine the new point IDs
  QVector <vtkIdType> OffSet(N_points);
  
  //copy undead points
  vtkIdType dst_id_node=0;
  for (vtkIdType src_id_node = 0; src_id_node < N_points; src_id_node++) {//loop through src points
    if(!DeadNode_vector.contains(src_id_node))//if the node isn't dead, copy it
    {
      vec3_t x;
      src->GetPoints()->GetPoint(src_id_node, x.data());
      dst->GetPoints()->SetPoint(dst_id_node, x.data());
      copyNodeData(src, src_id_node, dst, dst_id_node);
      OffSet[src_id_node]=src_id_node-dst_id_node;
      dst_id_node++;
    }
    else
    {
      if(DebugLevel>0) cout<<"src_id_node="<<src_id_node<<" dst_id_node="<<dst_id_node<<endl;
    }
  };
  if(DebugLevel>10) {
    cout<<"DeadCells="<<DeadCells<<endl;
    cout<<"MutatedCells="<<MutatedCells<<endl;
    cout<<"MutilatedCells="<<MutilatedCells<<endl;
  }
  //Copy undead cells
  for (vtkIdType id_cell = 0; id_cell < src->GetNumberOfCells(); ++id_cell) {//loop through src cells
    if(!DeadCells.contains(id_cell))//if the cell isn't dead
    {
      vtkIdType src_N_pts, *src_pts;
      vtkIdType dst_N_pts, *dst_pts;
      src->GetCellPoints(id_cell, src_N_pts, src_pts);
      
      vtkIdType type_cell = src->GetCellType(id_cell);
      if(DebugLevel>10) cout<<"-->id_cell="<<id_cell<<endl;
      if(DebugLevel>10) for(int i=0;i<src_N_pts;i++) cout<<"src_pts["<<i<<"]="<<src_pts[i]<<endl;
//       src->GetCellPoints(id_cell, dst_N_pts, dst_pts);
      dst_N_pts=src_N_pts;
      dst_pts=new vtkIdType[dst_N_pts];
      if(MutatedCells.contains(id_cell))//mutated cell
      {
        if(DebugLevel>10) cout<<"processing mutated cell "<<id_cell<<endl;
        for(int i=0;i<src_N_pts;i++)
        {
          int DeadIndex = DeadNode_vector.indexOf(src_pts[i]);
          if(DeadIndex!=-1) {
            if(DebugLevel>10) {
              cout<<"SnapPoint="<<SnapPoint[DeadIndex]<<endl;
              cout<<"OffSet[SnapPoint]="<<OffSet[SnapPoint[DeadIndex]]<<endl;
              cout<<"src_pts["<<i<<"]="<<src_pts[i]<<endl;
            }
            dst_pts[i]=SnapPoint[DeadIndex]-OffSet[SnapPoint[DeadIndex]];
          }
          else dst_pts[i]=src_pts[i]-OffSet[src_pts[i]];
        }
        if(DebugLevel>10) cout<<"--->dst_pts:"<<endl;
        if(DebugLevel>10) for(int i=0;i<dst_N_pts;i++) cout<<"dst_pts["<<i<<"]="<<dst_pts[i]<<endl;
        
      }
      else if(MutilatedCells.contains(id_cell))//mutilated cell (ex: square becoming triangle) (WARNING: Not fully functional yet)
      {
        cout<<"FATAL ERROR: Quads not supported yet."<<endl;EG_BUG;
        
        if(DebugLevel>10) cout<<"processing mutilated cell "<<id_cell<<endl;
        
        if(type_cell==VTK_QUAD) {
          type_cell=VTK_TRIANGLE;
          dst_N_pts-=1;
        }
        else {cout<<"FATAL ERROR: Unknown mutilated cell detected! It is not a quad! Potential xenomorph infestation!"<<endl;EG_BUG;}
        //merge points
        int j=0;
        for(int i=0;i<src_N_pts;i++)
        {
/*          if(src_pts[i]==SnapPoint) { dst_pts[j]=SnapPoint-OffSet[SnapPoint];j++; }//SnapPoint
          else if(src_pts[i]!=DeadNode_vector[i]) { dst_pts[j]=src_pts[i]-OffSet[src_pts[i]];j++; }//pre-snap/dead + post-snap/dead*/
          //do nothing in case of DeadNode_vector[i]
        }
      }
      else//normal cell
      {
        if(DebugLevel>10) cout<<"processing normal cell "<<id_cell<<endl;
        for(int i=0;i<src_N_pts;i++)
        {
          dst_pts[i]=src_pts[i]-OffSet[src_pts[i]];
        }
      }
      //copy the cell
      vtkIdType id_new_cell = dst->InsertNextCell(type_cell, dst_N_pts, dst_pts);
      copyCellData(src, id_cell, dst, id_new_cell);
      if(DebugLevel>10) {
        cout<<"===Copying cell "<<id_cell<<" to "<<id_new_cell<<"==="<<endl;
        cout<<"src_pts:"<<endl;
        for(int i=0;i<src_N_pts;i++) cout<<"src_pts["<<i<<"]="<<src_pts[i]<<endl;
        cout<<"dst_pts:"<<endl;
        for(int i=0;i<dst_N_pts;i++) cout<<"dst_pts["<<i<<"]="<<dst_pts[i]<<endl;
        cout<<"OffSet="<<OffSet<<endl;
        cout<<"===Copying cell end==="<<endl;
      }
      delete dst_pts;
    }
  };
  
//   cout_grid(cout,dst,true,true,true,true);
  makeCopy(dst, src);
  return(true);
}
//End of DeleteSetOfPoints

void Operation::TxtSave(QString a_filename)
{
  cout << a_filename.toAscii().data() << endl;
  ofstream file;
  file.open(a_filename.toAscii().data());
  cout_grid(file,grid,true,true,true,true);
  file.close();
}

void Operation::DualSave(QString a_filename)
{
  TxtSave(a_filename+".txt");
  GuiMainWindow::pointer()->QuickSave(a_filename+".vtu");
}

//---------------------------------------------------
//Utility functions used in Roland's formulas
//Should be renamed to be more explicit
//Some could be moved into geometrytools
//Some are pretty useless

///perimeter
double Operation::Um(vtkIdType D) {
  double ret=0;
  vtkIdType N_pts, *pts;
  grid->GetCellPoints(D, N_pts, pts);
  for(int i=0;i<N_pts;i++)
  {
    vec3_t A,B;
    grid->GetPoints()->GetPoint(pts[i], A.data());
    grid->GetPoints()->GetPoint(pts[(i+1)%N_pts], B.data());
    ret+=(B-A).abs();
  }
  return(ret);
}

/// area of the circumscribed circle of the triangle
double Operation::A_U(vtkIdType D) {
  vtkIdType N_pts, *pts;
  grid->GetCellPoints(D, N_pts, pts);
  vec3_t A,B,C;
  grid->GetPoints()->GetPoint(pts[0], A.data());
  grid->GetPoints()->GetPoint(pts[1], B.data());
  grid->GetPoints()->GetPoint(pts[2], C.data());
  double a=(C-B).abs();
  double alpha=angle((B-A),(C-A));
  double R=a/(2*sin(alpha));
  return(M_PI*R*R);
}

/// triangle area
double Operation::A_D(vtkIdType D) {
  return(cellVA(grid,D));
}

/// triangle neighbours
double Operation::DN(int i,vtkIdType D) {
  return(c2c[D][i]);
}

/// number of edges
double Operation::nk(vtkIdType P) {
  return(n2n[P].size());
}

double Operation::G_k(vtkIdType node) {
  EG_VTKDCN(vtkDoubleArray, node_meshdensity, grid, "node_meshdensity");
  return(1.0/node_meshdensity->GetValue(node));
}

/// triangle nodes
double Operation::DK(int i,vtkIdType D) {
  vtkIdType N_pts, *pts;
  grid->GetCellPoints(D, N_pts, pts);
  return(pts[i]);
}

vtkIdType Operation::KK(int i,vtkIdType j,vtkIdType K) {//i=1 or 2, j=node2, K=node1
  if(i==1) return(K);
  else return(j);
}

double Operation::L_k(vtkIdType j,vtkIdType K)// node1 K, node2 j
{
  vec3_t A;
  vec3_t B;
  grid->GetPoints()->GetPoint(K, A.data());
  grid->GetPoints()->GetPoint(j, B.data());
  return((B-A).abs());
}

double Operation::Q_L(vtkIdType D)
{
      // Um(D)/sum(G_k(DK(i,D)),i,1,3)
  double denom_sum=0;
  for(int i=0;i<3;i++)
  {
    denom_sum += G_k(DK(i,D));
  }
      /*if(DebugLevel>0) cout<<"D="<<D<<" Um(D)="<<Um(D)<<" denom_sum="<<denom_sum<<endl;*/
  return(Um(D)/denom_sum);
}

double Operation::Q_L1(vtkIdType P)
{
      // [2*sum(L_k(i~),i,1,nk(P))]/[sum(G_k(KK(1,i~))+G_k(KK(2,i~)),i,1,nk(P))]
  double num_sum=0;
  double denom_sum=0;
  foreach(vtkIdType j,n2n[P])
  {
    num_sum += 2*L_k(j,P);
    denom_sum += G_k(KK(1,j,P))+G_k(KK(2,j,P));
  }
  return(num_sum/denom_sum);
}

double Operation::Q_L2(vtkIdType P)
{
      // min([2*L_k(i~)]/[G_k(KK(1,i~))+G_k(KK(2,i~))])
  QVector <double> V;
  double num,denom;
  foreach(vtkIdType j,n2n[P])
  {
    num = 2*L_k(j,P);
    denom = G_k(KK(1,j,P))+G_k(KK(2,j,P));
    V.push_back(num/denom);
  }
  qSort(V.begin(),V.end());
  return(V[0]);
}

double Operation::T_min(int w)
{
      // sum([A_U(i)]/[A_D(i)^w]*[G_k(i)^(2*(w-1))],i,1,Nd)
  int N_cells=grid->GetNumberOfCells();
  double T=0;
  for(int i=0;i<N_cells;i++)
  {
    T += A_U(i)/pow(A_D(i),w)*pow(G_k(i),2*(w-1));
  }
  return(T);
}
//---------------------------------------------------
//These functions are not optimized. Avoid using them when possible.
QSet<vtkIdType> Operation::n2c_func(vtkIdType idx)
{
  QSet<int> tmp = n2c[_nodes[idx]];
  
  QSet<vtkIdType> ret;
  foreach(int i,tmp){
    if(i!=-1) ret.insert(cells[i]);
  }
  return(ret);
}

QSet<vtkIdType> Operation::n2n_func(vtkIdType idx)
{
  QSet<int> tmp = n2n[_nodes[idx]];

  QSet<vtkIdType> ret;
  foreach(int i,tmp){
    if(i!=-1) ret.insert(nodes[i]);
  }
  return(ret);
}

QVector<vtkIdType> Operation::c2c_func(vtkIdType idx)
{
  QVector<int> tmp = c2c[_cells[idx]];

  QVector<vtkIdType> ret;
  foreach(int i,tmp){
    if(i!=-1) ret.push_back(cells[i]);
  }
  return(ret);
}

int Operation::N_neighbour_BCs(vtkIdType a_node)
{
  EG_VTKDCC(vtkIntArray, cell_code, grid, "cell_code");
  QSet <vtkIdType> neighbours=n2c_func(a_node);
  QSet <int> bc;
  foreach(vtkIdType C, neighbours)
  {
    bc.insert(cell_code->GetValue(C));
  }
  return(bc.size());
}
