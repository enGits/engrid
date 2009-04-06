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

#include "geometrytools.h"
using namespace GeometryTools;

#include <QApplication>

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
  err = NULL;
  autoset = true;
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
};

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
};

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
    operate();
  };
  /*
  checkGrid();
  operate();
  */
};

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

stencil_t Operation::getStencil(vtkIdType id_cell1, int j1)
{
  stencil_t S;
  S.valid = true;
  S.id_cell1 = id_cell1;
  if (c2c[_cells[id_cell1]][j1] != -1) {
    S.id_cell2 = cells[c2c[_cells[id_cell1]][j1]];
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
double CurrentVertexAvgDist(vtkIdType a_vertex,QVector< QSet< int > > &n2n,vtkUnstructuredGrid *a_grid)
{
  double total_dist=0;
  double avg_dist=0;
  int N=n2n[a_vertex].size();
  vec3_t C;
  a_grid->GetPoint(a_vertex, C.data());
  foreach(int i,n2n[a_vertex])
  {
    vec3_t M;
    a_grid->GetPoint(i, M.data());
    total_dist+=(M-C).abs();
  }
  avg_dist=total_dist/(double)N;
  return(avg_dist);
}

double CurrentMeshDensity(vtkIdType a_vertex,QVector< QSet< int > > &n2n,vtkUnstructuredGrid *a_grid)
{
  double total_dist=0;
  double avg_dist=0;
  int N=n2n[a_vertex].size();
  vec3_t C;
  a_grid->GetPoint(a_vertex, C.data());
  foreach(int i,n2n[a_vertex])
  {
    vec3_t M;
    a_grid->GetPoint(i, M.data());
    total_dist+=(M-C).abs();
  }
  avg_dist=total_dist/(double)N;
  double avg_density=1./avg_dist;
  return(avg_density);
}

double DesiredVertexAvgDist(vtkIdType a_vertex,QVector< QSet< int > > &n2n,vtkUnstructuredGrid *a_grid)
{
  double total_dist=0;
  double avg_dist=0;
  EG_VTKDCN(vtkDoubleArray, node_meshdensity, a_grid, "node_meshdensity");
  int N=n2n[a_vertex].size();
  foreach(int i,n2n[a_vertex])
  {
    total_dist+=1./node_meshdensity->GetValue(i);
  }
  avg_dist=total_dist/(double)N;
  return(avg_dist);
}

double DesiredMeshDensity(vtkIdType a_vertex,QVector< QSet< int > > &n2n,vtkUnstructuredGrid *a_grid)
{
  double total_density=0;
  double avg_density=0;
  EG_VTKDCN(vtkDoubleArray, node_meshdensity, a_grid, "node_meshdensity");
  int N=n2n[a_vertex].size();
  foreach(int i,n2n[a_vertex])
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

// vtkIdType Operation::FindSnapPoint(vtkUnstructuredGrid *src, vtkIdType DeadNode)
// {
//   return(0);
// }

bool Operation::DeletePoint(vtkUnstructuredGrid *src, vtkIdType DeadNode)
{
  EG_VTKSP(vtkUnstructuredGrid, dst);
  
    //src grid info
  int N_points=src->GetNumberOfPoints();
  int N_cells=src->GetNumberOfCells();
  int N_newpoints=-1;
  int N_newcells=0;
  
  if(DeadNode<0 || DeadNode>=N_points)
  {
    cout<<"Warning: Point out of range: DeadNode="<<DeadNode<<" N_points="<<N_points<<endl;
    return(false);
  }
  
  QSet <vtkIdType> DeadCells;
  QSet <vtkIdType> MutatedCells;
  QSet <vtkIdType> MutilatedCells;
  
  vtkIdType SnapPoint=-1;
  //Find closest point to DeadNode
//   vtkIdType SnapPoint = getClosestNode(DeadNode,src);//DeadNode moves to SnapPoint
  
  foreach(vtkIdType PSP, n2n[DeadNode])
  {
    bool IsValidSnapPoint=true;
    
    cout<<"====>PSP="<<PSP<<endl;
    bool IsTetra=true;
    if(NumberOfCommonPoints(DeadNode,PSP,IsTetra)>2)//common point check
    {
      cout<<"Sorry, but you are not allowed to move point "<<DeadNode<<" to point "<<PSP<<"."<<endl;
      IsValidSnapPoint=false;
    }
    if(IsTetra)//tetra check
    {
      cout<<"Sorry, but you are not allowed to move point "<<DeadNode<<" to point "<<PSP<<"."<<endl;
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
      cout<<"C="<<C<<endl;
      //get points around cell
      vtkIdType N_pts, *pts;
      src->GetCellPoints(C, N_pts, pts);
      
      bool ContainsSnapPoint=false;
      bool invincible=false;
      for(int i=0;i<N_pts;i++)
      {
        cout<<"pts["<<i<<"]="<<pts[i]<<" and PSP="<<PSP<<endl;
        if(pts[i]==PSP) {ContainsSnapPoint=true;}
        if(pts[i]!=DeadNode && pts[i]!=PSP &&  n2c[pts[i]].size()<=1) invincible=true;
      }
      if(ContainsSnapPoint)
      {
        if(N_pts==3)//dead cell
        {
          if(invincible)//Check that empty lines aren't left behind when a cell is killed
          {
            cout<<"Sorry, but you are not allowed to move point "<<DeadNode<<" to point "<<PSP<<"."<<endl;
            IsValidSnapPoint=false;
          }
          else
          {
            DeadCells.insert(C);
            N_newcells-=1;
            cout<<"cell "<<C<<" has been pwned!"<<endl;
          }
        }
  /*      else if(N_pts==4)//mutilated cell
        {
          MutilatedCells.insert(C);
          cout<<"cell "<<C<<" has lost a limb!"<<endl;
        }*/
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
        
        cout<<"OldArea="<<OldArea<<endl;
        cout<<"NewArea="<<NewArea<<endl;
        cout<<"scal="<<scal<<endl;
        cout<<"cross="<<cross<<endl;
        
        if(Old_N*New_N<Old_N*Old_N*1./100.)//area + inversion check
        {
          cout<<"Sorry, but you are not allowed to move point "<<DeadNode<<" to point "<<PSP<<"."<<endl;
          IsValidSnapPoint=false;
        }
        
  /*      if(NewArea<OldArea*1./100.)
        {
          cout<<"Sorry, but you are not allowed to move point "<<DeadNode<<" to point "<<PSP<<"."<<endl;
          IsValidSnapPoint=false;
        }
        
        if(abs(cross)>10e-4)
        {
          cout<<"Sorry, but you are not allowed to move point "<<DeadNode<<" to point "<<PSP<<"."<<endl;
          IsValidSnapPoint=false;
        }*/
        
        //mutated cell
        MutatedCells.insert(C);
        cout<<"cell "<<C<<" has been infected!"<<endl;
      }
    }
    
    if(N_cells+N_newcells<=0)//survivor check
    {
      cout<<"Sorry, but you are not allowed to move point "<<DeadNode<<" to point "<<PSP<<"."<<endl;
      IsValidSnapPoint=false;
    }
/*    if(EmptyVolume(DeadNode,PSP))//simplified volume check
    {
      cout<<"Sorry, but you are not allowed to move point "<<DeadNode<<" to point "<<PSP<<"."<<endl;
      IsValidSnapPoint=false;
    }*/
    if(IsValidSnapPoint) {SnapPoint=PSP; break;}
  }//end of loop through potential SnapPoints
  
  cout<<"===>SNAPPOINT="<<SnapPoint<<endl;
  if(SnapPoint<0) {cout<<"Sorry no possible SnapPoint found."<<endl; return(false);}
  
  //allocate
  cout<<"N_points="<<N_points<<endl;
  cout<<"N_cells="<<N_cells<<endl;
  cout<<"N_newpoints="<<N_newpoints<<endl;
  cout<<"N_newcells="<<N_newcells<<endl;
  allocateGrid(dst,N_cells+N_newcells,N_points+N_newpoints);
  
  //vector used to redefine the new point IDs
  QVector <vtkIdType> OffSet(N_points);
  
  //copy undead points
  vtkIdType dst_id_node=0;
  for (vtkIdType src_id_node = 0; src_id_node < N_points; src_id_node++) {//loop through src points
    if(src_id_node!=DeadNode)//if the node isn't dead, copy it
    {
      vec3_t x;
      src->GetPoints()->GetPoint(src_id_node, x.data());
      dst->GetPoints()->SetPoint(dst_id_node, x.data());
      copyNodeData(src, src_id_node, dst, dst_id_node);
      OffSet[src_id_node]=src_id_node-dst_id_node;
      dst_id_node++;
    }
  };
  
  cout<<"DeadCells="<<DeadCells<<endl;
  cout<<"MutatedCells="<<MutatedCells<<endl;
  cout<<"MutilatedCells="<<MutilatedCells<<endl;
  
  //Copy undead cells
  for (vtkIdType id_cell = 0; id_cell < src->GetNumberOfCells(); ++id_cell) {//loop through src cells
    if(!DeadCells.contains(id_cell))//if the cell isn't dead
    {
      vtkIdType src_N_pts, *src_pts;
      vtkIdType dst_N_pts, *dst_pts;
      src->GetCellPoints(id_cell, src_N_pts, src_pts);
      
      vtkIdType type_cell = src->GetCellType(id_cell);
      cout<<"-->id_cell="<<id_cell<<endl;
      for(int i=0;i<src_N_pts;i++) cout<<"src_pts["<<i<<"]="<<src_pts[i]<<endl;
//       src->GetCellPoints(id_cell, dst_N_pts, dst_pts);
      dst_N_pts=src_N_pts;
      dst_pts=new vtkIdType[dst_N_pts];
      if(MutatedCells.contains(id_cell))//mutated cell
      {
        cout<<"processing mutated cell "<<id_cell<<endl;
        for(int i=0;i<src_N_pts;i++)
        {
          if(src_pts[i]==DeadNode) {
            cout<<"SnapPoint="<<SnapPoint<<endl;
            cout<<"OffSet[SnapPoint]="<<OffSet[SnapPoint]<<endl;
            cout<<"src_pts["<<i<<"]="<<src_pts[i]<<endl;
            dst_pts[i]=SnapPoint-OffSet[SnapPoint];
          }
          else dst_pts[i]=src_pts[i]-OffSet[src_pts[i]];
        }
        cout<<"--->dst_pts:"<<endl;
        for(int i=0;i<dst_N_pts;i++) cout<<"dst_pts["<<i<<"]="<<dst_pts[i]<<endl;
        
      }
      else if(MutilatedCells.contains(id_cell))//mutilated cell
      {
        cout<<"processing mutilated cell "<<id_cell<<endl;
        
        if(type_cell==VTK_QUAD) {
          type_cell=VTK_TRIANGLE;
          dst_N_pts-=1;
        }
        else {cout<<"FATAL ERROR: Unknown mutilated cell detected! It is not a quad! Potential xenomorph infestation!"<<endl;EG_BUG;}
        //merge points
        int j=0;
        for(int i=0;i<src_N_pts;i++)
        {
          if(src_pts[i]==SnapPoint) { dst_pts[j]=SnapPoint-OffSet[SnapPoint];j++; }//SnapPoint
          else if(src_pts[i]!=DeadNode) { dst_pts[j]=src_pts[i]-OffSet[src_pts[i]];j++; }//pre-snap/dead + post-snap/dead
          //do nothing in case of DeadNode
        }
      }
      else//normal cell
      {
        cout<<"processing normal cell "<<id_cell<<endl;
        for(int i=0;i<src_N_pts;i++)
        {
          dst_pts[i]=src_pts[i]-OffSet[src_pts[i]];
        }
      }
      //copy the cell
      vtkIdType id_new_cell = dst->InsertNextCell(type_cell, dst_N_pts, dst_pts);
      copyCellData(src, id_cell, dst, id_new_cell);
      cout<<"===Copying cell "<<id_cell<<" to "<<id_new_cell<<"==="<<endl;
      cout<<"src_pts:"<<endl;
      for(int i=0;i<src_N_pts;i++) cout<<"src_pts["<<i<<"]="<<src_pts[i]<<endl;
      cout<<"dst_pts:"<<endl;
      for(int i=0;i<dst_N_pts;i++) cout<<"dst_pts["<<i<<"]="<<dst_pts[i]<<endl;
      cout<<"OffSet="<<OffSet<<endl;
      cout<<"===Copying cell end==="<<endl;
      delete dst_pts;
    }
  };
//   cout_grid(cout,dst,true,true,true,true);
  makeCopy(dst, src);
  return(true);
}

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
