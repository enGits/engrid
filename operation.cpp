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
