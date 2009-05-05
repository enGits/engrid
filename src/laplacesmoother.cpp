#include "laplacesmoother.h"
#include <vtkCellLocator.h>
#include <vtkCharArray.h>
#include <vtkGenericCell.h>

LaplaceSmoother::LaplaceSmoother()
{
   DebugLevel=0;
}

LaplaceSmoother::~LaplaceSmoother()
{
}

void LaplaceSmoother::operate()
{
  if(DebugLevel>10) cout<<"LaplaceSmoother reporting in."<<endl;
  
  QVector<vtkIdType> AllCells;
  getAllSurfaceCells(AllCells, m_grid);
  QVector<vtkIdType> SelectedCells;
  getSurfaceCells(m_bcs, SelectedCells, m_grid);
  
  cout<<"setCells START"<<endl;
  setCells(AllCells);
  cout<<"setCells END"<<endl;
  
  QSet <vtkIdType> SelectedNodes;
  getSurfaceNodes(m_bcs,SelectedNodes,m_grid);
  
  EG_VTKDCC(vtkIntArray, cell_code, m_grid, "cell_code");
  
  EG_VTKSP(vtkUnstructuredGrid,grid_tmp);
  EG_VTKSP(vtkUnstructuredGrid,m_grid_orig);
  makeCopy(m_grid, m_grid_orig);
  
  double closestPoint[3];
  vtkIdType cellId;
  int subId;
  double dist2;
  vtkCellLocator* terminator=vtkCellLocator::New();
  terminator->SetDataSet(m_grid_orig);
  terminator->BuildLocator();
//   terminator->CacheCellBoundsOn();
  cout<<"terminator->GetNumberOfBuckets()="<<terminator->GetNumberOfBuckets()<<endl;
  cout<<"terminator->GetNumberOfCellsPerBucket()="<<terminator->GetNumberOfCellsPerBucket()<<endl;
  cout<<"terminator->GetCacheCellBounds()="<<terminator->GetCacheCellBounds()<<endl;
  
  vtkGenericCell * cell=vtkGenericCell::New();
  
  UpdateNodeType_all();
  EG_VTKDCN(vtkCharArray, node_type, m_grid, "node_type");
  int moved_points=0;
  
  cout<<"makeCopy(m_grid, grid_tmp); START"<<endl;
  makeCopy(m_grid, grid_tmp);
  cout<<"makeCopy(m_grid, grid_tmp); END"<<endl;
  for(int i_iter=0;i_iter<NumberOfIterations;i_iter++)
  {
//     if(DebugLevel>10) 
      cout<<"i_iter="<<i_iter<<endl;
    
    foreach(vtkIdType id_G,SelectedNodes)
    {
      if(node_type->GetValue(id_G)==VTK_SIMPLE_VERTEX)
      {
        vec3_t G(0,0,0);
        foreach(int id_M,n2n[id_G])
        {
          vec3_t M;
          grid_tmp->GetPoint(id_M, M.data());
          G+=M;
        }
        G=(1./n2n[id_G].size())*G;
        vec3_t P;
//         cout<<"Searching for target "<<id_G<<"..."<<endl;
        terminator->FindClosestPoint(G.data(),P.data(),cellId,subId,dist2);
//         terminator->FindClosestPoint(G.data(),P.data(),cell,cellId,subId,dist2);
//         cout<<"Target destroyed."<<endl;
        grid_tmp->GetPoints()->SetPoint(id_G, G.data());
        moved_points++;
      }
    }
  }
  cout<<"makeCopy(grid_tmp,m_grid); START"<<endl;
  makeCopy(grid_tmp,m_grid);
  cout<<"makeCopy(grid_tmp,m_grid); END"<<endl;
  
  if(DebugLevel>10) cout << "SelectedNodes.size()=" << SelectedNodes.size() << endl;
  if(DebugLevel>10) cout << "moved_points=" << moved_points << endl;
  if(DebugLevel>10) cout_grid(cout,m_grid);
  
}

// void LaplaceSmoother::operate()
// {
//   if(DebugLevel>10) cout<<"LaplaceSmoother reporting in."<<endl;
//   
//   QVector<vtkIdType> AllCells;
//   getAllSurfaceCells(AllCells, m_grid);
//   QVector<vtkIdType> SelectedCells;
//   getSurfaceCells(m_bcs, SelectedCells, m_grid);
//   
//   EG_VTKDCC(vtkIntArray, cell_code, m_grid, "cell_code");
//   createCellToCell(AllCells, c2c, m_grid);
//   
//   QSet <vtkIdType> SelectedNodes;
//   QSet <vtkIdType> InternalNodes;
//   QSet <vtkIdType> ExternalNodes;
//   
//   foreach(vtkIdType id_cell, SelectedCells)
//   {
//     vtkIdType N_pts, *pts;
//     m_grid->GetCellPoints(id_cell, N_pts, pts);
//     for(int i=0;i<N_pts;i++)
//     {
//       QSet <int> bc;
//       foreach(vtkIdType C, n2c[pts[i]])
//       {
//         bc.insert(cell_code->GetValue(C));
//       }
//       if(DebugLevel>10) cout<<"pts[i]="<<pts[i]<<" and bc="<<bc<<endl;
//       SelectedNodes.insert(pts[i]);
//       if(bc.size()>1) ExternalNodes.insert(pts[i]);
//       else
//       {
//         vtkIdType point=pts[i];
//         QSet< int > NeighbourCells=n2c[point];
//         vtkIdType start=*(NeighbourCells.begin());
//         vtkIdType current=start;
//         do
//         {
//           vtkIdType next=nextcell(current,point,c2c,m_grid);
//           current=next;
//         } while (current!=start && current!=-1);
//         if(current==-1) ExternalNodes.insert(point);
//         if(current==start) InternalNodes.insert(point);
//       }
//     }
//   }
//   
//   createNodeToNode(cells, nodes, _nodes, n2n, m_grid);
//   
//   EG_VTKSP(vtkUnstructuredGrid,grid_tmp);
//   EG_VTKSP(vtkUnstructuredGrid,m_grid_orig);
//   makeCopy(m_grid, m_grid_orig);
//   
//   double closestPoint[3];
//   vtkIdType cellId;
//   int subId;
//   double dist2;
//   vtkCellLocator* terminator=vtkCellLocator::New();
//   terminator->SetDataSet(m_grid_orig);
//   terminator->BuildLocator();
//   
//   for(int i_iter=0;i_iter<NumberOfIterations;i_iter++)
//   {
//     if(DebugLevel>10) cout<<"i_iter="<<i_iter<<endl;
//     makeCopy(m_grid, grid_tmp);
//     
//     foreach(vtkIdType id_G,InternalNodes)
//     {
//       vec3_t G(0,0,0);
//       foreach(int id_M,n2n[id_G])
//       {
//         vec3_t M;
//         m_grid->GetPoint(id_M, M.data());
//         G+=M;
//       }
//       G=(1./n2n[id_G].size())*G;
//       vec3_t P;
//       terminator->FindClosestPoint(G.data(),P.data(),cellId,subId,dist2);
//       grid_tmp->GetPoints()->SetPoint(id_G, P.data());
//     }
//     
//     if(DebugLevel>10) cout << "SelectedNodes.size()=" << SelectedNodes.size() << endl;
//     if(DebugLevel>10) cout << "InternalNodes.size()=" << InternalNodes.size() << endl;
//     if(DebugLevel>10) cout << "ExternalNodes.size()=" << ExternalNodes.size() << endl;
//     if(DebugLevel>10) cout << "InternalNodes=" << InternalNodes << endl;
//     
//     makeCopy(grid_tmp,m_grid);
//   }
//   if(DebugLevel>10) cout_grid(cout,m_grid);
//   
// }
