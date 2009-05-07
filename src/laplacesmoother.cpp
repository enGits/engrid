#include "laplacesmoother.h"
#include <vtkCellLocator.h>
#include <vtkCharArray.h>
#include <vtkGenericCell.h>
#include "guimainwindow.h"

using namespace GeometryTools;

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
          m_grid->GetPoint(id_M, M.data());
          G+=M;
        }
        G=(1./n2n[id_G].size())*G;
        vec3_t P;
//         cout<<"Searching for target "<<id_G<<"..."<<endl;
        terminator->FindClosestPoint(G.data(),P.data(),cellId,subId,dist2);
//         terminator->FindClosestPoint(G.data(),P.data(),cell,cellId,subId,dist2);
//         cout<<"Target destroyed."<<endl;
        
        //check that no cell gets flipped!
        while(FlippedCells(id_G,P))
        {
          vec3_t x0_old;
          m_grid->GetPoint(id_G, x0_old.data());
          P=x0_old+0.5*(P-x0_old);
        };
        
        m_grid->GetPoints()->SetPoint(id_G, P.data());
        
        int save=GuiMainWindow::pointer()->QuickSave();
        cout<<"save="<<save<<" : Moving "<<id_G<<" to "<<P<<endl;
        moved_points++;
      }
    }
  }
  
  if(DebugLevel>10) cout << "SelectedNodes.size()=" << SelectedNodes.size() << endl;
  if(DebugLevel>10) cout << "moved_points=" << moved_points << endl;
  if(DebugLevel>10) cout_grid(cout,m_grid);
  
}

bool LaplaceSmoother::FlippedCells(vtkIdType id_G, vec3_t P)
{
  vec3_t x0_old, x0_new;
  m_grid->GetPoint(id_G, x0_old.data());
  x0_new=P;
  
  cout_grid(cout,grid,true,true,true,true);
  foreach(vtkIdType id_cell,n2c[id_G])
  {
    vtkIdType N_pts, *pts;
    m_grid->GetCellPoints(id_cell, N_pts, pts);
    int i;
    for(i=0;i<N_pts;i++)
    {
      if(pts[i]==id_G) break;
    }
    vec3_t x2, x3;
    m_grid->GetPoint(pts[(i+1)%N_pts], x2.data());
    m_grid->GetPoint(pts[(i+2)%N_pts], x3.data());
    cout<<"Testing cell "<<id_cell<<": x0="<<id_G<<" x2="<<pts[(i+1)%N_pts]<<" x3="<<pts[(i+2)%N_pts]<<endl;
    vec3_t v2_old=x2-x0_old;
    vec3_t v3_old=x3-x0_old;
    
    //top point
    vec3_t S=v2_old.cross(v3_old);
    double V_old=tetraVol(x0_old, S, x2, x3, true);
    double V_new=tetraVol(x0_new, S, x2, x3, true);
    double prod=V_old*V_new;
    if( prod<0 ) {
      int save=GuiMainWindow::pointer()->QuickSave();
      cout<<"save="<<save<<" : Moving "<<id_G<<" to "<<P<<endl;
      cout<<"EPIC FAIL for id_G="<<id_G<<"!"<<endl;
      return(true);
    }
  }
  return(false);
}
