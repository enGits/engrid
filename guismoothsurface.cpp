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
#include "guismoothsurface.h"
#include "swaptriangles.h"
#include <vtkSmoothPolyDataFilter.h>
#include <vtksmoothpolydatafilter2.h>
#include <vtkWindowedSincPolyDataFilter.h>

#include <vtkLongArray.h>
#include <QtGui>

#include <QTextStream>
#include <stdio.h>

#include "vtkEgNormalExtrusion.h"
#include "containertricks.h"

#include <iostream>
#include <fstream>
#include <vtkIdList.h>
#include <vtkCell.h>
#include <cmath>
#include <vtkCellLocator.h>
#include <vtkFloatArray.h>
#include <vtkCellArray.h>

///////////////////////////////////////////
/* Here is how we we get QTextStreams that look like iostreams */
QTextStream Qcin(stdin, QIODevice::ReadOnly);
QTextStream Qcout(stdout, QIODevice::WriteOnly);
QTextStream Qcerr(stderr, QIODevice::WriteOnly);
///////////////////////////////////////////

//simple_vertices
//interior_edge_vertices
//fixed_vertices
int CreateSpecialMapping(QSet <vtkIdType> &simple_vertices, QSet <vtkIdType> &interior_edge_vertices, QSet <vtkIdType> &fixed_vertices)
{

  return(0);
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
//////////////////////////////////////////////
vtkIdType nextcell(vtkIdType a_cell, vtkIdType a_node, QVector< QVector< int > > &a_c2c, vtkUnstructuredGrid *a_grid)
{
  vtkIdType N_pts, *pts;
  a_grid->GetCellPoints(a_cell, N_pts, pts);
  
  int i;
  for(i=0;i<N_pts;i++)
  {
    if(pts[i]==a_node) break;
  }
  return a_c2c[a_cell][i];
}

pair<vtkIdType,vtkIdType> OrderedPair(vtkIdType a, vtkIdType b)
{
  vtkIdType x=min(a,b);
  vtkIdType y=max(a,b);
  return(pair<vtkIdType,vtkIdType>(x,y));
}

int cout_vtkWindowedSincPolyDataFilter(vtkWindowedSincPolyDataFilter* smooth)
{
  cout<<"NumberOfIterations="<<smooth->GetNumberOfIterations()<<endl;
  cout<<"PassBand="<<smooth->GetPassBand()<<endl;
  cout<<"FeatureEdgeSmoothing="<<smooth->GetFeatureEdgeSmoothing()<<endl;
  cout<<"FeatureAngle="<<smooth->GetFeatureAngle()<<endl;
  cout<<"EdgeAngle="<<smooth->GetEdgeAngle()<<endl;
  cout<<"BoundarySmoothing="<<smooth->GetBoundarySmoothing()<<endl;
  cout<<"GenerateErrorScalars="<<smooth->GetGenerateErrorScalars()<<endl;
  cout<<"GenerateErrorVectors="<<smooth->GetGenerateErrorVectors()<<endl;
  return(0);
}

int cout_vtkSmoothPolyDataFilter(vtkSmoothPolyDataFilter* smooth)
{
  cout<<"GetConvergence="<<smooth->GetConvergence()<<endl;
  cout<<"GetNumberOfIterations="<<smooth->GetNumberOfIterations()<<endl;
  cout<<"GetRelaxationFactor="<<smooth->GetRelaxationFactor()<<endl;
  cout<<"GetFeatureEdgeSmoothing="<<smooth->GetFeatureEdgeSmoothing()<<endl;
  cout<<"GetFeatureAngle="<<smooth->GetFeatureAngle()<<endl;
  cout<<"GetEdgeAngle="<<smooth->GetEdgeAngle()<<endl;
  cout<<"GetBoundarySmoothing="<<smooth->GetBoundarySmoothing()<<endl;
  cout<<"GetGenerateErrorScalars="<<smooth->GetGenerateErrorScalars()<<endl;
  cout<<"GetGenerateErrorVectors="<<smooth->GetGenerateErrorVectors()<<endl;
  return(0);
}
///////////////////////////////////////////

void GuiSmoothSurface::before()
{
  populateBoundaryCodes(ui.listWidget);
  populateBoundaryCodes(ui.listWidget_Source);
  
  ui.SmoothMethod->addItem("Method 0: vtkSmoothPolyDataFilter smoothing");
  ui.SmoothMethod->addItem("Method 1: vtkWindowedSincPolyDataFilter smoothing");
  ui.SmoothMethod->addItem("Method 2: edge subdivision");
  ui.SmoothMethod->addItem("Method 3: swap triangles");
  ui.SmoothMethod->addItem("Method 4: center subdivision");
  ui.SmoothMethod->addItem("Method 5: boundary refinement");
  ui.SmoothMethod->addItem("Method 6: Laplacian smoothing");
  ui.SmoothMethod->addItem("Method 7");
  ui.SmoothMethod->addItem("Method 8");
  ui.SmoothMethod->addItem("Method 9");
  ui.SmoothMethod->addItem("Method 10");
  ui.SmoothMethod->addItem("Method 11");
  ui.SmoothMethod->addItem("Method 12");
  ui.SmoothMethod->addItem("Method 13");
  ui.SmoothMethod->addItem("Method 14");
  ui.SmoothMethod->addItem("Method 15");
  ui.SmoothMethod->addItem("Method 16");
  ui.SmoothMethod->addItem("Method 17");
  ui.SmoothMethod->addItem("Method 18");
  
  ui.SmoothMethod->setCurrentIndex(10);
  
  if(ui.listWidget->count()>0) ui.lineEdit_BoundaryCode-> setText(ui.listWidget->item(0)->text());
  else ui.lineEdit_BoundaryCode-> setText("42");
  ui.spinBox_NumberOfSubdivisions->setValue(1);
  
  vtkSmoothPolyDataFilter* smooth=vtkSmoothPolyDataFilter::New();
  vtkWindowedSincPolyDataFilter* smooth2=vtkWindowedSincPolyDataFilter::New();
  
  cout_vtkSmoothPolyDataFilter(smooth);
  cout_vtkWindowedSincPolyDataFilter(smooth2);
  
  ui.doubleSpinBox_Convergence->setValue(smooth->GetConvergence());
//   ui.spinBox_NumberOfIterations->setValue(smooth->GetNumberOfIterations());
  ui.spinBox_NumberOfIterations->setValue(1000);
  QString tmp;
  ui.lineEdit_RelaxationFactor->setText(tmp.setNum(smooth->GetRelaxationFactor()));
  ui.doubleSpinBox_PassBand->setValue(smooth2->GetPassBand());
//   ui.checkBox_FeatureEdgeSmoothing->setCheckState(int2CheckState(smooth->GetFeatureEdgeSmoothing()));
  ui.checkBox_FeatureEdgeSmoothing->setCheckState(Qt::Checked);
  ui.doubleSpinBox_FeatureAngle->setValue(smooth->GetFeatureAngle());
  ui.doubleSpinBox_EdgeAngle->setValue(smooth->GetEdgeAngle());
  ui.checkBox_BoundarySmoothing->setCheckState(int2CheckState(smooth->GetBoundarySmoothing()));
/*  ui.checkBox_GenerateErrorScalars->setCheckState(int2CheckState(smooth->GetGenerateErrorScalars()));
  ui.checkBox_GenerateErrorVectors->setCheckState(int2CheckState(smooth->GetGenerateErrorVectors()));*/
  ui.checkBox_GenerateErrorScalars->setCheckState(int2CheckState(1));
  ui.checkBox_GenerateErrorVectors->setCheckState(int2CheckState(1));
  
};

void GuiSmoothSurface::operate()
{
  cout<<"METHOD "<<ui.SmoothMethod->currentIndex()<<endl;
  //can't use switch case because dynamic variables seem to be forbidden inside case statements
  //////////////////////////////////////////////////////////////////////////////////////////////
  if(ui.SmoothMethod->currentIndex()==0)//vtkSmoothPolyDataFilter smoothing
  {
    QSet<int> bcs;
    getSelectedItems(ui.listWidget, bcs);
    QVector<vtkIdType> cells;
    getSurfaceCells(bcs, cells, grid);
    EG_VTKSP(vtkPolyData, pdata);
    addToPolyData(cells, pdata, grid);
    EG_VTKSP(vtkSmoothPolyDataFilter, smooth);
    
    cout_vtkSmoothPolyDataFilter(smooth);
    
    smooth->SetInput(pdata);
    
    smooth->SetConvergence (ui.doubleSpinBox_Convergence->value());
    smooth->SetNumberOfIterations (ui.spinBox_NumberOfIterations->value());
    smooth->SetRelaxationFactor (ui.lineEdit_RelaxationFactor->text().toDouble());
    smooth->SetFeatureEdgeSmoothing (ui.checkBox_FeatureEdgeSmoothing->checkState());
    smooth->SetFeatureAngle (ui.doubleSpinBox_FeatureAngle->value());
    smooth->SetEdgeAngle (ui.doubleSpinBox_EdgeAngle->value());
    smooth->SetBoundarySmoothing (ui.checkBox_BoundarySmoothing->checkState());
    smooth->SetGenerateErrorScalars (ui.checkBox_GenerateErrorScalars->checkState());
    smooth->SetGenerateErrorVectors (ui.checkBox_GenerateErrorVectors->checkState());
    
    QSet<int> bcs_Source;
    getSelectedItems(ui.listWidget,bcs_Source);
    QVector<vtkIdType> cells_Source;
    getSurfaceCells(bcs_Source, cells_Source, grid);
    EG_VTKSP(vtkPolyData, pdata_Source);
    addToPolyData(cells_Source, pdata_Source, grid);
    smooth->SetSource (pdata_Source);
    
    cout_vtkSmoothPolyDataFilter(smooth);
    
    smooth->Update();
    
//     grid->GetPointCells();
    cout<<"==============="<<endl;
    int N1,N2;
    double x1[3], x2[3], x3[3], l1[3], l2[3];
    double dist;
    int numPts=0;
//     int indexArray[10];
    cout<<"ErrorScalars:"<<endl;
/*    int idx = output->GetPointData()->AddArray(newScalars);
    grid->GetPointData()->SetActiveAttribute(idx, vtkDataSetAttributes::SCALARS);*/
//     grid->GetPointData()->GetAttribute();
    N1=smooth->GetOutput()->GetPointData()->GetNumberOfArrays();
    cout<<"nb of arrays="<<N1<<endl;
//     indexArray=new int[N1];
    smooth->GetOutput();//vtkPolyData
    cout<<smooth->GetOutput()->GetPointData()<<endl;//vtkPointData*
//     cout<<"getting indices"<<endl;
//     smooth->GetOutput()->GetPointData()->GetAttributeIndices(indexArray);
//     for(int i=0;i<N1;i++) cout<<"indexArray["<<i<<"]="<<indexArray[i]<<endl;
//     delete[] indexArray;
    vtkFloatArray *newScalars = vtkFloatArray::New();
    newScalars=(vtkFloatArray *)smooth->GetOutput()->GetPointData()->GetArray(1);
//     newScalars=smooth->GetOutput()->GetPointData()->GetAttribute(vtkDataSetAttributes::SCALARS);
    N1=newScalars->GetNumberOfComponents();
    N2=newScalars->GetNumberOfTuples();
    cout<<"N1="<<N1<<endl;
    cout<<"N2="<<N2<<endl;
    for (int i=0; i<N2; i++)
    {
      dist=newScalars->GetComponent(i-1,1);//strange, but works. O.o
      cout<<"dist="<<dist<<endl;
    }
    
    cout<<"ErrorVectors:"<<endl;
    int index;
    vtkPointData* toto=vtkPointData::New();
    toto=grid->GetPointData();
    vtkDataArray* titi=(vtkDataArray* ) vtkDataArray::New();
    titi=toto->GetVectors();
//     titi->GetTuple(0,x3);
    cout<<smooth->GetOutput()->GetPointData()<<endl;
    cout<<smooth->GetOutput()->GetPoints()<<endl;
    cout<<smooth->GetOutput()->GetPointData()->GetVectors()<<endl;
//     cout<<smooth->GetOutput()->GetVectors()<<endl;
    vec3_t xx;
    smooth->GetOutput()->GetPoint(0, xx.data());
    N1=smooth->GetOutput()->GetPointData()->GetVectors()->GetNumberOfComponents();
    N2=smooth->GetOutput()->GetPointData()->GetVectors()->GetNumberOfTuples();
    cout<<"N1="<<N1<<endl;
    cout<<"N2="<<N2<<endl;
    for(vtkIdType i=0;i<N2;i++)
    {
      double tuple[4];
      smooth->GetOutput()->GetPointData()->GetTuple(i,tuple);
      cout<<"tuple["<<tuple[0]<<"]=("<<tuple[1]<<","<<tuple[2]<<","<<tuple[3]<<")"<<endl;
    }
//     ->GetArray(,index);
    cout<<"index="<<index<<endl;
    cout<<"==============="<<endl;
    
    EG_VTKDCN(vtkLongArray_t, node_index, pdata, "node_index");
    for (vtkIdType i = 0; i < smooth->GetOutput()->GetNumberOfPoints(); ++i) {
      vec3_t x;
      smooth->GetOutput()->GetPoints()->GetPoint(i, x.data());
      vtkIdType nodeId = node_index->GetValue(i);
      grid->GetPoints()->SetPoint(nodeId, x.data());
    };
    updateActors();
  }
  //////////////////////////////////////////////////////////////////////////////////////////////
  else if(ui.SmoothMethod->currentIndex()==1)//vtkWindowedSincPolyDataFilter smoothing
  {
    QSet<int> bcs;
    getSelectedItems(ui.listWidget, bcs);
    QVector<vtkIdType> cells;
    getSurfaceCells(bcs, cells, grid);
    EG_VTKSP(vtkPolyData, pdata);
    addToPolyData(cells, pdata, grid);
    EG_VTKSP(vtkWindowedSincPolyDataFilter, smooth);
    
    cout_vtkWindowedSincPolyDataFilter(smooth);
    
    smooth->SetInput(pdata);
    
    smooth->SetNumberOfIterations(ui.spinBox_NumberOfIterations->value());
    smooth->SetPassBand(ui.doubleSpinBox_PassBand->value());
    smooth->SetFeatureEdgeSmoothing(ui.checkBox_FeatureEdgeSmoothing->checkState());
    smooth->SetFeatureAngle(ui.doubleSpinBox_FeatureAngle->value());
    smooth->SetEdgeAngle(ui.doubleSpinBox_EdgeAngle->value());
    smooth->SetBoundarySmoothing(ui.checkBox_BoundarySmoothing->checkState());
    smooth->SetGenerateErrorScalars(ui.checkBox_GenerateErrorScalars->checkState());
    smooth->SetGenerateErrorVectors(ui.checkBox_GenerateErrorVectors->checkState());
    
    cout_vtkWindowedSincPolyDataFilter(smooth);
    
    smooth->Update();
    EG_VTKDCN(vtkLongArray_t, node_index, pdata, "node_index");
    for (vtkIdType i = 0; i < smooth->GetOutput()->GetNumberOfPoints(); ++i) {
      vec3_t x;
      smooth->GetOutput()->GetPoints()->GetPoint(i, x.data());
      vtkIdType nodeId = node_index->GetValue(i);
      grid->GetPoints()->SetPoint(nodeId, x.data());
    };
    updateActors();
  }
  //////////////////////////////////////////////////////////////////////////////////////////////
  else if(ui.SmoothMethod->currentIndex()==2)//edge subdivision
  {
    cout_grid(cout,grid);
    
    QSet<int> bcs;
    getSelectedItems(ui.listWidget, bcs);
    
    int N_iter=ui.spinBox_NumberOfSubdivisions->value();
    for(int i_iter=0;i_iter<N_iter;i_iter++)
    {
      cout<<"i_iter="<<i_iter<<endl;
      
      QVector<vtkIdType> SelectedCells;
      getSurfaceCells(bcs, SelectedCells, grid);
      QVector<vtkIdType> AllCells;
      getAllSurfaceCells(AllCells,grid);
      
      createCellToCell(AllCells, c2c, grid);
      
      int N_points=grid->GetNumberOfPoints();
      int N_cells=grid->GetNumberOfCells();
      
      QMap< pair<vtkIdType,vtkIdType>, vtkIdType> midpoint_map;
//       QMap<double, int> midpoint_map;
      int N_extmidpoints=0;
      int N_intmidpoints=0;
      
      int N_newpoints=0;
      int N_newcells=0;
      vtkIdType nodeId = N_points;
      foreach(vtkIdType id_cell, SelectedCells)
      {
        vtkIdType type_cell = grid->GetCellType(id_cell);
        int N_neighbours=c2c[id_cell].size();
        for(int i=0;i<N_neighbours;i++)
        {
          vtkIdType id_neighbour=c2c[id_cell][i];
          if(id_neighbour<0)
          {
            N_extmidpoints++;
          }
          else
          {
            midpoint_map[OrderedPair(id_cell,id_neighbour)]=nodeId; nodeId++;
          }
        }
        if (type_cell == VTK_TRIANGLE)
        {
          N_newcells+=3;
        }
        if (type_cell == VTK_QUAD)
        {
          N_newcells+=3;
          N_intmidpoints++;
        }
      }
      
      int N_c2cmidpoints=midpoint_map.size();
      N_newpoints=N_c2cmidpoints+N_extmidpoints+N_intmidpoints;
      
      cout<<"N_c2cmidpoints="<<N_c2cmidpoints<<endl;
      cout<<"N_extmidpoints="<<N_extmidpoints<<endl;
      cout<<"N_intmidpoints="<<N_intmidpoints<<endl;
      cout<<"N_newpoints="<<N_newpoints<<endl;
      cout<<"N_newcells="<<N_newcells<<endl;
      cout<<"N_cells="<<N_cells<<endl;
      cout<<"N_points="<<N_points<<endl;
      
      EG_VTKSP(vtkUnstructuredGrid,grid_tmp);
      allocateGrid(grid_tmp,N_cells+N_newcells,N_points+N_newpoints);
      makeCopyNoAlloc(grid, grid_tmp);
      
      EG_VTKDCC(vtkIntArray, cell_code, grid_tmp, "cell_code");
      
      midpoint_map.clear();//clear midpoint_map
      nodeId=N_points;//reset nodeId
      
      foreach(vtkIdType id_cell, SelectedCells)
      {
        vtkIdType N_pts, *pts;
        grid->GetCellPoints(id_cell, N_pts, pts);
        vtkIdType intmidpoint;
        vtkIdType edgemidpoints[4];
        vec3_t M[4];
        
        vtkIdType type_cell = grid->GetCellType(id_cell);
        int N_neighbours=c2c[id_cell].size();
        vec3_t corner[4];
        for(int i=0;i<N_neighbours;i++)
        {
          grid->GetPoints()->GetPoint(pts[i], corner[i].data());
        }
        
        for(int i=0;i<N_neighbours;i++)
        {
          vtkIdType id_neighbour=c2c[id_cell][i];
          if(id_neighbour<0)
          {
            M[i]=0.5*(corner[i]+corner[(i+1)%N_neighbours]);
            addPoint(grid_tmp,nodeId,M[i].data());
            edgemidpoints[i]=nodeId;
            nodeId++;
          }
          else
          {
            if(midpoint_map.contains(OrderedPair(id_cell,id_neighbour)))
            {
              //pt already exists!
/*              cout<<"pt already exists!: i="<<i<<" midpoint_map[OrderedPair(id_cell,id_neighbour)]="<<midpoint_map[OrderedPair(id_cell,id_neighbour)]<<endl;*/
              edgemidpoints[i]=midpoint_map[OrderedPair(id_cell,id_neighbour)];
            }
            else
            {
              M[i]=0.5*(corner[i]+corner[(i+1)%N_neighbours]);
              addPoint(grid_tmp,nodeId,M[i].data());
              midpoint_map[OrderedPair(id_cell,id_neighbour)]=nodeId;
              edgemidpoints[i]=nodeId;
              nodeId++;
            }
          }
        }
        if (type_cell == VTK_TRIANGLE)
        {
          vtkIdType pts_triangle[4][3];
          pts_triangle[0][0]=pts[0];//A;
          pts_triangle[0][1]=edgemidpoints[0];
          pts_triangle[0][2]=edgemidpoints[2];
          pts_triangle[1][0]=pts[1];//B;
          pts_triangle[1][1]=edgemidpoints[1];
          pts_triangle[1][2]=edgemidpoints[0];
          pts_triangle[2][0]=pts[2];//C;
          pts_triangle[2][1]=edgemidpoints[2];
          pts_triangle[2][2]=edgemidpoints[1];
          pts_triangle[3][0]=edgemidpoints[0];
          pts_triangle[3][1]=edgemidpoints[1];
          pts_triangle[3][2]=edgemidpoints[2];
          
          vtkIdType newCellId;
          newCellId = grid_tmp->InsertNextCell(VTK_TRIANGLE,3,pts_triangle[0]);
          cell_code->SetValue(newCellId, ui.lineEdit_BoundaryCode->text().toInt());
          newCellId = grid_tmp->InsertNextCell(VTK_TRIANGLE,3,pts_triangle[1]);
          cell_code->SetValue(newCellId, ui.lineEdit_BoundaryCode->text().toInt());
          newCellId = grid_tmp->InsertNextCell(VTK_TRIANGLE,3,pts_triangle[2]);
          cell_code->SetValue(newCellId, ui.lineEdit_BoundaryCode->text().toInt());
          grid_tmp->ReplaceCell(id_cell , 3, pts_triangle[3]);
          cell_code->SetValue(id_cell, ui.lineEdit_BoundaryCode->text().toInt());
        }
        if (type_cell == VTK_QUAD)
        {
          vec3_t C=0.25*(corner[0]+corner[1]+corner[2]+corner[3]);
          addPoint(grid_tmp,nodeId,C.data());
          intmidpoint=nodeId;
          nodeId++;
          
          vtkIdType pts_quad[4][4];
          pts_quad[0][0]=pts[0];
          pts_quad[0][1]=edgemidpoints[0];
          pts_quad[0][2]=intmidpoint;
          pts_quad[0][3]=edgemidpoints[3];
          pts_quad[1][0]=pts[1];
          pts_quad[1][1]=edgemidpoints[1];
          pts_quad[1][2]=intmidpoint;
          pts_quad[1][3]=edgemidpoints[0];
          pts_quad[2][0]=pts[2];
          pts_quad[2][1]=edgemidpoints[2];
          pts_quad[2][2]=intmidpoint;
          pts_quad[2][3]=edgemidpoints[1];
          pts_quad[3][0]=pts[3];
          pts_quad[3][1]=edgemidpoints[3];
          pts_quad[3][2]=intmidpoint;
          pts_quad[3][3]=edgemidpoints[2];
          
          vtkIdType newCellId;
          newCellId = grid_tmp->InsertNextCell(VTK_QUAD,4,pts_quad[0]);
          cell_code->SetValue(newCellId, ui.lineEdit_BoundaryCode->text().toInt());
          newCellId = grid_tmp->InsertNextCell(VTK_QUAD,4,pts_quad[1]);
          cell_code->SetValue(newCellId, ui.lineEdit_BoundaryCode->text().toInt());
          newCellId = grid_tmp->InsertNextCell(VTK_QUAD,4,pts_quad[2]);
          cell_code->SetValue(newCellId, ui.lineEdit_BoundaryCode->text().toInt());
          grid_tmp->ReplaceCell(id_cell , 4, pts_quad[3]);
          cell_code->SetValue(id_cell, ui.lineEdit_BoundaryCode->text().toInt());
        }
      }
      
//       cout_grid(cout,grid_tmp,true,true,true,true);
      cout<<"Copying..."<<endl;
      makeCopy(grid_tmp,grid);
      cout<<"Copy successful"<<endl;
    }
    cout_grid(cout,grid);
    updateActors();
  }
  //////////////////////////////////////////////////////////////////////////////////////////////
  else if(ui.SmoothMethod->currentIndex()==3)//swap triangles
  {
    cout_grid(cout,grid);
    
    QSet<int> bcs;
    getSelectedItems(ui.listWidget, bcs);
    
    QSet<int> bcs_complement=complementary_bcs(bcs,grid,cells);
    
    cout<<"bcs="<<bcs<<endl;
    cout<<"bcs_complement="<<bcs_complement<<endl;
    
    SwapTriangles swap;
    swap.setGrid(grid);
    swap.setBoundaryCodes(bcs_complement);
    swap();
    
    cout_grid(cout,grid);
    updateActors();
  }
  //////////////////////////////////////////////////////////////////////////////////////////////
  else if(ui.SmoothMethod->currentIndex()==4)//center subdivision
  {
    cout_grid(cout,grid);
    
    QSet<int> bcs;
    getSelectedItems(ui.listWidget, bcs);
    QSet<int> bcs_complement=complementary_bcs(bcs,grid,cells);
    cout<<"bcs="<<bcs<<endl;
    cout<<"bcs_complement="<<bcs_complement<<endl;
    
    int N_iter=ui.spinBox_NumberOfSubdivisions->value();
    for(int i_iter=0;i_iter<N_iter;i_iter++)
    {
      int N_points=grid->GetNumberOfPoints();
      int N_cells=grid->GetNumberOfCells();
      
      QVector<vtkIdType> cells;
      getSurfaceCells(bcs, cells, grid);
      
      int N_newcells=0;
      foreach(int id_cell, cells)
      {
        vtkIdType type_cell = grid->GetCellType(id_cell);
        if (type_cell == VTK_TRIANGLE) N_newcells+=2;
        if (type_cell == VTK_QUAD) N_newcells+=3;
      }
      
      cout<<"N_newcells="<<N_newcells<<endl;
      
      EG_VTKSP(vtkUnstructuredGrid,grid_tmp);
      allocateGrid(grid_tmp,N_cells+N_newcells,N_points+1*N_cells);
      makeCopyNoAlloc(grid, grid_tmp);
      EG_VTKDCC(vtkIntArray, cell_code, grid_tmp, "cell_code");
      
      vtkIdType nodeId=N_points;
      foreach(int id_cell, cells)
      {
        vtkIdType N_pts, *pts;
        grid->GetCellPoints(id_cell, N_pts, pts);
        vec3_t C(0,0,0);
        
        vtkIdType type_cell = grid->GetCellType(id_cell);
        int N_neighbours=N_pts;
        cout<<"N_neighbours="<<N_neighbours<<endl;
        vec3_t corner[4];
        vtkIdType pts_triangle[4][3];
        for(int i=0;i<N_neighbours;i++)
        {
          grid->GetPoints()->GetPoint(pts[i], corner[i].data());
          C+=(1/(double)N_neighbours)*corner[i];
        }
        addPoint(grid_tmp,nodeId,C.data());
        vtkIdType intmidpoint=nodeId;
        nodeId++;
        
        for(int i=0;i<N_neighbours;i++)
        {
          pts_triangle[i][0]=pts[i];
          pts_triangle[i][1]=pts[(i+1)%N_neighbours];
          pts_triangle[i][2]=intmidpoint;
          if(i==0)
          {
            grid_tmp->ReplaceCell(id_cell , 3, pts_triangle[0]);
            cell_code->SetValue(id_cell, ui.lineEdit_BoundaryCode->text().toInt());
          }
          else
          {
            vtkIdType newCellId = grid_tmp->InsertNextCell(VTK_TRIANGLE,3,pts_triangle[i]);
            cell_code->SetValue(newCellId, ui.lineEdit_BoundaryCode->text().toInt());
          }
        }
      }
      
      makeCopy(grid_tmp,grid);
    }
    cout_grid(cout,grid,true,true,true,true);
    updateActors();
  }
  //////////////////////////////////////////////////////////////////////////////////////////////
  else if(ui.SmoothMethod->currentIndex()==5)//boundary refinement
  {
    cout_grid(cout,grid);
    
    QSet<int> bcs;
    getSelectedItems(ui.listWidget, bcs);
    QSet<int> bcs_complement=complementary_bcs(bcs,grid,cells);
    cout<<"bcs="<<bcs<<endl;
    cout<<"bcs_complement="<<bcs_complement<<endl;
    
    int N_points=grid->GetNumberOfPoints();
    int N_cells=grid->GetNumberOfCells();
    
    QVector<vtkIdType> SelectedCells;
    getSurfaceCells(bcs, SelectedCells, grid);
    QVector<vtkIdType> AllCells;
    getAllSurfaceCells(AllCells,grid);
    
    EG_VTKDCC(vtkIntArray, cell_code, grid, "cell_code");
    
    QSet <int> cells_to_split;
    QVector <stencil_t> StencilVector;
    
    QMap <int,bool> marked;
//     marked.resize(SelectedCells.size());
    foreach(vtkIdType id_cell, SelectedCells)
    {
      marked[id_cell] = false;
    }
    createCellToCell(AllCells, c2c, grid);
/*    for (int i_cells = 0; i_cells < SelectedCells.size(); ++i_cells) {
      marked[i_cells] = false;
    };*/
    
    cout<<"AllCells.size()="<<AllCells.size()<<endl;
    cout<<"SelectedCells.size()="<<SelectedCells.size()<<endl;
    
    int N_newcells=0;
    int N_newpoints=0;
    foreach(vtkIdType id_cell, SelectedCells)
    {
      cout<<"==>id_cell="<<id_cell<<endl;
      int bc0=cell_code->GetValue(id_cell);
      if(!marked[id_cell])
      {
        vtkIdType N_pts, *pts;
        grid->GetCellPoints(id_cell, N_pts, pts);
        int count=0;
        for(int i=0;i<N_pts;i++)
        {
          int bc1=cell_code->GetValue(c2c[id_cell][i]);
          if(bc0!=bc1) count++;
        }
        if(count>0)
        {
          int SideToSplit = getLongestSide(id_cell,grid);
          cout<<"SideToSplit="<<SideToSplit<<endl;
          cout<<"c2c[id_cell][SideToSplit]="<<c2c[id_cell][SideToSplit]<<endl;
          for(int i=0;i<3;i++) cout<<"c2c[id_cell]["<<i<<"]="<<c2c[id_cell][i]<<endl;
          stencil_t S=getStencil(id_cell,SideToSplit);
          if(!marked[S.id_cell2])
          {
            cells_to_split.insert(id_cell);
            cout<<"marked["<<S.id_cell1<<"]=true;"<<endl;
            cout<<"marked["<<S.id_cell2<<"]=true;"<<endl;
            marked[S.id_cell1]=true;
            marked[S.id_cell2]=true;
            StencilVector.push_back(S);
            N_newpoints++;
            if(S.valid) N_newcells+=2;
            else N_newcells+=1;
          }
        }
      }
    }
    
    cout<<"cells_to_split.size()="<<cells_to_split.size()<<endl;
    cout<<cells_to_split<<endl;
    cout<<"StencilVector.size()="<<StencilVector.size()<<endl;
    cout<<"N_newpoints="<<N_newpoints<<endl;
    cout<<"N_newcells="<<N_newcells<<endl;
    
    EG_VTKSP(vtkUnstructuredGrid,grid_tmp);
    allocateGrid(grid_tmp,N_cells+N_newcells,N_points+N_newpoints);
    makeCopyNoAlloc(grid, grid_tmp);
    EG_VTKDCC(vtkIntArray, cell_code_tmp, grid_tmp, "cell_code");
    
    vtkIdType nodeId=N_points;
    foreach(stencil_t S,StencilVector)
    {
      cout<<S<<endl;
      vtkIdType N_pts, *pts;
      vec3_t A,B;
      grid_tmp->GetPoint(S.p[1],A.data());
      grid_tmp->GetPoint(S.p[3],B.data());
      vec3_t M=0.5*(A+B);
      addPoint(grid_tmp,nodeId,M.data());
      
      vtkIdType pts_triangle[4][3];
      
      for(int i=0;i<4;i++)
      {
        pts_triangle[i][0]=S.p[i];
        pts_triangle[i][1]=S.p[(i+1)%4];
        pts_triangle[i][2]=nodeId;
      }
      
      int bc1=cell_code_tmp->GetValue(S.id_cell1);
      int bc2=cell_code_tmp->GetValue(S.id_cell2);
      
      grid_tmp->ReplaceCell(S.id_cell1 , 3, pts_triangle[0]);
      cell_code_tmp->SetValue(S.id_cell1, bc1);
      
      grid_tmp->ReplaceCell(S.id_cell2 , 3, pts_triangle[1]);
      cell_code_tmp->SetValue(S.id_cell2, bc2);
      
      vtkIdType newCellId;
      newCellId = grid_tmp->InsertNextCell(VTK_TRIANGLE,3,pts_triangle[2]);
      cell_code_tmp->SetValue(newCellId, bc2);
      newCellId = grid_tmp->InsertNextCell(VTK_TRIANGLE,3,pts_triangle[3]);
      cell_code_tmp->SetValue(newCellId, bc1);
      
      nodeId++;
    }
    
    makeCopy(grid_tmp,grid);
//     cout_grid(cout,grid,true,true,true,true);
    updateActors();
  }
  //////////////////////////////////////////////////////////////////////////////////////////////
  else if(ui.SmoothMethod->currentIndex()==6)
  {
    QSet<int> bcs;
    getSelectedItems(ui.listWidget, bcs);
    
    QVector<vtkIdType> AllCells;
    getAllSurfaceCells(AllCells, grid);
    QVector<vtkIdType> SelectedCells;
    getSurfaceCells(bcs, SelectedCells, grid);
    
    EG_VTKDCC(vtkIntArray, cell_code, grid, "cell_code");
    createCellToCell(AllCells, c2c, grid);
    
    QSet <vtkIdType> SelectedNodes;
    QSet <vtkIdType> InternalNodes;
    QSet <vtkIdType> ExternalNodes;
    
    foreach(vtkIdType id_cell, SelectedCells)
    {
      vtkIdType N_pts, *pts;
      grid->GetCellPoints(id_cell, N_pts, pts);
      for(int i=0;i<N_pts;i++)
      {
        QSet <int> bc;
        foreach(int C, n2c[pts[i]])
        {
          bc.insert(cell_code->GetValue(C));
        }
        cout<<"pts[i]="<<pts[i]<<" and bc="<<bc<<endl;
        SelectedNodes.insert(pts[i]);
        if(bc.size()>1) ExternalNodes.insert(pts[i]);
        else
        {
          vtkIdType point=pts[i];
          QSet< int > NeighbourCells=n2c[point];
          vtkIdType start=*(NeighbourCells.begin());
          vtkIdType current=start;
          do
          {
            vtkIdType next=nextcell(current,point,c2c,grid);
            current=next;
          } while (current!=start && current!=-1);
          if(current==-1) ExternalNodes.insert(point);
          if(current==start) InternalNodes.insert(point);
        }
      }
    }
    
    createNodeToNode(cells, nodes, _nodes, n2n, grid);
    
    EG_VTKSP(vtkUnstructuredGrid,grid_tmp);
    EG_VTKSP(vtkUnstructuredGrid,grid_orig);
    makeCopy(grid, grid_orig);
    
    double closestPoint[3];
    vtkIdType cellId;
    int subId;
    double dist2;
    vtkCellLocator* terminator=vtkCellLocator::New();
    terminator->SetDataSet(grid_orig);
    terminator->BuildLocator();
    
    int N_iter=ui.spinBox_NumberOfIterations->value();
    for(int i_iter=0;i_iter<N_iter;i_iter++)
    {
      cout<<"i_iter="<<i_iter<<endl;
      makeCopy(grid, grid_tmp);
      
      foreach(vtkIdType id_G,InternalNodes)
      {
        vec3_t G(0,0,0);
        foreach(int id_M,n2n[id_G])
        {
          vec3_t M;
          grid->GetPoint(id_M, M.data());
          G+=M;
        }
        G=(1./n2n[id_G].size())*G;
        vec3_t P;
        terminator->FindClosestPoint(G.data(),P.data(),cellId,subId,dist2);
        grid_tmp->GetPoints()->SetPoint(id_G, P.data());
      }
      
      cout << "SelectedNodes.size()=" << SelectedNodes.size() << endl;
      cout << "InternalNodes.size()=" << InternalNodes.size() << endl;
      cout << "ExternalNodes.size()=" << ExternalNodes.size() << endl;
      cout << "InternalNodes=" << InternalNodes << endl;
      
/*      QSet <vtkIdType> SelectedNodes2;
      getSurfaceNodes(bcs,SelectedNodes2,grid);
      cout << "SelectedNodes=" << SelectedNodes << endl;
      cout << "SelectedNodes2=" << SelectedNodes2 << endl;*/
      
      makeCopy(grid_tmp,grid);
    }
    cout_grid(cout,grid);
    updateActors();
  }
  //////////////////////////////////////////////////////////////////////////////////////////////
  else if(ui.SmoothMethod->currentIndex()==7)//VertexAvgDist test
  {
    cout_grid(cout,grid);
    
    QSet<int> bcs;
    getSelectedItems(ui.listWidget, bcs);
    QSet<int> bcs_complement=complementary_bcs(bcs,grid,cells);
    cout<<"bcs="<<bcs<<endl;
    cout<<"bcs_complement="<<bcs_complement<<endl;
    
    int N_points=grid->GetNumberOfPoints();
    int N_cells=grid->GetNumberOfCells();
    
    QVector<vtkIdType> SelectedCells;
    getSurfaceCells(bcs, SelectedCells, grid);
    QVector<vtkIdType> AllCells;
    getAllSurfaceCells(AllCells,grid);
    
    QSet <vtkIdType> SelectedNodes;
    getSurfaceNodes(bcs,SelectedNodes,grid);
    createNodeToNode(cells, nodes, _nodes, n2n, grid);
    
    foreach(vtkIdType node,SelectedNodes)
    {
      cout<<"node="<<node<<" VertexAvgDist="<<CurrentVertexAvgDist(node,n2n,grid)<<endl;
    }
    
    updateActors();
  }
  //////////////////////////////////////////////////////////////////////////////////////////////
  else if(ui.SmoothMethod->currentIndex()==8)//boundary refinement v2
  {
    cout_grid(cout,grid);
    
    QSet<int> bcs;
    getSelectedItems(ui.listWidget, bcs);
    QSet<int> bcs_complement=complementary_bcs(bcs,grid,cells);
    cout<<"bcs="<<bcs<<endl;
    cout<<"bcs_complement="<<bcs_complement<<endl;
    
    int N_points=grid->GetNumberOfPoints();
    int N_cells=grid->GetNumberOfCells();
    
    QVector<vtkIdType> SelectedCells;
    getSurfaceCells(bcs, SelectedCells, grid);
    QVector<vtkIdType> AllCells;
    getAllSurfaceCells(AllCells,grid);
    
    EG_VTKDCC(vtkIntArray, cell_code, grid, "cell_code");
    EG_VTKDCN(vtkDoubleArray, node_meshdensity, grid, "node_meshdensity");
    
    QSet <vtkIdType> SelectedNodes;
    getSurfaceNodes(bcs,SelectedNodes,grid);
    createNodeToNode(cells, nodes, _nodes, n2n, grid);
    
    foreach(vtkIdType node,SelectedNodes)
    {
      double L=CurrentVertexAvgDist(node,n2n,grid);
      double D=1./L;
      cout<<"node="<<node<<" VertexAvgDist="<<L<<" Net density="<<D<<endl;
      node_meshdensity->SetValue(node, D);
    }
    
    int N_iter=ui.spinBox_NumberOfIterations->value();
    for(int i_iter=0;i_iter<N_iter;i_iter++)
    {
      foreach(vtkIdType node,SelectedNodes)
      {
        double D=DesiredMeshDensity(node,n2n,grid);
        double L=1./D;
        cout<<"node="<<node<<" VertexAvgDist="<<L<<" Net density="<<D<<endl;
        node_meshdensity->SetValue(node, D);
      }
    }
    
/*    QSet <int> cells_to_split;
    QVector <stencil_t> StencilVector;
    
    QMap <int,bool> marked;
//     marked.resize(SelectedCells.size());
    foreach(vtkIdType id_cell, SelectedCells)
    {
      marked[id_cell] = false;
    }
    createCellToCell(AllCells, c2c, grid);
    
    cout<<"AllCells.size()="<<AllCells.size()<<endl;
    cout<<"SelectedCells.size()="<<SelectedCells.size()<<endl;
    
    int N_newcells=0;
    int N_newpoints=0;
    foreach(vtkIdType id_cell, SelectedCells)
    {
      cout<<"==>id_cell="<<id_cell<<endl;
      int bc0=cell_code->GetValue(id_cell);
      if(!marked[id_cell])
      {
        vtkIdType N_pts, *pts;
        grid->GetCellPoints(id_cell, N_pts, pts);
        int count=0;
        for(int i=0;i<N_pts;i++)
        {
          int bc1=cell_code->GetValue(c2c[id_cell][i]);
          if(bc0!=bc1) count++;
        }
        if(count>0)
        {
          int SideToSplit = getLongestSide(id_cell,grid);
          cout<<"SideToSplit="<<SideToSplit<<endl;
          cout<<"c2c[id_cell][SideToSplit]="<<c2c[id_cell][SideToSplit]<<endl;
          for(int i=0;i<3;i++) cout<<"c2c[id_cell]["<<i<<"]="<<c2c[id_cell][i]<<endl;
          stencil_t S=getStencil(id_cell,SideToSplit);
          if(!marked[S.id_cell2])
          {
            cells_to_split.insert(id_cell);
            cout<<"marked["<<S.id_cell1<<"]=true;"<<endl;
            cout<<"marked["<<S.id_cell2<<"]=true;"<<endl;
            marked[S.id_cell1]=true;
            marked[S.id_cell2]=true;
            StencilVector.push_back(S);
            N_newpoints++;
            if(S.valid) N_newcells+=2;
            else N_newcells+=1;
          }
        }
      }
    }
    
    cout<<"cells_to_split.size()="<<cells_to_split.size()<<endl;
    cout<<cells_to_split<<endl;
    cout<<"StencilVector.size()="<<StencilVector.size()<<endl;
    cout<<"N_newpoints="<<N_newpoints<<endl;
    cout<<"N_newcells="<<N_newcells<<endl;
    
    EG_VTKSP(vtkUnstructuredGrid,grid_tmp);
    allocateGrid(grid_tmp,N_cells+N_newcells,N_points+N_newpoints);
    makeCopyNoAlloc(grid, grid_tmp);
    EG_VTKDCC(vtkIntArray, cell_code_tmp, grid_tmp, "cell_code");
    
    vtkIdType nodeId=N_points;
    foreach(stencil_t S,StencilVector)
    {
      cout<<S<<endl;
      vtkIdType N_pts, *pts;
      vec3_t A,B;
      grid_tmp->GetPoint(S.p[1],A.data());
      grid_tmp->GetPoint(S.p[3],B.data());
      vec3_t M=0.5*(A+B);
      addPoint(grid_tmp,nodeId,M.data());
      
      vtkIdType pts_triangle[4][3];
      
      for(int i=0;i<4;i++)
      {
        pts_triangle[i][0]=S.p[i];
        pts_triangle[i][1]=S.p[(i+1)%4];
        pts_triangle[i][2]=nodeId;
      }
      
      int bc1=cell_code_tmp->GetValue(S.id_cell1);
      int bc2=cell_code_tmp->GetValue(S.id_cell2);
      
      grid_tmp->ReplaceCell(S.id_cell1 , 3, pts_triangle[0]);
      cell_code_tmp->SetValue(S.id_cell1, bc1);
      
      grid_tmp->ReplaceCell(S.id_cell2 , 3, pts_triangle[1]);
      cell_code_tmp->SetValue(S.id_cell2, bc2);
      
      vtkIdType newCellId;
      newCellId = grid_tmp->InsertNextCell(VTK_TRIANGLE,3,pts_triangle[2]);
      cell_code_tmp->SetValue(newCellId, bc2);
      newCellId = grid_tmp->InsertNextCell(VTK_TRIANGLE,3,pts_triangle[3]);
      cell_code_tmp->SetValue(newCellId, bc1);
      
      nodeId++;
    }
    
    makeCopy(grid_tmp,grid);
//     cout_grid(cout,grid,true,true,true,true);*/
    updateActors();
  }
  //////////////////////////////////////////////////////////////////////////////////////////////
  else if(ui.SmoothMethod->currentIndex()==9)//vtkWindowedSincPolyDataFilter smoothing
  {
	cout<<"HUHUHHHUHUHUHUH"<<endl;
    QSet<int> bcs;
    getSelectedItems(ui.listWidget, bcs);
    QVector<vtkIdType> cells;
    getSurfaceCells(bcs, cells, grid);
    EG_VTKSP(vtkPolyData, input);
    addToPolyData(cells, input, grid);

    EG_VTKSP(vtkSmoothPolyDataFilter, smooth);
    smooth->SetInput(input);


  int j, k;
  vtkIdType npts = 0;
  vtkIdType *pts = 0;

  vtkCellArray *inVerts, *inLines, *inPolys, *inStrips;

	cout<<"input->GetVerts()="<<input->GetVerts()<<endl;
	cout<<"input->GetVerts()->GetSize()="<<input->GetVerts()->GetSize()<<endl;
inVerts=input->GetVerts();
inVerts->InitTraversal();
cout<<"inVerts->GetSize()="<<inVerts->GetSize()<<endl;
cout<<"inVerts->GetNextCell(npts,pts)="<<inVerts->GetNextCell(npts,pts)<<endl;
cout<<"inVerts->GetNumberOfCells()="<<inVerts->GetNumberOfCells()<<endl;

cout<<"input->GetVerts()->GetSize()="<<input->GetVerts()->GetSize()<<endl;
cout<<"input->GetVerts()->GetNextCell(npts,pts)="<<input->GetVerts()->GetNextCell(npts,pts)<<endl;
cout<<"input->GetVerts()->GetNumberOfCells()="<<input->GetVerts()->GetNumberOfCells()<<endl;

/*inVerts->GetNumberOfTuples();
inVerts->GetNumberOfComponents();*/
  // check vertices first. Vertices are never smoothed_--------------
  for (inVerts=input->GetVerts(), inVerts->InitTraversal(); 
  inVerts->GetNextCell(npts,pts); )
    {
	cout<<"npts="<<npts<<endl;
    for (j=0; j<npts; j++)
      {
	cout<<"pts["<<j<<"]="<<pts[j]<<endl;
//       Verts[pts[j]].type = VTK_FIXED_VERTEX;
      }
    }

/*    EG_VTKSP(vtkWindowedSincPolyDataFilter, smooth);
    
    cout_vtkWindowedSincPolyDataFilter(smooth);
    
    smooth->SetInput(pdata);
    
    smooth->SetNumberOfIterations(ui.spinBox_NumberOfIterations->value());
    smooth->SetPassBand(ui.doubleSpinBox_PassBand->value());
    smooth->SetFeatureEdgeSmoothing(ui.checkBox_FeatureEdgeSmoothing->checkState());
    smooth->SetFeatureAngle(ui.doubleSpinBox_FeatureAngle->value());
    smooth->SetEdgeAngle(ui.doubleSpinBox_EdgeAngle->value());
    smooth->SetBoundarySmoothing(ui.checkBox_BoundarySmoothing->checkState());
    smooth->SetGenerateErrorScalars(ui.checkBox_GenerateErrorScalars->checkState());
    smooth->SetGenerateErrorVectors(ui.checkBox_GenerateErrorVectors->checkState());
    
    cout_vtkWindowedSincPolyDataFilter(smooth);
    
    smooth->Update();
    EG_VTKDCN(vtkLongArray_t, node_index, pdata, "node_index");
    for (vtkIdType i = 0; i < smooth->GetOutput()->GetNumberOfPoints(); ++i) {
      vec3_t x;
      smooth->GetOutput()->GetPoints()->GetPoint(i, x.data());
      vtkIdType nodeId = node_index->GetValue(i);
      grid->GetPoints()->SetPoint(nodeId, x.data());
    };*/
    updateActors();
  }
  //////////////////////////////////////////////////////////////////////////////////////////////
  if(ui.SmoothMethod->currentIndex()==10)//vtkSmoothPolyDataFilter2 smoothing
  {
    QSet<int> bcs;
    getSelectedItems(ui.listWidget, bcs);
    QVector<vtkIdType> cells;
    getSurfaceCells(bcs, cells, grid);
    EG_VTKSP(vtkPolyData, pdata);
    addToPolyData(cells, pdata, grid);
    EG_VTKSP(vtkSmoothPolyDataFilter2, smooth);
    
//     cout_vtkSmoothPolyDataFilter2(smooth);
    
    smooth->SetInput(pdata);
    
    smooth->SetConvergence (ui.doubleSpinBox_Convergence->value());
    smooth->SetNumberOfIterations (ui.spinBox_NumberOfIterations->value());
    smooth->SetRelaxationFactor (ui.lineEdit_RelaxationFactor->text().toDouble());
    smooth->SetFeatureEdgeSmoothing (ui.checkBox_FeatureEdgeSmoothing->checkState());
    smooth->SetFeatureAngle (ui.doubleSpinBox_FeatureAngle->value());
    smooth->SetEdgeAngle (ui.doubleSpinBox_EdgeAngle->value());
    smooth->SetBoundarySmoothing (ui.checkBox_BoundarySmoothing->checkState());
    smooth->SetGenerateErrorScalars (ui.checkBox_GenerateErrorScalars->checkState());
    smooth->SetGenerateErrorVectors (ui.checkBox_GenerateErrorVectors->checkState());
    
    QSet<int> bcs_Source;
    getSelectedItems(ui.listWidget,bcs_Source);
    QVector<vtkIdType> cells_Source;
    getSurfaceCells(bcs_Source, cells_Source, grid);
    EG_VTKSP(vtkPolyData, pdata_Source);
    addToPolyData(cells_Source, pdata_Source, grid);
    smooth->SetSource (pdata_Source);
    
//     cout_vtkSmoothPolyDataFilter2(smooth);

    smooth->DebugOn();
    smooth->Update();

    cout<<"==============="<<endl;
    int N_components,N_tuples,N_arrays;
    N_arrays=smooth->GetOutput()->GetPointData()->GetNumberOfArrays();
    cout<<"N_arrays="<<N_arrays<<endl;

    cout<<"ErrorScalars:"<<endl;
    vtkFloatArray *newScalars = vtkFloatArray::New();
    newScalars=(vtkFloatArray *)smooth->GetOutput()->GetPointData()->GetArray(1);
    N_components=newScalars->GetNumberOfComponents();
    N_tuples=newScalars->GetNumberOfTuples();
    cout<<"N_components="<<N_components<<endl;
    cout<<"N_tuples="<<N_tuples<<endl;
    for (int i=0; i<N_tuples; i++)
    {
      double dist=newScalars->GetComponent(i-1,1);//strange, but works. O.o
      cout<<"dist="<<dist<<endl;
    }
    
    cout<<"ErrorVectors:"<<endl;
    N_components=smooth->GetOutput()->GetPointData()->GetVectors()->GetNumberOfComponents();
    N_tuples=smooth->GetOutput()->GetPointData()->GetVectors()->GetNumberOfTuples();
    cout<<"N_components="<<N_components<<endl;
    cout<<"N_tuples="<<N_tuples<<endl;
    for(vtkIdType i=0;i<N_tuples;i++)
    {
      double tuple[4];
      smooth->GetOutput()->GetPointData()->GetTuple(i,tuple);
      cout<<"tuple["<<tuple[0]<<"]=("<<tuple[1]<<","<<tuple[2]<<","<<tuple[3]<<")"<<endl;
    }
    cout<<"==============="<<endl;

    updateActors();
  }
  //////////////////////////////////////////////////////////////////////////////////////////////
  else
  {
    cout<<"UNKNOWN METHOD"<<endl;
  }
};
