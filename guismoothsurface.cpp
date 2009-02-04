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

/* Here is how we we get QTextStreams that look like iostreams */

// QTextStream cin(stdin, QIODevice::ReadOnly);
// QTextStream cout(stdout, QIODevice::WriteOnly);
// QTextStream cerr(stderr, QIODevice::WriteOnly);

void GuiSmoothSurface::before()
{
  populateBoundaryCodes(ui.listWidget, grid);
  ui.SmoothMethod-> setCurrentIndex(1);
};

void GuiSmoothSurface::operate()
{
  //can't use switch case because dynamic variables seem to be forbidden inside case statements
  //////////////////////////////////////////////////////////////////////////////////////////////
  if(ui.SmoothMethod->currentIndex()==0)
  {
    cout<<"METHOD 0"<<endl;
    QSet<int> bcs;
    getSelectedItems(ui.listWidget, bcs);
    QVector<vtkIdType> cells;
    getSurfaceCells(bcs, cells, grid);
    EG_VTKSP(vtkPolyData, pdata);
    addToPolyData(cells, pdata, grid);
    EG_VTKSP(vtkSmoothPolyDataFilter, smooth);
    smooth->SetInput(pdata);
    smooth->FeatureEdgeSmoothingOn();
    smooth->SetFeatureAngle(ui.doubleSpinBox->value());
    smooth->SetNumberOfIterations(ui.spinBox->value());
    smooth->Update();
    EG_VTKDCN(vtkLongArray_t, node_index, pdata, "node_index");
    for (vtkIdType i = 0; i < smooth->GetOutput()->GetNumberOfPoints(); ++i) {
      vec3_t x;
//       smooth->GetOutput() => returns vtkPolyData
      //smooth->GetOutput()->GetPoints() returns vtkPoints
      smooth->GetOutput()->GetPoints()->GetPoint(i, x.data());
      vtkIdType nodeId = node_index->GetValue(i);
      grid->GetPoints()->SetPoint(nodeId, x.data());
    };
    updateActors();
  }
  //////////////////////////////////////////////////////////////////////////////////////////////
  else if(ui.SmoothMethod->currentIndex()==1)
  {
    ofstream myfile;
    myfile.open ("/data1/home/mtaverne/debug.log");
    
    myfile<<"METHOD 1"<<endl;
    QSet<int> bcs;
    getSelectedItems(ui.listWidget, bcs);
    QVector<vtkIdType> cells;
    getSurfaceCells(bcs, cells, grid);
    
    foreach (int value, cells)
      myfile << value << endl;
    
    vtkIdType new_pts1[3], new_pts2[3];
    new_pts1[0] = 3;
    new_pts1[1] = 0;
    new_pts1[2] = 1;
    new_pts2[0] = 3;
    new_pts2[1] = 1;
    new_pts2[2] = 2;
    
    grid->ReplaceCell(0 , 3, new_pts1);
    grid->ReplaceCell(1 , 3, new_pts2);

    //Allocate temp grid
    EG_VTKSP(vtkUnstructuredGrid,grid_tmp);
    allocateGrid(grid_tmp,14,9);
    myfile << "ALLOCATED" ;
    myfile<<"grid_tmp:"<<endl;cout_grid(myfile, grid_tmp);
    makeCopy_noAlloc(grid, grid_tmp);
    myfile<<"grid_tmp:"<<endl;cout_grid(myfile, grid_tmp);
    
    vec3_t newpoint(5,5,5);
    myfile << grid <<endl;
    myfile<<"grid_tmp:"<<endl;cout_grid(myfile, grid_tmp);
    myfile<<"grid:"<<endl;cout_grid(myfile, grid);
    
    grid_tmp->GetPoints()->SetPoint(8, newpoint.data());
    myfile << "new point added!!!!" << endl;
    cout_grid(myfile,grid_tmp,true,true,true,true);
/*    vtkIdList* ptIds;
    ptIds->Allocate(3);
    ptIds[0]=0;
    ptIds[1]=1;
    ptIds[2]=9;*/
//     grid_tmp->InsertNextCell(VTK_TRIANGLE,ptIds);
    vtkIdType npts=3;
    vtkIdType pts[3];
    pts[0]=0;
    pts[1]=1;
    pts[2]=8;
    vtkIdType newCellId = grid_tmp->InsertNextCell(VTK_TRIANGLE,npts,pts);
    myfile << "new cell added!!!! : " << newCellId << endl;
    cout_grid(myfile,grid_tmp,true,true,true,true);
    
    // boundary conditions
    EG_VTKDCC(vtkIntArray, cell_code, grid_tmp, "cell_code");
    cell_code->SetValue(newCellId, ui.lineEdit_BoundaryCode->text().toInt());
    
    //Copy temp grid to grid
    makeCopy(grid_tmp, grid);
    myfile<<"COPIED!!!!"<<endl;
    myfile<<"grid_tmp:"<<endl;cout_grid(myfile, grid_tmp);
    myfile<<"grid:"<<endl;cout_grid(myfile, grid,true,true);
    myfile<<"grid:"<<endl;cout_grid(myfile, grid,true,true,true,true);
    cout<<"OUTPUT: grid:"<<endl;cout_grid(cout, grid,true,true,true,true);
    
    myfile.close();
    updateActors();
  }
  //////////////////////////////////////////////////////////////////////////////////////////////
  else if(ui.SmoothMethod->currentIndex()==2)
  {
    cout<<"METHOD 2"<<endl;
    EG_VTKSP(vtkEgNormalExtrusion, extr);
    QVector<double> y;
    cout <<"y.size()="<<y.size()<<endl;
      y.resize(3 + 1);
    cout <<"y.size()="<<y.size()<<endl;
    double h = 1;
      double f = 1;
      y[0] = 0.0;
      for (int i = 1; i < y.size(); ++i) {
        y[i] = y[i-1] + h;
        h *= f;
      };
    cout << "y=" ;
    for(int i=0;i<y.size();i++) cout << y[i] << ",";
    cout << endl;
    extr->SetLayers(y);
    
    QSet<int> bcs;
    getSelectedItems(ui.listWidget, bcs);
    extr->SetBoundaryCodes(&bcs);
    EG_VTKSP(vtkUnstructuredGrid,ug);
    makeCopy(grid, ug);
    extr->SetInput(ug);
    cout << "Pre-update: extr->GetOutput()->GetNumberOfPoints()=" << extr->GetOutput()->GetNumberOfPoints() << endl;
    extr->Update();
    cout << "Post-update: extr->GetOutput()->GetNumberOfPoints()=" << extr->GetOutput()->GetNumberOfPoints() << endl;
    
    extr->GetOutput();
    
    for (vtkIdType i = 0; i < extr->GetOutput()->GetNumberOfPoints(); ++i) {
      vec3_t x;
      extr->GetOutput()->GetPoint(i, x.data());
//       vtkIdType nodeId = node_index->GetValue(i);
      cout << "Vertex " << i << " = " << x << endl;
      extr->GetOutput()->GetPoints()->SetPoint(i, x.data());
    };
    
    cout << "Pre-copy: grid->GetNumberOfPoints()=" << grid->GetNumberOfPoints() << endl;
    makeCopy(extr->GetOutput(), grid);
    cout << "Post-copy: grid->GetNumberOfPoints()=" << grid->GetNumberOfPoints() << endl;
    
  }
  //////////////////////////////////////////////////////////////////////////////////////////////
  else if(ui.SmoothMethod->currentIndex()==3)
  {
    
    ofstream myfile;
    myfile.open ("/data1/home/mtaverne/debug.log");
    cout_grid(myfile,grid,true,true,true,true);
    myfile.close();
    
    cout_grid(cout,grid,true,true,true,true);
    
  }
  //////////////////////////////////////////////////////////////////////////////////////////////
  else if(ui.SmoothMethod->currentIndex()==4)
  {
    ofstream myfile;
    myfile.open ("/data1/home/mtaverne/debug.log");
    
    EG_VTKSP(vtkUnstructuredGrid,grid_tmp);
    allocateGrid(grid_tmp,14,9);
    myfile<<"grid_tmp:"<<endl;cout_grid(myfile, grid_tmp);
    
    myfile.close();
    updateActors();
  }
  //////////////////////////////////////////////////////////////////////////////////////////////
  else
  {
    cout<<"UNKNOWN METHOD"<<endl;
  }
};
