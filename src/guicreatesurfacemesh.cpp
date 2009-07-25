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
#include "guicreatesurfacemesh.h"

#include "swaptriangles.h"
#include "surfacemesher.h"
#include "vertexdelegate.h"
#include "settingssheet.h"
#include "laplacesmoother.h"
#include "guimainwindow.h"
#include "containertricks.h"
#include "updatedesiredmeshdensity.h"

#include <vtkSmoothPolyDataFilter.h>
#include <vtkWindowedSincPolyDataFilter.h>
#include <vtkLongArray.h>
#include <vtkEgNormalExtrusion.h>
#include <vtkIdList.h>
#include <vtkCell.h>
#include <vtkCellLocator.h>
#include <vtkFloatArray.h>
#include <vtkCellArray.h>

#include <QtGui>
#include <QTextStream>

#include <iostream>
#include <fstream>
#include <cmath>

#include <stdio.h>

//////////////////////////////////////////////

GuiCreateSurfaceMesh::GuiCreateSurfaceMesh()
{
  EG_TYPENAME;
  
  setQuickSave(true);

  m_tableWidget = new SettingsSheet();
  ui.verticalLayout_SettingsSheet->addWidget(m_tableWidget);
  
  populateBoundaryCodes(ui.listWidget);
  ui.lineEditMaximalEdgeLength->setText("1000");

  //Load settings
  readSettings();
  
  int Nrow,Ncol;
  Nrow=m_tableWidget->rowCount();
  Ncol=m_tableWidget->columnCount();
  
  QList<QString> list;
  list
    <<"VTK_SIMPLE_VERTEX"
    <<"VTK_FIXED_VERTEX"
    <<"VTK_FEATURE_EDGE_VERTEX"
    <<"VTK_BOUNDARY_EDGE_VERTEX"
    <<"any";
  
  QList<QString> list2;
  list2
    << "yes"
    <<"no"
    <<"any";
  
  Nbc=ui.listWidget-> count ();
  m_tableWidget->setColumnCount(Nbc+3);
  VertexDelegate* item_delegate = new VertexDelegate(Nbc, list);
  m_tableWidget->setItemDelegate(item_delegate);
  
  QStringList L;
  for(int i = 0; i < Nbc; ++i) {
    L << ui.listWidget->item(i)->text().split(":")[0];
  }
  L<<"Vertex Type";
  L<<"Nodelist";
  L<<"Mesh Density";
  m_tableWidget->setHorizontalHeaderLabels(L);
  m_tableWidget->resizeColumnsToContents();
  
  
  current_filename= GuiMainWindow::pointer()->getFilename();
  qDebug()<<"current_filename="<<current_filename;
  qDebug()<<"Loading settings from "+current_filename+".sp...";
  
  if (!m_tableWidget->readFile(0)) {
    cout<<"Loading settingssheet failed"<<endl;
  }
  
  connect(ui.pushButton_AddSet, SIGNAL(clicked()), this, SLOT(AddSet()));
  connect(ui.pushButton_RemoveSet, SIGNAL(clicked()), this, SLOT(RemoveSet()));
  connect(ui.pushButton_TestSet, SIGNAL(clicked()), this, SLOT(TestSet()));
  connect(ui.pushButton_SelectAll_BC, SIGNAL(clicked()), this, SLOT(SelectAll_BC()));
  connect(ui.pushButton_ClearAll_BC, SIGNAL(clicked()), this, SLOT(ClearAll_BC()));
  connect(ui.pushButtonSave, SIGNAL(clicked()), this, SLOT(writeSettings()));
  connect(ui.pushButtonSave, SIGNAL(clicked()), m_tableWidget, SLOT(writeFile()));

}

///////////////////////////////////////////

int GuiCreateSurfaceMesh::readSettings()
{
  QString buffer = GuiMainWindow::pointer()->getXmlSection("engrid/surface/settings");
  QTextStream in(&buffer, QIODevice::ReadOnly);
  QString str;
  in >> str;
  ui.lineEditMaximalEdgeLength->setText(str);
  double nodes_per_quarter_circle;
  in >> nodes_per_quarter_circle;
  ui.doubleSpinBoxCurvature->setValue(nodes_per_quarter_circle);
  int num_bcs;
  in >> num_bcs;
  if (num_bcs == ui.listWidget->count()) {
    int check_state;
    for (int i = 0; i < ui.listWidget->count(); ++i) {
      in >> check_state;
      if (check_state == 1) {
        ui.listWidget->item(i)->setCheckState(Qt::Checked);
      } else {
        ui.listWidget->item(i)->setCheckState(Qt::Unchecked);
      }
    }
  }
  return(0);
}

int GuiCreateSurfaceMesh::writeSettings()
{
  QString buffer = "";
  {
    QTextStream out(&buffer, QIODevice::WriteOnly);
    out << "\n";
    out << ui.lineEditMaximalEdgeLength->text() << "\n";
    out << ui.doubleSpinBoxCurvature->value() << "\n";
    out << ui.listWidget->count() << "\n";
    for (int i = 0; i < ui.listWidget->count(); ++i) {
      if (ui.listWidget->item(i)->checkState() == Qt::Checked) {
        out << "1 \n";
      } else {
        out << "0 \n";
      }
    }
  }
  GuiMainWindow::pointer()->setXmlSection("engrid/surface/settings", buffer);
  return(0);
}

///////////////////////////////////////////

void GuiCreateSurfaceMesh::SelectAll_BC()
{
  for (int i = 0; i < ui.listWidget->count(); ++i) {
    ui.listWidget->item(i)->setCheckState(Qt::Checked);
  }
}

void GuiCreateSurfaceMesh::ClearAll_BC()
{
  for (int i = 0; i < ui.listWidget->count(); ++i) {
    ui.listWidget->item(i)->setCheckState(Qt::Unchecked);
  }
}

void GuiCreateSurfaceMesh::TestSet()
{
  cout<<"Testing set"<<endl;
  getSet();
}

//This is where we get the user defined mesh densities
QVector <VertexMeshDensity> GuiCreateSurfaceMesh::getSet()
{
  cout<<"Getting set"<<endl;
  QVector <VertexMeshDensity> VMDvector;
  
  cout<<"VMDvector:"<<VMDvector<<endl;
  
  int N_VMD=m_tableWidget->rowCount();
  VMDvector.resize(N_VMD);
  cout<<"VMDvector.size()="<<VMDvector.size()<<endl;
  for(int i=0;i<N_VMD;i++)
  {
    for(int j=0;j<Nbc;j++)
    {
      int bc = m_tableWidget->horizontalHeaderItem(j)->text().toInt();
      int state = CheckState2int( m_tableWidget->item(i,j)->checkState() );
      VMDvector[i].BCmap[bc]=state;
    }
    VMDvector[i].type=Str2VertexType(m_tableWidget->item(i,Nbc)->text());
    VMDvector[i].setNodes(m_tableWidget->item(i,Nbc+1)->text());
    VMDvector[i].density=m_tableWidget->item(i,Nbc+2)->text().toDouble();
  }
  cout<<"VMDvector:"<<VMDvector<<endl;
  return(VMDvector);
}

void GuiCreateSurfaceMesh::AddSet()
{
  cout<<"Adding set"<<endl;
  int row=m_tableWidget->rowCount();
  m_tableWidget->insertRow(row);
  
  int Nbc=ui.listWidget->count();
  for(int i=0;i<Nbc;i++)
  {
    TriStateTableWidgetItem* newBC = new TriStateTableWidgetItem();
    newBC->setFlags(Qt::ItemIsTristate | Qt::ItemIsEnabled | Qt::ItemIsUserCheckable);
    m_tableWidget->setItem(row, i, newBC);
  }
  QTableWidgetItem* item1 = new QTableWidgetItem("any");
  m_tableWidget->setItem(row, Nbc, item1);
  QTableWidgetItem* item2 = new QTableWidgetItem("");
  m_tableWidget->setItem(row, Nbc+1, item2);
  QTableWidgetItem* item3 = new QTableWidgetItem("-1");
  m_tableWidget->setItem(row, Nbc+2, item3);
  m_tableWidget->resizeColumnsToContents();
}

void GuiCreateSurfaceMesh::RemoveSet()
{
  cout<<"Removing set"<<endl;
  m_tableWidget->removeRow(m_tableWidget->currentRow());
  m_tableWidget->resizeColumnsToContents();
}

int GuiCreateSurfaceMesh::DisplayErrorScalars(vtkPolyDataAlgorithm* algo)
{
  return(0);
  cout<<"==============="<<endl;
  cout<<"ErrorScalars:"<<endl;
  int N1,N2;
  double dist;
  N1=algo->GetOutput()->GetPointData()->GetNumberOfArrays();
//   cout<<"nb of arrays="<<N1<<endl;
  algo->GetOutput();//vtkPolyData
//   cout<<algo->GetOutput()->GetPointData()<<endl;//vtkPointData*
  vtkFloatArray *newScalars = vtkFloatArray::New();
  newScalars=(vtkFloatArray *)algo->GetOutput()->GetPointData()->GetArray(1);
  N1=newScalars->GetNumberOfComponents();
  N2=newScalars->GetNumberOfTuples();
  cout<<"Number of components=N1="<<N1<<endl;
  cout<<"Number of tuples=N2="<<N2<<endl;
  for (int i=0; i<N2; i++)
  {
    dist=newScalars->GetComponent(i-1,1);//strange, but works. O.o
    cout<<"dist["<<i<<"]="<<dist<<endl;
  }
  
  cout<<"==============="<<endl;
  newScalars->Delete();
  return(0);
}

void GuiCreateSurfaceMesh::operate()
{
  writeSettings();

  m_tableWidget->writeFile();
  return;


  QSet<int> bcs;
  getSelectedItems(ui.listWidget, bcs);

  QVector <VertexMeshDensity> VMDvector = getSet();

  SurfaceMesher surfacemesher;
  surfacemesher.setGrid(grid);
  surfacemesher.setBoundaryCodes(bcs);
  surfacemesher.setVertexMeshDensityVector(VMDvector);
  surfacemesher.setMaxEdgeLength(ui.lineEditMaximalEdgeLength->text().toDouble());
  surfacemesher.setNodesPerQuarterCircle(ui.doubleSpinBoxCurvature->value());

  surfacemesher();

  grid->Modified();
}
